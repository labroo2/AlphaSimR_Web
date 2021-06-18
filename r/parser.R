library(rjson)
library(AlphaSimR)



# import the breeding scheme
scheme <- fromJSON(file = "json_MRL_20211606.json")



# make the founder population
nInd <- scheme$Nodes[which(sapply(scheme$Nodes, `[[`, "Node Type") == "Crossing Block")][[1]]["Number of Individuals"][[1]] # from crossing block node
nChr <- scheme$Genome$`Number of Chromosomes`
segSites <- max(sapply(scheme$`Trait Info`, `[[`, "QTL")) + scheme$Genome$nSnp # nQTL + nSNPs
inbred <- as.logical(scheme$Genome$Inbred)
Ne <- scheme$Genome$Ne
ploidy <- scheme$Genome$Ploidy

founderPop <- runMacs2(nInd = nInd,            
                       nChr = nChr,
                       segSites = segSites,
                       inbred = inbred,
                       Ne = Ne,
                       ploidy = ploidy)



# make the simulation parameters
SP <- SimParam$new(founderPop)

nTraits <- length(scheme$`Trait Info`)

# add the traits
for(i in nTraits){
  
  # addTraitA
  if(scheme$`Trait Info`[[i]]$`Trait Type` == "A"){
    SP$addTraitA(nQtlPerChr = scheme$`Trait Info`[[i]]$QTL, 
                 mean = scheme$`Trait Info`[[i]]$Mean,
                 var = scheme$`Trait Info`[[i]]$`Genetic Variance`
                 )
  }
  
  # addTraitAD
  if(scheme$`Trait Info`[[i]]$`Trait Type` == "AD"){
    SP$addTraitAD(nQtlPerChr = scheme$`Trait Info`[[i]]$QTL,
                  mean = scheme$`Trait Info`[[i]]$Mean,
                  var = scheme$`Trait Info`[[i]]$`Genetic Variance`,
                  meanDD = scheme$`Trait Info`[[i]]$`Mean Dominance Degree`,
                  varDD = scheme$`Trait Info`[[i]]$`Variance of Dominance Degrees`
                  )
  }
    
    # addTraitAG
    if(scheme$`Trait Info`[[i]]$`Trait Type` == "AG"){
      SP$addTraitAG(nQtlPerChr = scheme$`Trait Info`[[i]]$QTL,
                    mean = scheme$`Trait Info`[[i]]$Mean,
                    var = scheme$`Trait Info`[[i]]$`Genetic Variance`,
                    varGxE = scheme$`Trait Info`[[i]]$`GxE Variance`
                    # any need for the varEnv argument?
                    )
    }
    
    # addTraitADG
    if(scheme$`Trait Info`[[i]]$`Trait Type` == "ADG"){
      SP$addTraitADG(nQtlPerChr = scheme$`Trait Info`[[i]]$QTL,
                     mean = scheme$`Trait Info`[[i]]$Mean,
                     var = scheme$`Trait Info`[[i]]$`Genetic Variance`,
                     # any need for the varEnv argument?
                     varGxE = scheme$`Trait Info`[[i]]$`GxE Variance`,
                     meanDD = scheme$`Trait Info`[[i]]$`Mean Dominance Degree`,
                     varDD = scheme$`Trait Info`[[i]]$`Variance of Dominance Degrees`
                     )
    }
    
    # addTraitAE
    if(scheme$`Trait Info`[[i]]$`Trait Type` == "AE"){
      SP$addTraitAE(nQtlPerChr = scheme$`Trait Info`[[i]]$QTL,
                    mean = scheme$`Trait Info`[[i]]$Mean,
                    var = scheme$`Trait Info`[[i]]$`Genetic Variance`,
                    relAA = scheme$`Trait Info`[[i]]$`Relative Epistasis Variance`
                    )
    }
    
    # addTraitADE
      if(scheme$`Trait Info`[[i]]$`Trait Type` == "ADE"){
        SP$addTraitADE(nQtlPerChr = scheme$`Trait Info`[[i]]$QTL,
                       mean = scheme$`Trait Info`[[i]]$Mean,
                       var = scheme$`Trait Info`[[i]]$`Genetic Variance`,
                       meanDD = scheme$`Trait Info`[[i]]$`Mean Dominance Degree`,
                       varDD = scheme$`Trait Info`[[i]]$`Variance of Dominance Degrees`,
                       relAA= scheme$`Trait Info`[[i]]$`Relative Epistasis Variance`
                       )
      }
    
    # addTraitAEG
    if(scheme$`Trait Info`[[i]]$`Trait Type` == "AEG"){
      SP$addTraitAEG(nQtlPerChr = scheme$`Trait Info`[[i]]$QTL,
                     mean = scheme$`Trait Info`[[i]]$Mean,
                     var = scheme$`Trait Info`[[i]]$`Genetic Variance`,
                     relAA= scheme$`Trait Info`[[i]]$`Relative Epistasis Variance`,
                     varGxE = scheme$`Trait Info`[[i]]$`GxE Variance`
                     # is varEnv option needed?
                     )
    }
    
    # addTraitADEG
    if(scheme$`Trait Info`[[i]]$`Trait Type` == "ADEG"){
      SP$addTraitADEG(nQtlPerChr = scheme$`Trait Info`[[i]]$QTL,
                      mean = scheme$`Trait Info`[[i]]$Mean,
                      var = scheme$`Trait Info`[[i]]$`Genetic Variance`,
                      # is varEnv option needed?
                      varGxE = scheme$`Trait Info`[[i]]$`GxE Variance`,
                      meanDD = scheme$`Trait Info`[[i]]$`Mean Dominance Degree`,
                      varDD = scheme$`Trait Info`[[i]]$`Variance of Dominance Degrees`,
                      relAA= scheme$`Trait Info`[[i]]$`Relative Epistasis Variance`
                      )
    }
} # end trait addition loop

SP$addSnpChip(nSnpPerChr = scheme$Genome$nSnp) # add this to the JSON file if 




# make the breeding scheme

startPop <- newPop(founderPop) #pull the initial individuals
nCrosses <- scheme$Edges[which(sapply(scheme$Edges, `[[`, "Edge Type") == "Cross")][[1]]$nCrosses
nProg <- scheme$Edges[which(sapply(scheme$Edges, `[[`, "Edge Type") == "Cross")][[1]]$nProgenyPerCross

# fill breeding pipeline with unique individuals from initial parents
crossingBlock <- randCross(startPop,
                           nCrosses = nCrosses,
                           nProgeny = nProg)


# just don't rename object/modify same object

# store all populations forever

  
  



#traitMaker <- matrix(nrow = nTraits, ncol = 4); colnames(traitMaker) <- c("A", "D", "E", "G")

#for(i in nTraits){
  
#  # if GxE is not zero, add a G
#  if(scheme$`Trait Info`[[i]]$`GxE Variance` != 0){
#    traitMaker[i, "G"] <- "G"
#  }
  
#  # if there is dominance, add a D
#  if(scheme$`Trait Info`[[i]]$Dominance == "TRUE"){
#    traitMaker[i, "D"] <- "D"
#  }
#  
#}









lapply(scheme$Nodes,
       function(x) grep("Crossing Block", x))
mapply(scheme$Nodes,
       function(x) grep("Crossing Block", x))
