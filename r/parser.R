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

# add the traits
nTraits <- length(scheme$`Trait Info`)
traitTypes <- c()

# define the trait types
for(i in 1:nTraits){
  
  type <- "A"
 
  # check for D
  if(as.logical(scheme$`Trait Info`[[i]]$Dominance) == TRUE){
    type <- paste(type, "D", sep = "")
  }
  
  # check for E
  if(as.logical(scheme$`Trait Info`[[i]]$`Additive x Additive Epistasis`) == TRUE){
    type <- paste(type, "E", sep = "")
  }
  
  # check for G
  if(as.logical(scheme$`Trait Info`[[i]]$GxE) == TRUE){
    type <- paste(type, "G", sep = "")
  }

  # paste to make the traitTypes
  traitTypes[i] <- type
    
} # end trait definition



# add the traits
for(i in nTraits){
  
  # addTraitA
  if(traitTypes[i] == "A"){
    SP$addTraitA(nQtlPerChr = scheme$`Trait Info`[[i]]$QTL, 
                 mean = scheme$`Trait Info`[[i]]$Mean,
                 var = scheme$`Trait Info`[[i]]$`Genetic Variance`
    )
  }
  
  # addTraitAD
  if(traitTypes[i] == "AD"){
    SP$addTraitAD(nQtlPerChr = scheme$`Trait Info`[[i]]$QTL,
                  mean = scheme$`Trait Info`[[i]]$Mean,
                  var = scheme$`Trait Info`[[i]]$`Genetic Variance`,
                  meanDD = scheme$`Trait Info`[[i]]$`Mean Dominance Degree`,
                  varDD = scheme$`Trait Info`[[i]]$`Variance of Dominance Degrees`
    )
  }
  
  # addTraitAG
  if(traitTypes[i] == "AG"){
    SP$addTraitAG(nQtlPerChr = scheme$`Trait Info`[[i]]$QTL,
                  mean = scheme$`Trait Info`[[i]]$Mean,
                  var = scheme$`Trait Info`[[i]]$`Genetic Variance`,
                  varGxE = scheme$`Trait Info`[[i]]$`GxE Variance`
                  # any need for the varEnv argument?
    )
  }
  
  # addTraitADG
  if(traitTypes[i] == "ADG"){
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
  if(traitTypes[i] == "AE"){
    SP$addTraitAE(nQtlPerChr = scheme$`Trait Info`[[i]]$QTL,
                  mean = scheme$`Trait Info`[[i]]$Mean,
                  var = scheme$`Trait Info`[[i]]$`Genetic Variance`,
                  relAA = scheme$`Trait Info`[[i]]$`Relative Epistasis Variance`
    )
  }
  
  # addTraitADE
  if(traitTypes[i] == "ADE"){
    SP$addTraitADE(nQtlPerChr = scheme$`Trait Info`[[i]]$QTL,
                   mean = scheme$`Trait Info`[[i]]$Mean,
                   var = scheme$`Trait Info`[[i]]$`Genetic Variance`,
                   meanDD = scheme$`Trait Info`[[i]]$`Mean Dominance Degree`,
                   varDD = scheme$`Trait Info`[[i]]$`Variance of Dominance Degrees`,
                   relAA= scheme$`Trait Info`[[i]]$`Relative Epistasis Variance`
    )
  }
  
  # addTraitAEG
  if(traitTypes[i] == "AEG"){
    SP$addTraitAEG(nQtlPerChr = scheme$`Trait Info`[[i]]$QTL,
                   mean = scheme$`Trait Info`[[i]]$Mean,
                   var = scheme$`Trait Info`[[i]]$`Genetic Variance`,
                   relAA= scheme$`Trait Info`[[i]]$`Relative Epistasis Variance`,
                   varGxE = scheme$`Trait Info`[[i]]$`GxE Variance`
                   # is varEnv option needed?
    )
  }
  
  # addTraitADEG
  if(traitTypes[i] == "ADEG"){
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

  
  






