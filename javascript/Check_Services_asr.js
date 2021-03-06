// this function will always be called when the body is focused on,
// it will check all possibly inconsistencies if anything of the following components are changed:
function checkEverything(id){
	console.log(id);
	

if(id == "General_Info"){  
	
	data_Vue.warnings1 = [];

	gen_warn_text = "Please specify a Project Name.";
	if(!data_Vue.geninfo['Project Name'] & data_Vue.warnings1.indexOf(gen_warn_text) == -1){
		data_Vue.warnings1.push(gen_warn_text);
	}

	var pos1 = data_Vue.geninfo['Project Name'].search(",")
	gen_warn_text = "Project name is not allowed to contain ',' ";
	if(pos1 >0 ){
		data_Vue.warnings1.push(gen_warn_text);
	}
	gen_warn_text = "Please specify a Species.";
	if(!data_Vue.geninfo['Species'] & data_Vue.warnings1.indexOf(gen_warn_text) == -1){
		data_Vue.warnings1.push(gen_warn_text);
	}
					
	gen_warn_text = "Please specify Time Unit.";
	if(!data_Vue.geninfo['Time Unit'] & data_Vue.warnings1.indexOf(gen_warn_text) == -1){
		data_Vue.warnings1.push(gen_warn_text);
	}			
	
	
	gen_warn_text = "Ensembl Dataset not available for Species 'Other'";
	if(data_Vue.genetic_data =="Ens" & data_Vue.geninfo['Species'] =="Other"){
		data_Vue.warnings1.push(gen_warn_text);
	} else{
		gen_warn_textEns = "Please specify an Ensembl Dataset.";
		if ((data_Vue.genetic_data == "Ens" & !data_Vue.geninfo['Ensembl Dataset'] & data_Vue.warnings1.indexOf(gen_warn_textEns) == -1)){
			data_Vue.warnings1.push(gen_warn_textEns);
		}
	}
	


	gen_warn_textSNP = "Please specify Max.Number of SNPs. SNPs are in positive integer";
	checkPosSNP = isPositiveInt(data_Vue.geninfo['Max Number of SNPs']);
	if (data_Vue.genetic_data == "Ens" & (isNaN(data_Vue.geninfo['Max Number of SNPs']) || data_Vue.geninfo['Max Number of SNPs'] < 0  || (data_Vue.geninfo['Max Number of SNPs'] != "" & checkPosSNP == false )) & data_Vue.warnings1.indexOf(gen_warn_textSNP) == -1){
		data_Vue.warnings1.push(gen_warn_textSNP);
	}	
		
	gen_warn_text = "Please upload Own Map.";
	if ((data_Vue.genetic_data == "Own" & !data_Vue.geninfo['Own Map Path']) & data_Vue.warnings1.indexOf(gen_warn_text) == -1){
		data_Vue.warnings1.push(gen_warn_text);
	}
		
	gen_warn_text = "Please enter number of Chromosomes. Chromosomes are in positive integer.";			
	checkCus = isPositiveInteger(data_Vue.geninfo['Number of Chromosomes']) ? "true" : "false";
	checkNumeric = isNumeric(data_Vue.geninfo['Number of Chromosomes']);
	if ((data_Vue.genetic_data == "Cus" & !data_Vue.geninfo['Number of Chromosomes'] & checkCus == "false") & data_Vue.warnings1.indexOf(gen_warn_text) == -1){
		data_Vue.warnings1.push(gen_warn_text);
	}
	else if ((data_Vue.genetic_data == "Cus" & (data_Vue.geninfo['Number of Chromosomes'] != "" & (checkNumeric == "false") || checkCus == "false")) & data_Vue.warnings1.indexOf(gen_warn_text) == -1){
		data_Vue.warnings1.push(gen_warn_text);
	}



	thisChromoLen = data_Vue.geninfo['Chromosomes of Equal Length'];
		thisChromoNum = data_Vue.geninfo['Chromosomes Info'].length;

		if (thisChromoNum > 0 ) {
			
		gen_warn_text = "Simulation contains more than 1.000.000 SNPs! Check if that is needed!";

		if(data_Vue.geninfo["Chromosomes Info"][0].Length * data_Vue.geninfo["Chromosomes Info"][0].MD >= 1000000){
			data_Vue.warnings1.push(gen_warn_text);
		} 
		
		for (i=0; i<thisChromoNum; i++) {
			gen_warn_text = "Chromo "+(i+1)+" : "+"Please enter Chromosomes Length and must be in number.";
			if((!data_Vue.geninfo['Chromosomes Info'][i]['Length'] || isNaN(data_Vue.geninfo['Chromosomes Info'][i]['Length']) || data_Vue.geninfo['Chromosomes Info'][i]['Length'] < 0 ) & data_Vue.warnings1.indexOf(gen_warn_text) == -1){
				data_Vue.warnings1.push(gen_warn_text);
			}

			gen_warn_text = "Chromo "+(i+1)+" : "+"Please enter Chromosomes MD and must be in number.";
			if((!data_Vue.geninfo['Chromosomes Info'][i]['MD'] || isNaN(data_Vue.geninfo['Chromosomes Info'][i]['MD']) || data_Vue.geninfo['Chromosomes Info'][i]['MD'] < 0) & data_Vue.warnings1.indexOf(gen_warn_text) == -1){
				data_Vue.warnings1.push(gen_warn_text);
			}

			gen_warn_text = "Chromo "+(i+1)+" : "+"Please enter Chromosomes Recombination and must be in number.";
			if((!data_Vue.geninfo['Chromosomes Info'][i]['Recombination'] || isNaN(data_Vue.geninfo['Chromosomes Info'][i]['Recombination']) || data_Vue.geninfo['Chromosomes Info'][i]['Recombination'] < 0) & data_Vue.warnings1.indexOf(gen_warn_text) == -1){
				data_Vue.warnings1.push(gen_warn_text);
			}
		}
	}
	
	gen_warn_text = "You are not allowed to run more than 5 simulations in parallel";
	if(data_Vue.geninfo['number-simulations-parallel']!= undefined & data_Vue.geninfo['number-simulations-parallel']!="" & data_Vue.geninfo['number-simulations-parallel'] > 5){
		data_Vue.warnings1.splice(data_Vue.warnings1.indexOf(gen_warn_text),1);
	}
	
	data_Vue.warnings = data_Vue.warnings1.concat(data_Vue.warnings2, data_Vue.warnings3, data_Vue.warnings4, data_Vue.warnings5, data_Vue.warnings6);
	data_Vue.warnings = data_Vue.warnings.filter(Boolean);
}             
// end general info validation


	
// phenotype validation
if(id == "Phenotype_Info_Div" || id == "General_Info"){
	

	data_Vue.warnings2 = [];
	count = 0;
	for(key in data_Vue.warnings1) {
		if(data_Vue.warnings1.hasOwnProperty(key)) {
			count++;
		}	
	}
	
	if(count==0){
		data_Vue.warnings1 = [];
	}
	
	var phenoLength = data_Vue.traitsinfo.length;

	if (phenoLength > 0 ) {
	for (i=0; i<phenoLength; i++)
		{
			
			gen_warn_text = "No trait name provided for phenotype "+(i+1);
			if(!data_Vue.traitsinfo[i]['Trait Name']  & data_Vue.warnings2.indexOf(gen_warn_text) == -1) {
				data_Vue.warnings2.push(gen_warn_text);
				}

			
			gen_warn_text = "No phenotypic mean provided for phenotype "+(i+1);
			if((!data_Vue.traitsinfo[i]['Trait Mean'] || isNaN(data_Vue.traitsinfo[i]['Trait Mean'])) & data_Vue.warnings2.indexOf(gen_warn_text) == -1) {
				data_Vue.warnings2.push(gen_warn_text);
			}


			gen_warn_text = "No phenotypic standard deviation provided for phenotype "+(i+1);
			thisStdDev = data_Vue.traitsinfo[i]['Trait Std Deviation'];

			if ((!thisStdDev || isNaN(thisStdDev) || thisStdDev < 0)  & data_Vue.warnings2.indexOf(gen_warn_text) == -1) {
				data_Vue.warnings2.push(gen_warn_text);
			}

	
			gen_warn_text = "Heritability for phenotype"+(i+1)+" must be between 0 and 1";
			if((!data_Vue.traitsinfo[i]['Trait Heritability'] || data_Vue.traitsinfo[i]['Trait Heritability'] < 0 ||  data_Vue.traitsinfo[i]['Trait Heritability'] > 1) & data_Vue.warnings2.indexOf(gen_warn_text) == -1) {
				data_Vue.warnings2.push(gen_warn_text);
			}


			if(data_Vue.geninfo.advanced_trait && data_Vue.geninfo.advanced){
				
				if(data_Vue.geninfo.advanced_trait_repeat){
					gen_warn_text = "Repeatability for phenotype"+(i+1)+" must be between Heritability and 1";
					if(data_Vue.geninfo.advanced_trait & data_Vue.traitsinfo[i]['Trait Repeatability'] != "" & ( data_Vue.traitsinfo[i]['Trait Heritability'] > data_Vue.traitsinfo[i]['Trait Repeatability'] ||  data_Vue.traitsinfo[i]['Trait Repeatability'] > 1) & data_Vue.warnings2.indexOf(gen_warn_text) == -1) {
						data_Vue.warnings2.push(gen_warn_text);
					}
					
					gen_warn_text = "Repeatability set to Heritability"; 
					if(data_Vue.geninfo.advanced_trait & data_Vue.traitsinfo[i]['Trait Repeatability'] != "" & ( data_Vue.traitsinfo[i]['Trait Heritability'] > data_Vue.traitsinfo[i]['Trait Repeatability'] ||  data_Vue.traitsinfo[i]['Trait Repeatability'] > 1) & data_Vue.warnings2.indexOf(gen_warn_text) == -1) {
						data_Vue.warnings2.push(gen_warn_text);
					}
				}

				if(data_Vue.geninfo.advanced_trait_epi){
					checkPoly = isPositiveInt(data_Vue.traitsinfo[i].dominant_qtl);
					gen_warn_text = "Number of dominant QTL not provided for phenotype"+(i+1)
					if(data_Vue.geninfo.advanced_trait & checkPoly == false & data_Vue.warnings2.indexOf(gen_warn_text) == -1) {
						data_Vue.warnings2.push(gen_warn_text);
					}

					checkPoly = isPositiveInt(data_Vue.traitsinfo[i].qualitative_qtl);
					gen_warn_text = "Number of qualitative epistatic QTL not provided for phenotype"+(i+1)
					if(data_Vue.geninfo.advanced_trait &  checkPoly == false & data_Vue.warnings2.indexOf(gen_warn_text) == -1) {
						data_Vue.warnings2.push(gen_warn_text);
					}

					
					checkPoly = isPositiveInt(data_Vue.traitsinfo[i].quantitative_qtl);
					gen_warn_text = "Number of quantitative epistatic QTL not provided for phenotype"+(i+1)
					if((data_Vue.geninfo.advanced_trait & checkPoly == false & data_Vue.warnings2.indexOf(gen_warn_text) == -1)) {
						data_Vue.warnings2.push(gen_warn_text);
					}
				}
			}


			
			checkPoly = isPositiveInt(data_Vue.traitsinfo[i]['Trait Number of Polygenic Loci']);
			gen_warn_text = "Number of purely additive QTL not provided for phenotype"+(i+1)
			if(checkPoly == false & data_Vue.warnings2.indexOf(gen_warn_text) == -1) {
				data_Vue.warnings2.push(gen_warn_text);
			}

			gen_warn_text = "Pheno"+(i+1)+" : "+"Please enter Value per Unit and must be in number.";
			if((data_Vue.traitsinfo[i]['Trait Value per Unit'] == "" & data_Vue.traitsinfo[i]['Trait Value per Unit']  != "0") & data_Vue.warnings2.indexOf(gen_warn_text) == -1) {
				data_Vue.warnings2.push(gen_warn_text);
			}

							
			gen_warn_text = "Pheno"+(i+1)+" : "+"Please enter Major QTL and must be in positive number or can be zero. ";
			thisMajorQTL = data_Vue.traitsinfo[i]['Trait Major QTL']
			if(thisMajorQTL < 0 & data_Vue.warnings2.indexOf(gen_warn_text) ==-1){  
				data_Vue.warnings2.push(gen_warn_text);
			}

			var qtlCnt = data_Vue.traitsinfo[i]['Trait Major QTL'];
			var qtlName = i;

			if (qtlCnt > 0) {
				qtlsforPheno(qtlName, qtlCnt);
			}
				
			matrix1();
			matrix2();				

		}
	}
	data_Vue.warnings = data_Vue.warnings1.concat(data_Vue.warnings2, data_Vue.warnings3, data_Vue.warnings4, data_Vue.warnings5, data_Vue.warnings6);
	data_Vue.warnings = data_Vue.warnings.filter(Boolean);
}	// end of phenotype validation


if(id == "Selectionindex_Div" || id == "General_Info"){
	
	data_Vue.warnings6 = []
	
	if(data_Vue.geninfo.advanced && data_Vue.geninfo.advanced_sub){
		for( var i = 1; i < data_Vue.subpopulation.subpopulation_list.length; i++){
			var ktemp = Number(data_Vue.subpopulation.subpopulation_list[i].share0) + Number(data_Vue.subpopulation.subpopulation_list[i].share1);
			
			var count = 0; 
			for(var prop in data_Vue.subpopulation.subpopulation_list[i]){
				if(data_Vue.subpopulation.subpopulation_list[i].hasOwnProperty(prop)){
					++count;
				}
			}
			

				if(isNaN(ktemp)){
					gen_warn_text = " No valid input for share of fixated markers for subpopulation " + data_Vue.subpopulation.subpopulation_list[i].Name;
					data_Vue.warnings6.push(gen_warn_text);
				} else if(ktemp>1){
					gen_warn_text = " Share of fixated markers should be lower than 1 for subpopulation " + data_Vue.subpopulation.subpopulation_list[i].Name;

					data_Vue.warnings6.push(gen_warn_text);
				} else if(ktemp<0 ){
					gen_warn_text = " Share of fixated markers should be higher than 0 for subpopulation " + data_Vue.subpopulation.subpopulation_list[i].Name;
					data_Vue.warnings6.push(gen_warn_text);
				}

			
			if(count > 7){
				if(isNaN(ktemp)){
					gen_warn_text = " No valid input for share of fixated markers for subpopulation " + data_Vue.subpopulation.subpopulation_list[i].Name;
					data_Vue.warnings6.push(gen_warn_text);
				} else if(ktemp==1){
					gen_warn_text = " Share of fixated markers should be lower than 1 for subpopulation " + data_Vue.subpopulation.subpopulation_list[i].Name + " when applying deviation from mean";

					data_Vue.warnings6.push(gen_warn_text);
				} else if(ktemp== 0 ){
					gen_warn_text = " Share of fixated markers should be higher than 0 for subpopulation " + data_Vue.subpopulation.subpopulation_list[i].Name + " when applying deviation from mean";;
					data_Vue.warnings6.push(gen_warn_text);
				}
			}

		


		}
	}

	
	
	var traitsCnt = data_Vue.traitsinfo.length;
	var siCnt = data_Vue.selection_index.length

	if (siCnt > 0) {
		for(var i=0; i<siCnt; i++){
			for(var j=0; j < traitsCnt; j++)
				{

					var objSI = pVar =	data_Vue.selection_index[i];
					var keys = Object.keys(objSI);
					constP = keys[j+1];			
					pVar =	data_Vue.selection_index[i][constP];
					pName = data_Vue.selection_index[i]['Name'];
					checkpVar = isNumeric(pVar);	

					gen_warn_text = pName+(i+1)+":Pheno"+(j+1)+":Please enter selection index and it must be a number.";
					if ((!pVar || isNaN(pVar) || checkpVar == false) & data_Vue.warnings6.indexOf(gen_warn_text) == -1) {
						data_Vue.warnings6.push(gen_warn_text);
					}
				
			}
		}	
	}			

	var traitsCnt = data_Vue.traitsinfo.length;
	var piCnt = data_Vue.phenotyping_class.length

	if (piCnt > 0) {
	   for(var i=0; i<piCnt; i++){
		for(var j=0; j < traitsCnt; j++)
		{	
			var objPC = pVar =	data_Vue.phenotyping_class[i];
			var keys = Object.keys(objPC);
			
			constPC = keys[j+2];			
					
			pVar =  data_Vue.phenotyping_class[i][constPC];
			pName = data_Vue.phenotyping_class[i]['Name'];
			Cost = data_Vue.phenotyping_class[i]['Cost of phenotyping'];

			checkpVar = isNumeric(pVar);  
			checkPospVar = isPositiveInt(pVar);
				
			gen_warn_text = pName+":Phenotyping Cost"+(i+1)+": Please enter phenotyping classes and it must be a positive number or zero.";
			if ((isNaN(Cost) || Cost <0) & data_Vue.warnings6.indexOf(gen_warn_text) == -1) {
				data_Vue.warnings6.push(gen_warn_text);
			}

			gen_warn_text = pName+(i+1)+":Pheno"+(j+1)+": Please enter phenotyping classes and it must be a positive number or zero.";
			if ((pVar == null || isNaN(pVar) || pVar < 0 || checkpVar == false || checkPospVar == false) & data_Vue.warnings6.indexOf(gen_warn_text) == -1) {
				data_Vue.warnings6.push(gen_warn_text);
			}

		}
		}             
	}
		
	gen_warn_text = "Fixed costs are negative!";
	
	curFixedCost = data_Vue.economy['Fixed Cost'];
	console.log(curFixedCost);

	if(((!curFixedCost & curFixedCost !=0) || curFixedCost < 0 || isNaN(curFixedCost)) & data_Vue.warnings6.indexOf(gen_warn_text) == -1){
		data_Vue.warnings6.push(gen_warn_text);
		console.log('empty'+curFixedCost);
	}
	
	gen_warn_text = "Invalid Interest Rate provided. Assume interest of 0 percent.";
	curInt = data_Vue.economy['Interest Rate'];	

	if(((!curInt & curInt !=0) || curInt < (-100) || isNaN(curInt)) & data_Vue.warnings6.indexOf(gen_warn_text) == -1){
		data_Vue.warnings6.push(gen_warn_text);
		console.log('empty'+curInt);
	}
	
	gen_warn_text = "Genotyping costs are negative / not entered!";
	curGenoType = data_Vue.economy['Genotyping Cost'];

	if(((!curGenoType & curGenoType !=0) || curGenoType < 0 || isNaN(curGenoType)) & data_Vue.warnings6.indexOf(gen_warn_text) == -1){
		data_Vue.warnings6.push(gen_warn_text);
		console.log('empty'+curInt);
	}
	
	var animal_housing_costCnt = data_Vue.economy['Animal Housing Costs'].length;
	var thisAnimalHousingCosts = data_Vue.economy['Animal Housing Costs'];
	
	for(var i=0; i < animal_housing_costCnt; i++){			
		thisCost = thisAnimalHousingCosts[i]['Cost'];
		nameOfCost = thisAnimalHousingCosts[i]['Name'];

		gen_warn_text = "Housing costs for class " + nameOfCost + " are negative / not entered!";

		if(((thisCost =="" & thisCost !=0) || thisCost <0 || isNaN(thisCost)) & data_Vue.warnings6.indexOf(gen_warn_text) == -1){
			data_Vue.warnings6.push(gen_warn_text);
		}

		
	}
	
	data_Vue.warnings = data_Vue.warnings1.concat(data_Vue.warnings2, data_Vue.warnings3, data_Vue.warnings4, data_Vue.warnings5 , data_Vue.warnings6);
	data_Vue.warnings = data_Vue.warnings.filter(Boolean);
}


if(id == "Variables_Info" || id == "General_Info"){

	data_Vue.warnings3 = [];

	count = 0;
	for(key in data_Vue.warnings1) {
		if(data_Vue.warnings1.hasOwnProperty(key)) {
			count++;
		}	
	}
	
	if(count==0){
		data_Vue.warnings1 = [];
	}
	var ownVariableCnt = data_Vue.individualsVar_options.length;
	var ownVariable = data_Vue.individualsVar_options;

	for(var i=0; i < ownVariableCnt; i++){			
		thisName = ownVariable[i]['name'];
		thisValue = ownVariable[i]['value'];

		gen_warn_text = thisName+" : "+"Please enter Own variable and must be a positive number.";

		if((thisValue == "" || thisValue =="0" || thisValue <0 || isNaN(thisValue)) & data_Vue.warnings3.indexOf(gen_warn_text) == -1){
			data_Vue.warnings2.push(gen_warn_text);
		}
		else if((thisValue != "" & thisValue >0) & data_Vue.warnings3.indexOf(gen_warn_text) > -1){
			data_Vue.warnings3.splice(data_Vue.warnings3.indexOf(gen_warn_text),1); 
		}		
			
	}
	data_Vue.warnings = data_Vue.warnings1.concat(data_Vue.warnings2, data_Vue.warnings3, data_Vue.warnings4, data_Vue.warnings5 , data_Vue.warnings6);
	data_Vue.warnings = data_Vue.warnings.filter(Boolean);
}



if(id == "edge-popUp" || id == "node-popUp"){
	
	data_Vue.warnings4 = [];
	count = 0;
	for(key in data_Vue.warnings1) {
		if(data_Vue.warnings1.hasOwnProperty(key)) {
			count++;
		}	
	}
	
	if(count==0){
		data_Vue.warnings1 = [];
	}
	nodes = data_Vue.nodes.get();
	nodes = Array.from(nodes);
	if(nodes.length > 0){
		for(var i=0; i < nodes.length; i++){
			active_node = nodes[i];
			
			gen_warn_text = "Please assign each node with a name.";
			if (!active_node['id'] & data_Vue.warnings4.indexOf(gen_warn_text) == -1){			
				data_Vue.warnings4.push(gen_warn_text);
			}
			
			gen_warn_text = "Number of individuals in Node " + active_node['id'] + " is not specified";
			if (!active_node['Number of Individuals'] & data_Vue.warnings4.indexOf(gen_warn_text) == -1){			
				data_Vue.warnings4.push(gen_warn_text);
			}
			
			gen_warn_text = "Share of genotyped individuals in " + active_node['id'] + " must be between 0 and 1.";
			if((!active_node['Proportion of genotyped individuals'] || active_node['Proportion of genotyped individuals'] < 0 || active_node['Proportion of genotyped individuals']  > 1) & data_Vue.warnings4.indexOf(gen_warn_text) == -1) {
				data_Vue.warnings4.push(gen_warn_text);
			}
			
			gen_warn_text = "Share of male individuals in " + active_node['id'] + " must be between 0 and 1.";
			if((active_node['Sex'] == "Both" & (!active_node['Proportion of Male'] || active_node['Proportion of Male'] < 0 || active_node['Proportion of Male']  > 1)) & data_Vue.warnings4.indexOf(gen_warn_text) == -1) {
				data_Vue.warnings4.push(gen_warn_text);
			}
			
			gen_warn_text = "Node " + active_node['id'] + " only contains male individuals. Adapt chosen sex!";
			if((active_node['Sex'] == "Both" & (!active_node['Proportion of Male'] || active_node['Proportion of Male']  == 1)) & data_Vue.warnings4.indexOf(gen_warn_text) == -1) {
				data_Vue.warnings4.push(gen_warn_text);
			}
			
			gen_warn_text = "Node " + active_node['id'] + " only contains female individuals. Adapt chosen sex!";
			if((active_node['Sex'] == "Both" & (!active_node['Proportion of Male'] || active_node['Proportion of Male']  == 0)) & data_Vue.warnings4.indexOf(gen_warn_text) == -1) {
				data_Vue.warnings4.push(gen_warn_text);
			}
			
		}
	}
	
	data_Vue.warnings = data_Vue.warnings1.concat(data_Vue.warnings2, data_Vue.warnings3, data_Vue.warnings4, data_Vue.warnings5 , data_Vue.warnings6);
	data_Vue.warnings = data_Vue.warnings.filter(Boolean);
}


if(id == "edge-popUp" || id == "node-popUp"){
	
	count = 0;
	for(key in data_Vue.warnings1) {
		if(data_Vue.warnings1.hasOwnProperty(key)) {
			count++;
		}	
	}
	
	if(count==0){
		data_Vue.warnings1 = [];
	}
	
	data_Vue.warnings5 = [];

			
	nodes = data_Vue.nodes.get();
	nodes = Array.from(nodes);
	edges = data_Vue.edges.get();
	edges = Array.from(edges);
	
	if(edges.length > 0){
		edge_from = [];
		edge_to = [];
		edge_nrfrom = [];
		edge_nrto = [];
		edge_type = [];
		node_name = [];
		node_founder = [];
		
		for(var i=0; i < nodes.length; i++){
			node_name.push(nodes[i]['id'])
			node_founder.push(nodes[i]['Founder'])
		}
		
		for(var i=0; i < edges.length; i++){
			edge_to.push(edges[i]['to']);
			edge_from.push(edges[i]['from']);
			edge_type.push(edges[i]['Breeding Type']);
			for(var j=0; j < nodes.length; j++){
				if(edge_from[i]==node_name[j]){
					edge_nrfrom.push(j);
				}
				if(edge_to[i]==node_name[j]){
					edge_nrto.push(j);
				}
			}
			gen_warn_text = "No Breeding Type selected for Edge between " + edges[i]['from'] + " and " + edges[i]['to'];
			if (edges[i]['Breeding Type']=="" & data_Vue.warnings5.indexOf(gen_warn_text) == -1){			
				data_Vue.warnings5.push(gen_warn_text);
			}
		}
		
		
		for(var i=0; i < edges.length; i++){
			if(edge_type[i] =="Selection" || edge_type[i] == "Aging" || edge_type[i] == "Split"){
				size1 = parseFloat(nodes[edge_nrfrom[i]]['Number of Individuals']);
				size2 = parseFloat(nodes[edge_nrto[i]]['Number of Individuals']);
				if(size1<size2){
					gen_warn_text = "More individuals selected in " + edge_to[i] + " than present in " + edge_from[i];
					data_Vue.warnings5.push(gen_warn_text);
				}
			}
			if(edge_type[i] =="Selection" || edge_type[i] == "Aging" || edge_type[i] == "Split"){
				if(nodes[edge_nrfrom[i]]['Sex']!=nodes[edge_nrto[i]]['Sex'] & nodes[edge_nrto[i]]['Sex'] != "Both"){
					gen_warn_text = "Different sex between nodes " + edge_to[i] + " and " + edge_from[i];
					data_Vue.warnings5.push(gen_warn_text);
				}
			}
			if(edge_type[i] =="Repeat"){
				if(nodes[edge_nrfrom[i]]['Number of Individuals']!=nodes[edge_nrto[i]]['Number of Individuals']){
					gen_warn_text = "Different number of individuals in nodes of the Repeat between " + edge_to[i] + " and " + edge_from[i];
					data_Vue.warnings5.push(gen_warn_text);
				}
			}
		}
		
		for(var i=0; i < nodes.length; i++){
			if(node_founder[i] == 'Yes'){
				for(var j=0; j <edges.length; j++){
					if(edge_nrto[j] == i && edge_type[j] != 'Repeat'){
						gen_warn_text = "Founder-Node " + node_name[i] + " has incoming Edges.";
						data_Vue.warnings5.push(gen_warn_text);	
					}
				}
			}
			else{
				count = 0;
				combine_count = 0;
				split_count = 0;
				types = [];
				for(var j=0; j < edges.length; j++){
					if(edge_nrto[j] == i){
						count++;
						if(edge_type[j] == "Combine"){
							size1 = parseFloat(nodes[edge_nrfrom[j]]['Number of Individuals']);
							combine_count = combine_count + size1;
						}
					}
					if(edge_nrfrom[j] == i){
						if(edge_type[j] == "Split"){
							size1 = parseFloat(nodes[edge_nrto[j]]['Number of Individuals']);
							split_count = split_count + size1;
						}
					}
					
					if(edge_nrto[j] == i){
						if(types.length==0 & edge_type[j] != "Repeat"){
							types = edge_type[j];
						}

						else if(types != edge_type[j] & edge_type[j] != "Repeat"){
							data_Vue.warnings5.push(types)
							data_Vue.warnings5.push(edge_type[j])
							gen_warn_text = "Node "+ node_name[i] +" is generated by different breeding types. This is no intended/supported structure";
							data_Vue.warnings5.push(gen_warn_text);	
						
						}
					}
					
				}

				
				if(count==0){
					gen_warn_text = "Node "+ node_name[i] +" is no Founder but has no incoming Edges";
					data_Vue.warnings5.push(gen_warn_text);	
				}
				size2 = nodes[i]['Number of Individuals'];
				if(combine_count>0 && combine_count != size2){
					gen_warn_text = "Individual number for Combining nodes to "+ node_name[i] +" do not add up";
					data_Vue.warnings5.push(gen_warn_text);	
				}
				if(split_count>0 && split_count != size2){
					gen_warn_text = "Individual numbers for Splitting node "+ node_name[i] +" do not add up";
					data_Vue.warnings5.push(gen_warn_text);	
				}
			}
		}
		
		
		for(var j=0; j < edges.length; j++){
			if(edges[j]['Selection Type'] =="BVE" && (edges[j]['Cohorts used in BVE'] != undefined && edges[j]['Cohorts used in BVE'] != "Manual select")){
				gen_warn_text = "It's highly suggested to use Manual Select to selected cohorts used for BVE (Revisit Edge between: " + edges[j]['from'] + " and " + edges[j]['to'] +")";
				data_Vue.warnings5.push(gen_warn_text);	
			}
		}
		
		
	}
	


	
	data_Vue.warnings = data_Vue.warnings1.concat(data_Vue.warnings2, data_Vue.warnings3, data_Vue.warnings4, data_Vue.warnings5, data_Vue.warnings6);
	data_Vue.warnings = data_Vue.warnings.filter(Boolean);

	
	if(data_Vue.edges._data.ZuchthengsteCopy4598_Zuchthengste != undefined && data_Vue.edges._data.ZuchthengsteCopy4598_Zuchthengste["Manuel selected cohorts"] != undefined){
		for(var j = 0; j < data_Vue.edges._data.ZuchthengsteCopy4598_Zuchthengste["Manuel selected cohorts"].length; j++){
			if(data_Vue.edges._data.ZuchthengsteCopy4598_Zuchthengste["Manuel selected cohorts"][j]== "stallion licensingOCD:-4"){
				data_Vue.edges._data.ZuchthengsteCopy4598_Zuchthengste["Manuel selected cohorts"][j] = "stallion licensing OCD:-4"
			}
			if(data_Vue.edges._data.ZuchthengsteCopy4598_Zuchthengste["Manuel selected cohorts"][j]== "breeding marespr??fung:-5"){
				data_Vue.edges._data.ZuchthengsteCopy4598_Zuchthengste["Manuel selected cohorts"][j] = "mare performance test:-5"
			}
			if(data_Vue.edges._data.ZuchthengsteCopy4598_Zuchthengste["Manuel selected cohorts"][j]== "breeding marespr??fung:-6"){
				data_Vue.edges._data.ZuchthengsteCopy4598_Zuchthengste["Manuel selected cohorts"][j] = "mare performance test:-6"
			}
			if(data_Vue.edges._data.ZuchthengsteCopy4598_Zuchthengste["Manuel selected cohorts"][j]== "breeding marespr??fung:-7"){
				data_Vue.edges._data.ZuchthengsteCopy4598_Zuchthengste["Manuel selected cohorts"][j] = "mare performance test:-7"
			}
		}

	}

}


data_Vue.project_saved = false;
}


function isNumeric(val) {
  return !isNaN(parseFloat(val)) && isFinite(val);
}

function isPositiveInteger(val) {
	var isNumber = Math.floor(Number(val));
	var thisNumber = (String(isNumber) == val && isNumber >= 0);
	return thisNumber;
}


function isPositiveInt(val) {
	return val >>> 0 === parseFloat(val);
}


function qtlsforPheno(qtlName, qtlCnt) {

	var thisQTL = data_Vue.traitsinfo[qtlName]["Trait QTL Info"];
	var thisTrait = data_Vue.traitsinfo[qtlName]["Trait Name"];
	
	for (k=0; k<qtlCnt; k++) {	
		curBP = thisQTL[k]['QTL BP'];
		curID = thisQTL[k]['QTL ID'];
		gen_warn_text = thisTrait+"-bp"+(k+1)+":"+"Please enter bp and must be a Positive Number .";	
		checkBP = isPositiveInt(thisQTL[k]['QTL BP']);
		if ((!thisQTL[k]['QTL BP'] || thisQTL[k]['QTL BP'] <= 0 || checkBP == false )& data_Vue.warnings.indexOf(gen_warn_text) == -1) {
			data_Vue.warnings.push(gen_warn_text);
		}
		else if ((thisQTL[k]['QTL BP'] != "" & thisQTL[k]['QTL BP'] > 0 & checkBP == true) & data_Vue.warnings.indexOf(gen_warn_text) > -1){                                                                    
			data_Vue.warnings.splice(data_Vue.warnings.indexOf(gen_warn_text),1); 
		}		

		gen_warn_text = thisTrait+"-Chromo "+(k+1)+":"+"Please enter chromosomes and must be a Positive Number.";
		checkchromo = isPositiveInt(thisQTL[k]['QTL Chromosome']);
		if ((!thisQTL[k]['QTL Chromosome'] || thisQTL[k]['QTL Chromosome'] <= 0 || checkchromo == false) & data_Vue.warnings.indexOf(gen_warn_text) == -1) {
			data_Vue.warnings.push(gen_warn_text);
		}
		else if ((thisQTL[k]['QTL Chromosome'] != "" & thisQTL[k]['QTL Chromosome'] > 0 & checkchromo ==true) & data_Vue.warnings.indexOf(gen_warn_text) > -1){
			data_Vue.warnings.splice(data_Vue.warnings.indexOf(gen_warn_text),1); 
		}		
		
		checkAANumeric = isNumeric(thisQTL[k]['QTL Effect AA']);
		gen_warn_text = thisTrait+"-EffectAA"+(k+1)+":"+"Please enter Effect AA and must be a Number.";
		if ((!thisQTL[k]['QTL Effect AA'] || checkAANumeric == false ) & data_Vue.warnings.indexOf(gen_warn_text) == -1) {
			data_Vue.warnings.push(gen_warn_text);
		}
		else if (thisQTL[k]['QTL Effect AA'] != "" & checkAANumeric == true & data_Vue.warnings.indexOf(gen_warn_text) > -1){                                                                       
			data_Vue.warnings.splice(data_Vue.warnings.indexOf(gen_warn_text),1); 
		}	
		
		checkABNumeric = isNumeric(thisQTL[k]['QTL Effect AB']);	
		gen_warn_textAB = thisTrait+"-EffectAB"+(k+1)+":"+"Please enter Effect AB and must be a Number.";
		if ((!thisQTL[k]['QTL Effect AB'] || checkABNumeric == false) & data_Vue.warnings.indexOf(gen_warn_textAB) == -1) {
			data_Vue.warnings.push(gen_warn_textAB);
		}
		else if ((thisQTL[k]['QTL Effect AB'] != "" ) & checkABNumeric == true & data_Vue.warnings.indexOf(gen_warn_textAB) > -1){                                                                       
			data_Vue.warnings.splice(data_Vue.warnings.indexOf(gen_warn_textAB),1); 
		}

		
		checkBBNumeric = isNumeric(thisQTL[k]['QTL Effect BB']);	
		gen_warn_textBB = thisTrait+"-EffectBB"+(k+1)+":"+"Please enter Effect BB and must be a Number.";
		if ((!thisQTL[k]['QTL Effect BB'] || checkBBNumeric == false) & data_Vue.warnings.indexOf(gen_warn_textBB) == -1) {
			data_Vue.warnings.push(gen_warn_textBB);
		}
		else if ((thisQTL[k]['QTL Effect BB'] != "") & checkBBNumeric == true & data_Vue.warnings.indexOf(gen_warn_textBB) > -1){                                                                       
			data_Vue.warnings.splice(data_Vue.warnings.indexOf(gen_warn_textBB),1); 
		}
				
		checkAllele = isNumeric(thisQTL[k]['QTL Allele Frequency']);		
		gen_warn_text = thisTrait+"-Allele Freq "+(k+1)+":"+"Please enter Allele Freq and it must be between 0 and 1.";
		if ((!thisQTL[k]['QTL Allele Frequency']  || thisQTL[k]['QTL Allele Frequency'] <0 || thisQTL[k]['QTL Allele Frequency'] >1 ) & data_Vue.warnings.indexOf(gen_warn_text) == -1) {
			data_Vue.warnings.push(gen_warn_text);
		}
		else if ((thisQTL[k]['QTL Allele Frequency'] != "" & thisQTL[k]['QTL Allele Frequency'] >=0 & thisQTL[k]['QTL Allele Frequency'] <=1) & data_Vue.warnings.indexOf(gen_warn_text) > -1){                                                                       
			data_Vue.warnings.splice(data_Vue.warnings.indexOf(gen_warn_text),1); 
		}	
				
	}
}

function matrix1() {
	var traitsCnt = data_Vue.traitsinfo.length;
	for(var i=0; i < traitsCnt; i++){			
		for(var j=0; j < data_Vue.matrix[i].row.length; j++){
			curMatVal = data_Vue.matrix[i].row[j].val;	
				
			gen_warn_text = data_Vue.traitsinfo[i]['Trait Name']+"-"+data_Vue.traitsinfo[j]['Trait Name']+": Phenotypic Correlation must be between -1 and 1 ";

			if (((!curMatVal & curMatVal == "" & curMatVal !="0") || curMatVal <-1 || curMatVal > 1 || isNaN(curMatVal)) & data_Vue.warnings.indexOf(gen_warn_text) == -1) {
				data_Vue.warnings.push(gen_warn_text);
			}
			else if ((curMatVal != "" & curMatVal >=-1 & curMatVal <=1 & isNaN(curMatVal) == false ) & data_Vue.warnings.indexOf(gen_warn_text) > -1){
				data_Vue.warnings.splice(data_Vue.warnings.indexOf(gen_warn_text),1); 
			}
			
		}			
	}			
}	


function matrix2() {
		var traitsCnt = data_Vue.traitsinfo.length;
		for(var i=0; i < traitsCnt; i++){			
			for(var j=0; j < data_Vue.matrix2[i].row.length; j++){
				curMatVal2 = data_Vue.matrix2[i].row[j].val;	

				gen_warn_text = data_Vue.traitsinfo[i]['Trait Name']+"-"+data_Vue.traitsinfo[j]['Trait Name']+": Genetic Correlation must be between -1 and 1 ";
				
				if (((!curMatVal2 & curMatVal2 == "" & curMatVal2 !="0") || curMatVal2 <-1 || curMatVal2 > 1 || isNaN(curMatVal2)) & data_Vue.warnings.indexOf(gen_warn_text) == -1) {
					data_Vue.warnings.push(gen_warn_text);
				}
				else if ((curMatVal2 != "" & curMatVal2 >=-1 & curMatVal2 <=1 & isNaN(curMatVal2) == false ) & data_Vue.warnings.indexOf(gen_warn_text) > -1){
					data_Vue.warnings.splice(data_Vue.warnings.indexOf(gen_warn_text),1); 
				}
				
			}			
		}			
}

			




// this funtion will always be called after an edge has been created,
// you can use this function to add further services, e.g add the secong parent automatically if Breeding Type == Reproduction etc.
function addServices(){
if(data_Vue.active_edge['Breeding Type'] == "Reproduction"){
	var val = document.getElementById("parents_input").value ;
	if(val =='') return;
	if(data_Vue.nodes.get(val) == null){  	// create new node by copying:
		var p1_node = data_Vue.nodes.get(data_Vue.active_edge.from);
		p1_node.Sex = p1_node.Sex == "Male" ? "Female" : "Male" ;
		p1_node.id = val;
		p1_node.x = p1_node.x+50;
		addNode_extern(p1_node);
	}								
	var child_node = data_Vue.nodes.get(data_Vue.active_edge.to);
	var new_edge = new myEdge(val, child_node.id);
	new_edge['Breeding Type'] = "Reproduction";
	new_edge['Time Needed'] = document.getElementById("time_n").value;
	for(var i=0; i < data_Vue.active_edge.useVar.length; i++){
		new_edge['useVar'].push(data_Vue.active_edge.useVar[i]);
	}
	
	addEdge_extern(new_edge);
}
if(data_Vue.active_edge['Breeding Type'] == "Split"){
	var val = document.getElementById("split_input").value ;
	if(val =='') return;
	if(data_Vue.nodes.get(val) == null){  	// create new node by copying:
		var p1_node = data_Vue.nodes.get(data_Vue.active_edge.to);
		p1_node.id = val;
		p1_node.x = p1_node.x+50;
		p1_node['Number of Individuals'] = data_Vue.nodes.get(data_Vue.active_edge.from)['Number of Individuals'] - p1_node['Number of Individuals'];
		p1_node.individualsVar ='';
		addNode_extern(p1_node);
	}								
	var parent_node = data_Vue.nodes.get(data_Vue.active_edge.from);
	var new_edge = new myEdge(parent_node.id, val);
	new_edge['Breeding Type'] = "Split";
	
	addEdge_extern(new_edge);
}
}
