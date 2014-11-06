This project contains 2 files oneDJury.R and util.R. It only contains the functions and does not contain the procedures for the real clustering analysis and the ploting. These two scripts were used for Xuan's master thesis, Chap 2 and Chap 3.
* __SomeDiscription denotes used show provide proper data

1. For Chap2(evaluating discretizing gene expression data), the following pipline is suggested:
	1.0 load the util.R into the environment:
			source("util.R");
	1.1 download the data from genomic portals(http://eh3.uc.edu/GenomicsPortals/), also need to know the types of samples, which should exactly match the colnames of the samples before the first "." .Read data into the environment: 
			sData = read.csv(sourceFilePass, sep = "\t", header=TRUE)
	1.2 process the header infomation: 
			hInfo = colnames(sData)
			hInfo = unlist(lapply(hInfo, processHeader, sampleTypes));
			colnames(sData) = hInfo
		select the data of interest:
			rownames(sData) = sData[,1]
			sData = sData[, hInfo > 0]
	1.3 normalize and discretize data:
		normalization(returning a TRANSPOSED matrix):
			sourceData = apply(sData, MARGIN = 1, __normalizationFunctionUsed)
		discretization(returning a TRANSPOSED matrix)
			sourceData = apply(sData, MARGIN = 1, __discretizationFunctionUsed, __numOfStatesToDiscretizeInto)
	1.4 calculate the dist object (for normalized and discretized data; for original data, need to transpose it first, because dist is calcuated between rows of samples)
		For Euclidean distance:
			distances = dist(__data)
		For Hamming distance:
			distances = hammingMatrix(__data)
		For invertedCosineSimilarity distance:
			distances = cosineMatrix(__data)
	1.5 hierarchical clustering and cut the tree
			fit = hclust(distances)
			groups = cutree(fit, __numOfClustersAsOutPut)
	1.6 evaluate clustering result against the true labeling:
			trueLabeling = hInfo[hInfo > 0]
			result = RRand(groups, trueLabeling)$adjRand
	1.7 plot the data as needed.
2. For Chap3(oneDJury dimension reduction):
	2.0 loaddata into environment:
			source("oneDJury.R");
		if combinedMatrix is already saved
			load("__savedCombinedMatrixFIle")
	2.1 discretize the combinedMatrix by gene
			sourceData = apply(combinedMatrix, MARGIN = 1, __discretizationFunctionOfChoice, __numOfStatesToDiscretizeInto)
			rownames(sourceData) = as.numeric(colnames(combinedMatrix))
	2.2 get the gene programs
			geneLists = gpLists(__geneProgramsFileName)
			sourceData = findCombinedList(sourceData, geneLists)
			geneIndexs = findIndexesFromCombined(sourceData, geneLists)
	2.3 calculate oneDJury:
			sourceData = oneDJuryScore(sData = sourceData, genesInGroups = geneIndexs, numOfStates = __theNumOfStatesTODiscretizedInto);
			rownames(sourceData) = as.numeric(colnames(combinedMatrix))
	2.4 clustering
			distances = dist(sourceData)
			fit = hclust(distances)
			groups = cutree(fit, __numOfClustersAsOutPut)
	2.5 evaluate clustering result
			trueLabel = as.numeric(colnames(combinedMatrix))
			result = RRand(trueLabel, groups) $adjRand
	2.6 plot the data as needed



