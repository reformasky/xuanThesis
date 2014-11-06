###############################################
## all the utility functions
###############################################

# for adjusted rand index
library(phyclust, quiet = TRUE)

###############################################
## process the header infomation for data downloaded from Genomic portals contains 
## extract the sample types from header file
## for pValues, return 0;
## for other meta infomations return -1;
##    str: the an element of the header; 
##    sampleTypes:types of samples need to kept(corespond to the substring before "." in str)
## usage : headInfo = unlist(lapply(header, processHeader, sampleTypes))
###############################################
processHeader = function(str, sampleTypes) {	
	stringVec = unlist(strsplit(str, split = ".", fixed = TRUE));
	if (length(stringVec) >= 2){
		if(stringVec[1] %in% sampleTypes)
				which(sampleTypes == stringVec[1])
		else
			-1
	}
	#only useful for selecting the features with small pValue
	else if (str == "Pvalues"){
		0;
	}
	# mark as discard data
	else {
		-1
	}
}

##############################################
## normalizations functions: 
##	normalizationLinear: feature scaling;
##  normalizationZScore: z score normalization
##############################################

normalizationLinear = function(vec) {
	vec = as.numeric(vec)
	minVec = min(vec);
	maxVec = max(vec);
	(vec - minVec)/ (maxVec - minVec);
}

normalizationZScore = function(vec) {
	vec = as.numeric(vec)
	sdVec = sd(vec);
	meanVec = mean(vec);
	(vec - meanVec) / sdVec;
}

##############################################
## discretization functions: 
##	discretizationQuantile: equal frequency
##	discretizationFloor: equal interval
##  discretizationZScore: Z Score
##############################################
discretizationQuantile = function(vec, numOfStates = 3) {
	vec = rank(vec);
	vec = floor(vec * numOfStates / length(vec));
	vec[vec == numOfStates] = numOfStates - 1;
	vec;
}

discretizationFloor = function(vec, numOfStates = 3) {
	vec = normalizationLinear(vec);
	vec = floor(vec * numOfStates);
	vec[vec == numOfStates] = numOfStates - 1;
	vec;
} 


discretizationZScore = function(vec, numOfStates = 3) {
	result = floor( pnorm(normalizationZScore(vec)) * numOfStates );
	result[result >= numOfStates] = numOfStates - 1;
	result;
}

############################################################
## Distance objects for hamming distance and invertedCosineSimilarity
############################################################

#calculate the hamming distance between TRANSPOSED of discretized expression data,
#returns a dist object
hammingMatrix = function(mtx) {
	hammingDistance = function(vec1, vec2) {
		sum(as.integer(vec1) != as.integer(vec2))
	}

	result = matrix(0, nrow = dim(mtx)[1], ncol = dim(mtx)[1])
	colnames(result) = rownames(result) = rownames(mtx);
	for( i in 1 : (dim(mtx)[1] - 1) ){
		for(j in (i + 1) : dim(mtx)[1]) {
			result[i,j] = result[j,i] = hammingDistance(mtx[i,], mtx[j,])
		}
	}
	as.dist(result)
}

#calculate the inverted cosine similarity between TRANSPOSED of discretized expression data,
#returns a dist object
cosineMatrix = function(mtx) {
	cosineSimilarity = function(vec1 , vec2) {
		(vec1 %*% vec2) / sqrt(vec1 %*% vec1) / sqrt(vec2 %*% vec2)
	}
	result = matrix(1, nrow = dim(mtx)[1], ncol = dim(mtx)[1])
	colnames(result) = rownames(result) = rownames(mtx);
	for( i in 1 : (dim(mtx)[1] - 1) ){
		for(j in (i + 1) : dim(mtx)[1]) {
			result[i,j] = result[j,i] =1/ cosineSimilarity(mtx[i,], mtx[j,])
		}
	}
	as.dist(result)
}