#####################################
# prepare configureations 
# make sure all folders exists!!!
#  ./thesis
#  ./thesis/data
####################################
set.seed(1024)
setwd("./thesis")


#####################################
#   Download data from cBioPortal
#####################################

#converte geneId to HUGO sambols
library("cgdsr")
library(org.Hs.eg.db)
library(annotate)

#GPLists.csv contains lists of 22 lists of genes
geneListFile = "GPLists.csv"
openedFile = readLines(file(geneListFile, open = "r"))
geneIds = c()
for(i in 1: length(openedFile)) {
	geneIds = union(geneIds, unlist(strsplit(openedFile[i], ",")))
}
geneIds = geneIds[!(geneIds == "")]
geneList = getSYMBOL(geneIds, data='org.Hs.eg')	

mycgds = CGDS("http://www.cbioportal.org/public-portal/")
#studies Id for chromphobe, rental clear cell, ovarian servious and prostate ad (28, 30,49,55)
studiesIndex = c(28, 30,49,55)
#studies discreption for chromphobe, rental clear cell, ovarian servious and prostate ad
studiesId = c("kich_tcga", "kirc_tcga_pub","ov_tcga","prad_tcga")
sampleSize = c(66, 392 + 66, 66+392+158, 66+392+158+ 246)
cummulated = c(0,66,66 + 392, 66+392+158) + 1


# download data
for(i in 1 : length(studiesId)){
	print("**********")
	print(i)
	print("**********")
	
	# [1, 1] returns a unique ID used to identify the case list ID in subsequent interface calls
	mycaselist = getCaseLists(mycgds,studiesId[i])[1, 1]
	#[3,1] returns mRNA seq result
	mygeneticprofile = getGeneticProfiles(mycgds,studiesId[i])[3,1]

	#down load all 6800 genes 100 genes a time
	for(j in 1 : 68){
		if(j == 1)
		downloadData = getProfileData(mycgds,geneList[(j-1) * 100 + 1 : 100],mygeneticprofile,mycaselist)
		else
		downloadData = cbind(downloadData, getProfileData(mycgds,geneList[(j-1) * 100 + 1 : 100],mygeneticprofile,mycaselist))
		print(j)
	}
	#remove duplicated data? "XK"
	downloadData = downloadData[, colnames(downloadData) != "XK"]
	downloadData = t(downloadData)

	#remove na, nan, inf
	downloadData = downloadData[complete.cases(downloadData * 0), ,drop = F]

	colnames(downloadData) =  replicate(dim(downloadData)[2], i)
	downloadData = data.frame(id = rownames(downloadData), downloadData)
	if(i == 1) {
		combinedMatrix = downloadData
	}
	else {
		#merge with existing data.
		combinedMatrix = merge(combinedMatrix, downloadData, by = "id")
	}
}
# save data for further use in ./thesis/data
save(file = file.path("data","combinedMatrix"), combinedMatrix)

##############################################
## adjust data from different studies into 
## zeroMeanOnesd
##############################################
rownamesMtx = rownames(combinedMatrix) 
zeroMeanOneSd = function(mtx) {	
	means = apply(mtx, MARGIN =2, mean)
	subtract = function(vec, means) {
		vec - means;
	}
	mtx = t(apply(mtx, MARGIN =1, subtract, means))

	suqardSum  = apply(mtx, MARGIN =2,  function(vec){sum(vec **2)})

	divide = function(vec, suqardSum) {
		vec / sqrt(suqardSum)
	}

	mtx = t(apply(mtx, MARGIN = 1, divide, suqardSum))
	mtx = mtx * sqrt(dim(mtx)[1] -1)
}
# adjust data
combinedMatrix = zeroMeanOneSd(combinedMatrix)
rownames(combinedMatrix) =  rownamesMtx

#################################################
## intersecting combined matrix for gene programs
#################################################

#read geneIds from file into lists of lists
gpLists = function(fileName) {
  lines = readLines(file(fileName, open = "r"))
  results = list()
  for( i in 1 : length(lines)) {
    tempList = unlist(strsplit(lines[i], ","))
    result = tempList[ tempList != "" ]
    results = c(results, list(result))
  }
  results
}

geneListFileName = "GPLists.csv"
geneLists = gpLists(geneListFileName)
# returns a matrix which is a subset of sData, containing all geneIds appear in lists
# and remove the duplicated genes in the geneLists
findCombinedList = function(sData, geneLists) {
  geneIds = rownames(sData)
  selectedGeneId = c()
  for(geneList in geneLists) {
    geneList = unlist(geneList)
    indexes = intersect(geneList, geneIds)
    selectedGeneId = union(selectedGeneId, indexes)
  }

  indexes = match(selectedGeneId, geneIds)
  sData[indexes,]
}
combinedMatrix = findCombinedList(combinedMatrix, geneLists)

# for given gene ids in lists, find their coresponding row indexes in sData
# return a list of lists of indexes
findIndexesFromCombined = function(sData, lists) {
  results = list();
  vec = rownames(sData)
  for( i in 1 : length(lists) ) {
    lhs = unlist(lists[i])
    temp = match(lhs, vec);
    result = temp[!is.na(temp)]
    if(length(result) > 0)
    	results = c(results, list(result))
  }
  results
}

# assuming mtx is a matrix with proper discretization; returns nStates x numOfTotalGenes;mtx is 
# TRANSPOSED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# because of discretization! numOfSamples x numOfGenes
oneDJuryCache = function(mtx, numOfStates) {
	result = matrix(0, ncol = dim(mtx)[2], nrow = numOfStates);
	rownames(result) = 1 : numOfStates ;
	colnames(result) = colnames(mtx);
	for(r in 1: dim(result)[1]) {
		for(c in 1 : dim(result)[2]) {
			result[r, c] = sum(mtx[,c] == r );
		}
	}
	result;
}

# sData: matrix of TRANSPOSED discretized expression ; samples X gene;
# genesInGroups: a list of lists of genes, each list is the column indexes for that particular lists of genes
#	in the sData
# returns the numOfSample x numOfGroups matrix, each row is a 1DJury score vector.
oneDJuryScore = function(sData, genesInGroups, numOfStates = 3) {

	numOfGroups = length(genesInGroups);
	sourceData = sData + 1
	cachedStates = oneDJuryCache(sourceData, numOfStates);

	#translate [state, gene] score -> [sample, gene] score
	cachedSamples = matrix(0, nrow = dim(sourceData)[1], ncol = dim(sourceData)[2])
	for(r in 1 : dim(cachedSamples)[1]){
		for( c in 1 : dim(cachedSamples)[2]){
			cachedSamples[r,c] = cachedStates[sourceData[r,c], c];
		}
	}
	# numOfSamples x numOfGroups
	result = matrix(0, nrow = dim(sourceData)[1], ncol = numOfGroups);
	rownames(result) = rownames(sourceData);
	colnames(result) = 1 : numOfGroups
	for(group in 1 : numOfGroups) {
		selectedGenes = unlist(genesInGroups[group])
		result[,group] = t(rowSums(cachedSamples[, selectedGenes]))/ length(selectedGenes)
	}	
	result;
}
