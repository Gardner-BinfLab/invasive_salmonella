library("randomForest")
library("DMwR")

load("finalmodel.Rdata")
lookup <- read.delim("lookup.tsv")
strain <- read.table("strain", comment.char="@")

# Reading in the eggNOG model scores, checking to make sure that top eggNOG model hit for each protein in the orthogroup is the same
read <- read.delim("hmmsearch.parse", stringsAsFactors=F, comment.char="@", header=F)
bitscores <- cbind.data.frame(strain="strain", read) 

testdata <- data.frame()

scores <- vector()
for (j in lookup$Gene) { # for each key indicator gene
  if(! j %in% bitscores$V1) {
    scores <- c(scores, NA)
  } else if(bitscores$V7[bitscores$V1==j&bitscores$V4==lookup$Model[lookup$Gene==j]][1]<0.0001) { # if the E-value looks reasonable
    scores <- c(scores, bitscores$V8[bitscores$V1==j&bitscores$V4==lookup$Model[lookup$Gene==j]][1])
  } else {
    scores <- c(scores, NA)
  }
}
testdata <- rbind(testdata, scores)

names(testdata) <- colnames(train5)[-197]
row.names(testdata) <- unique(bitscores$strain)

# impute scores for genes not covered by enough reads
merge <- rbind(testdata, train5[,-197])
merge2 <- knnImputation(merge, k=1)

testdata <- merge2[1,]
# testdata[is.na(testdata)] <- 0 # use with caution, only with high-coverage genomes

# if a gene is missing from your samples, it will be filled in based on the median observed in our training data
# for(i in 1:ncol(testdata)) {
#   if(sum(is.na(testdata[,i])) == nrow(testdata)) {
#     testdata[,i] <- median(train5[,names(testdata)[i]])
#   }
# }

# load the random forest model and score your isolates
votes <- data.frame(predict(model5, na.roughfix(testdata), type="vote"))
train_votes <- data.frame(predict(model5, train5, type="vote"))
oob_votes <- data.frame(model5$votes)
row.names(votes) <- strain[1,1]
write.csv(votes, file="votes.csv")
