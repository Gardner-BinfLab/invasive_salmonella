library("randomForest")
library("DMwR")

load("finalmodel.Rdata")
lookup <- read.delim("lookup.tsv")
strain <- read.table("strain", comment.char="@")

qc <- read.delim("NsinAssembly", header=F)
qcnum <- sum(qc[,1]>0)

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

# impute missing scores based on nearest training example (one per major serovar)
merge <- rbind(testdata, train5[,-197])
merge2 <- knnImputation(merge, k=1)
testdata <- merge2[1,]

# load the random forest model and score your isolates
votes <- data.frame(predict(model5, na.roughfix(testdata), type="vote"))
train_votes <- data.frame(predict(model5, train5, type="vote"))
oob_votes <- data.frame(model5$votes)
row.names(votes) <- strain[1,1]
write.csv(cbind(votes, qcnum), file="votes.csv", col.names = F)
