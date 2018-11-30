# ---------------------- reading shell command --------------------- 
args = (commandArgs(TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}
# ------------------------------------------------------------------ 

library(R.matlab)
library(reshape2)
library(abind)
library(glmnet)
library(grplasso)
library(caret)
library(doParallel)
library(glasso)
library(SGL)
library(e1071)

set.seed(myseed)

# mac
# source('~/Documents/Research/coding/R/fMRI_R/myfunctions.R')
# data = readMat("/Users/MonicaW/Documents/Research/graph_matlab/ADNI/AD_array.mat")$timeseries.AD
# ix = read.table("/Users/MonicaW/Documents/Research/graph_matlab/ADNI/ix.txt")$V1
# label = read.table("/Users/MonicaW/Documents/Research/graph_matlab/ADNI/label.txt")$V1

# longleaf
source('~/fMRI_R/myfunctions.R')
data = readMat("~/fMRI_data/AD_array.mat")$timeseries.AD
ix = read.table("~/fMRI_data/ix.txt")$V1
label = read.table("~/fMRI_data/label.txt")$V1

unique.ix = sapply(1:174, function(i) seq(1,563)[ix == i][1])
data = data[unique.ix,,]
label = as.factor(label[unique.ix])

# subset.ix = seq(1, 174, by = 5)
# data = data[subset.ix,,]
# label = label[subset.ix]

n.vec = as.vector(summary(label))
L = 2
data.list = lapply(1:L, function(l) data[label==l,,])
data.list = lapply(data.list, function(x) sapply(1:dim(x)[1], function(i) scale(x[i,,])))
data.list = lapply(data.list, function(x) aperm(array_reshape(x, dim = c(116,137,dim(x)[2])), c(3,2,1)))

data_array = abind(data.list[[1]], data.list[[2]], along = 1)

cl = makeCluster(4) # number of cores you can use
registerDoParallel(cl)
G.mtx.list = foreach(col_ix = 1:116, .packages = c("SGL", "reticulate", "glasso", "glmnet")) %dopar% {
  getGraph2.parallel(data_array, col_ix)
}
stopCluster(cl)


G.mtx.array = sapply(G.mtx.list, function(x) x, simplify = "array")
graphs = aperm(G.mtx.array, c(1,3,2))

# symmetric
graphs = array_reshape(apply(graphs, 3, function(x) {
  x = (x+t(x))/2
  x[x<0] = 0
  return(x)
}), c(116, 116, dim(graphs)[3]))

X.feature = t(apply(graphs, 3, graph.features))
Y.label = as.factor(c(rep(1, n.vec[1]), rep(2, n.vec[2])))
feature.ix = seq(1,116)[apply(X.feature, 2, sd) != 0]
X.feature = X.feature[, feature.ix]

N.vec = c(0, cumsum(n.vec))
ix.train.list = lapply(1:L, function(l) ((N.vec[l]+1):N.vec[l+1])[unlist(createDataPartition(1:n.vec[l], times = 1, p = 3/4))])
ix.train = unlist(ix.train.list)
X.feature.train = X.feature[ix.train,]
X.feature.test = X.feature[-ix.train,]
Y.label.train = Y.label[ix.train]
Y.label.test = Y.label[-ix.train]
n.train = length(Y.label.train)
n.test = length(Y.label.test)

cost.vec = exp(seq(log(.1), log(.01), length.out = 100))
#cost.vec = seq(.1,.01, length.out = 10)
svm.list = lapply(cost.vec, function(cost) svm(x = X.feature.train, y = Y.label.train, scale = T, kernel = "linear", cost = cost))
p.table = sapply(1:100, function(ix) sum(predict(svm.list[[ix]], X.feature.test) == Y.label.test)/n.test)
nSV.table = sapply(1:100, function(ix) svm.list[[ix]]$nSV)
print(nSV.table)
svm.cv.ml = cv.svm(X.feature.train, Y.label.train, 10, cost.vec)
acc = sum(predict(svm.cv.ml[[4]], X.feature.test) == Y.label.test)/n.test

file.name = 'accuracy.csv'
write.table(t(c(p.table, acc)), file = file.name, sep = ',', append = T, col.names = F, row.names = F)




