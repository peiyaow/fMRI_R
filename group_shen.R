library(R.matlab)
library(SGL)
#library(Matrix)
library(doParallel)
library(grplasso)
library(e1071)
source('~/Documents/Research/coding/R/fMRI_R/myfunctions.R')
data = readMat("/Users/MonicaW/Documents/Research/graph_matlab/ADNI/AD_array.mat")$timeseries.AD
ix = read.table("/Users/MonicaW/Documents/Research/graph_matlab/ADNI/ix.txt")$V1
label = read.table("/Users/MonicaW/Documents/Research/graph_matlab/ADNI/label.txt")$V1

unique.ix = sapply(1:174, function(i) seq(1,563)[ix == i][1])
data = data[unique.ix,,]
label = label[unique.ix]

subset.ix = seq(1, 174, by = 3)
data = data[subset.ix,,]
label = label[subset.ix]

n.vec = as.vector(summary(as.factor(label)))
data.list = lapply(1:4, function(l) data[label==l,,])
data.concat.list = scale_data(data.list)
data.concat = do.call(rbind, data.concat.list[1:2])

cl = makeCluster(4) # number of cores you can use
registerDoParallel(cl)
G.mtx.list = foreach(col_ix = 1:116, .packages = c("SGL", "reticulate", "grplasso")) %dopar% {
  getGraph.parallel(data.concat, col_ix)
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
X.feature = X.feature[,feature.ix]
cost.vec = exp(seq(log(.01),log(.03), length.out = 10))
svm.list = lapply(cost.vec, function(cost) svm(x = X.feature, y = Y.label, scale = T, kernel = "linear", cost = cost))










