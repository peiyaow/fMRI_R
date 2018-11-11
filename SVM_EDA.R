library(caret)
library(huge)
library(NetworkToolbox)
library(e1071)

library(R.matlab)
library(abind)
library(reticulate)

data = readMat("/Users/MonicaW/Documents/Research/graph_matlab/ADNI/AD_array.mat")$timeseries.AD
label = read.table("/Users/MonicaW/Documents/Research/graph_matlab/ADNI/label.txt")$V1
ix = read.table("/Users/MonicaW/Documents/Research/graph_matlab/ADNI/ix.txt")$V1

data.list = lapply(1:4, function(l) data[label==l,,])
data.list = lapply(data.list, function(x) sapply(1:dim(x)[1], function(i) scale(x[i,,])))
data.array.list = lapply(data.list, function(x) aperm(array_reshape(x, dim = c(116,137,dim(x)[2])), c(3,2,1)))

n.vec = sapply(1:4, function(l) dim(data.array.list[[l]])[1])
train_val_test.list = lapply(1:4, function(l) createFolds(1:n.vec[l], k = 3))
graph.train.list = list()
for (l in 1:2){
  graph.train.list[[l]] = lapply(train_val_test.list[[l]][[1]], function(train_ix) huge(data.array.list[[l]][train_ix,,], lambda.min.ratio = 0.05))
}
# lambda_ix = 3
# graph.train.list[[l]][[1]]$path[[lambda_ix]]
graph.val.list = list()
for (l in 1:2){
  graph.val.list[[l]] = lapply(train_val_test.list[[l]][[2]], function(train_ix) huge(data.array.list[[l]][train_ix,,], lambda.min.ratio = 0.05))
}

n.train.vec = sapply(1:4, function(l) length(train_val_test.list[[l]][[1]]))
n.val.vec = sapply(1:4, function(l) length(train_val_test.list[[l]][[2]]))

feature.train.list = list()
for (lambda_ix in 1:10){
  feature.train.list[[lambda_ix]] = list()
  for (l in 1:2){
    feature.train.list[[lambda_ix]][[l]] = lapply(1:n.train.vec[l], function(i) graph.features(graph.train.list[[l]][[i]]$path[[lambda_ix]]))
  }
}
feature.train.concat.list = lapply(1:10, function(lambda_ix) do.call(rbind, sapply(1:2, function(l) do.call(rbind, feature.train.list[[lambda_ix]][[l]])))) # length:lambda
y.train = as.factor(unlist(sapply(1:2, function(l) rep(l, n.train.vec[l]))))

feature.val.list = list()
for (lambda_ix in 1:10){
  feature.val.list[[lambda_ix]] = list()
  for (l in 1:2){
   feature.val.list[[lambda_ix]][[l]] = lapply(1:n.val.vec[l], function(i) graph.features(graph.val.list[[l]][[i]]$path[[lambda_ix]]))
  }
}
feature.val.concat.list = lapply(1:10, function(lambda_ix) do.call(rbind, sapply(1:2, function(l) do.call(rbind,feature.val.list[[lambda_ix]][[l]])))) # length:lambda
y.val = as.factor(unlist(sapply(1:2, function(l) rep(l, n.val.vec[l]))))

p.list = list()
i = 1
cost.vec = exp(seq(log(.01),log(.2), length.out = 10))# linear
cost.vec = exp(seq(log(.5),log(2), length.out = 10))# radial

for (cost in cost.vec){
  svm.fit.list = lapply(3:10, function(lambda_ix) svm(x = feature.train.concat.list[[lambda_ix]], y = y.train, kernel = "radial", cost = cost))
  pred = lapply(3:10, function(lambda_ix) predict(svm.fit.list[[lambda_ix-2]], feature.val.concat.list[[lambda_ix]]))
  p.list[[i]] = unlist(lapply(3:10, function(lambda_ix) sum(pred[[lambda_ix-2]] == y.val)/length(y.val)))
  i = i+1
}
do.call(rbind, p.list)




