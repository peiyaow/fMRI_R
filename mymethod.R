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
n.vec = as.vector(summary(label))

L = 2
data.list = lapply(1:L, function(l) data[label==l,,])
ix.train.list = lapply(1:L, function(l) unlist(createDataPartition(1:n.vec[l], times = 1, p = 3/4)))
data.train.list = lapply(1:L, function(l) data.list[[l]][ix.train.list[[l]],,])
data.test.list = lapply(1:L, function(l) data.list[[l]][-ix.train.list[[l]],,])
n.train.vec = sapply(ix.train.list, function(x) length(x))
n.test.vec = n.vec[1:L] - n.train.vec
data.concat.train.list = scale_data(data.train.list)
data.concat.test.list = scale_data(data.test.list)
label.train = as.factor(unlist(sapply(1:L, function(l) rep(l, n.train.vec[l]))))
label.test = as.factor(unlist(sapply(1:L, function(l) rep(l, n.test.vec[l]))))

lambda_ix = 10
cl = makeCluster(4) # number of cores you can use
registerDoParallel(cl)
X.train.test.group.list = foreach(col_ix = 1:116, .packages = c("grplasso", "reticulate", "glmnet")) %dopar% {
  getX.group.parallel(data.concat.train.list, data.concat.test.list, lambda_ix, col_ix)
}
stopCluster(cl)

X1.train.test.group.list = lapply(1:116, function(col_ix) X.train.test.group.list[[col_ix]][[1]])
X1.train.test.group.list = sapply(X1.train.test.group.list, function(list) as.array(list))
X1.train = do.call(cbind, X1.train.test.group.list[1,])
X1.test = do.call(cbind, X1.train.test.group.list[2,])

X2.train.test.group.list = lapply(1:116, function(col_ix) X.train.test.group.list[[col_ix]][[2]])
X2.train.test.group.list = sapply(X2.train.test.group.list, function(list) as.array(list))
X2.train = do.call(cbind, X2.train.test.group.list[1,])
X2.test = do.call(cbind, X2.train.test.group.list[2,])

lambda.vec = exp(seq(log(0.05), log(0.0005), length.out = 100))
ml1 = cv.logistic(X1.train, X2.train, label.train, 10, lambda.vec, 0.2)
logistic.list1 = glmnet(x = X1.train, y = label.train, family = "binomial", standardize = F, alpha = 0.2, lambda = lambda.vec)
logistic.list2 = glmnet(x = X2.train, y = label.train, family = "binomial", standardize = F, alpha = 0.2, lambda = lambda.vec)
prob1 = predict(logistic.list1, s = lambda.vec[ml1[[1]]], newx = X1.test, type = "response")
prob2 = predict(logistic.list2, s = lambda.vec[ml1[[1]]], newx = X2.test, type = "response")
prob = (prob1 + prob2)/2
pred1 = prob2pred(prob1)
pred2 = prob2pred(prob2)
pred = prob2pred(prob)
pred.mtx = cbind(pred1, pred2, pred, label.test)
prob.mtx = cbind(prob1, prob2, prob, label.test)
row.names(prob.mtx) = seq(1, length(label.test))
acc.ml1 = sum(pred == label.test)/length(label.test)

ml0 = cv.logistic0(X1.train, X2.train, label.train, 10, lambda.vec, 0.2)
logistic.list1 = glmnet(x = X1.train, y = label.train, family = "binomial", standardize = F, alpha = 0.2, lambda = lambda.vec)
logistic.list2 = glmnet(x = X2.train, y = label.train, family = "binomial", standardize = F, alpha = 0.2, lambda = lambda.vec)
prob1 = predict(logistic.list1, s = lambda.vec[ml0[[1]]], newx = X1.test, type = "response")
prob2 = predict(logistic.list2, s = lambda.vec[ml0[[2]]], newx = X2.test, type = "response")
prob = (prob1 + prob2)/2
pred1 = prob2pred(prob1)
pred2 = prob2pred(prob2)
pred = prob2pred(prob)
pred.mtx = cbind(pred1, pred2, pred, label.test)
prob.mtx = cbind(prob1, prob2, prob, label.test)
row.names(prob.mtx) = seq(1, length(label.test))
acc.ml0 = sum(pred == label.test)/length(label.test)

logistic.list1 = cv.glmnet(x = X1.train, y = label.train, family = "binomial", standardize = F, alpha = 0.2, lambda = lambda.vec)
logistic.list2 = cv.glmnet(x = X2.train, y = label.train, family = "binomial", standardize = F, alpha = 0.2, lambda = lambda.vec)
prob1 = predict(logistic.list1, s = logistic.list1$lambda.min, newx = X1.test, type = "response")
prob2 = predict(logistic.list2, s = logistic.list2$lambda.min, newx = X2.test, type = "response")
prob = (prob1 + prob2)/2
pred1 = prob2pred(prob1)
pred2 = prob2pred(prob2)
pred = prob2pred(prob)
pred.mtx = cbind(pred1, pred2, pred, label.test)
prob.mtx = cbind(prob1, prob2, prob, label.test)
row.names(prob.mtx) = seq(1, length(label.test))
acc.ml = sum(pred == label.test)/length(label.test)

file.name = 'accuracy.csv'
write.table(t(c(acc.ml1, acc.ml0, acc.ml)), file = file.name, sep = ',', append = T, col.names = F, row.names = F)










