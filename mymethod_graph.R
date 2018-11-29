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
# source('~/Documents/GitHub/fMRI_R/myfunctions.R')
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
n.train = sum(n.train.vec)
n.test = sum(n.test.vec)
data.concat.train.list = scale_data(data.train.list)
data.concat.test.list = scale_data(data.test.list)
label.train = as.factor(unlist(sapply(1:L, function(l) rep(l, n.train.vec[l]))))
label.test = as.factor(unlist(sapply(1:L, function(l) rep(l, n.test.vec[l]))))

lambda_group_ix = 10
lambda_ix = 10

cl = makeCluster(4) # number of cores you can use
registerDoParallel(cl)
X.group.list = foreach(col_ix = 1:116, .packages = c("grplasso", "reticulate", "glmnet", "SGL")) %dopar% {
  getX.group(data.concat.train.list, col_ix)
} # col_num by group by train and test
stopCluster(cl)

X.group.list = lapply(1:116, function(col_ix) X.group.list[[col_ix]][,,lambda_group_ix]) 
X1.group.list = lapply(1:116, function(col_ix) X.group.list[[col_ix]][, 1]) 
X2.group.list = lapply(1:116, function(col_ix) X.group.list[[col_ix]][, 2])

G1 = mtx2graph(do.call(cbind, X1.group.list))
G2 = mtx2graph(do.call(cbind, X2.group.list))

cl = makeCluster(4) # number of cores you can use
registerDoParallel(cl)
X.train.test.res.list = foreach(col_ix = 1:116, .packages = c("grplasso", "reticulate", "glmnet", "SGL")) %dopar% {
  getX.res.parallel(data.concat.train.list, data.concat.test.list, X.group.list, lambda_ix, col_ix, library = "SGL", alpha = 0.8)
} # col_num by group by train and test
stopCluster(cl)

X1.train.list = lapply(1:116, function(col_ix) X.train.test.res.list[[col_ix]][[1]]$X.train) # 116 by n.train by 115
X1.train.array = array(unlist(X1.train.list), dim = c(n.train, 115, 116))
X2.train.list = lapply(1:116, function(col_ix) X.train.test.res.list[[col_ix]][[2]]$X.train) # 116 by n.train by 115
X2.train.array = array(unlist(X2.train.list), dim = c(n.train, 115, 116))

X1.test.list = lapply(1:116, function(col_ix) X.train.test.res.list[[col_ix]][[1]]$X.test) # 116 by n.train by 115
X1.test.array = array(unlist(X1.test.list), dim = c(n.test, 115, 116))
X2.test.list = lapply(1:116, function(col_ix) X.train.test.res.list[[col_ix]][[2]]$X.test) # 116 by n.train by 115
X2.test.array = array(unlist(X2.test.list), dim = c(n.test, 115, 116))

G1.train = sapply(1:n.train, function(ix) mtx2graph(X1.train.array[ix,,]), simplify = "array")
G2.train = sapply(1:n.train, function(ix) mtx2graph(X2.train.array[ix,,]), simplify = "array")
G1.test = sapply(1:n.test, function(ix) mtx2graph(X1.test.array[ix,,]), simplify = "array")
G2.test = sapply(1:n.test, function(ix) mtx2graph(X2.test.array[ix,,]), simplify = "array")

# # using t-test to select features
# feature.ix.mtx1 = mtx.feature.ix(G1.train, label.train)$TF.mtx
# feature.ix.mtx2 = mtx.feature.ix(G2.train, label.train)$TF.mtx
# 
# X1.train = t(sapply(1:n.train, function(ix) G1.train[,,ix][feature.ix.mtx1]))
# X1.test = t(sapply(1:n.test, function(ix) G1.test[,,ix][feature.ix.mtx1]))
# 
# X2.train = t(sapply(1:n.train, function(ix) G2.train[,,ix][feature.ix.mtx2]))
# X2.test = t(sapply(1:n.test, function(ix) G2.test[,,ix][feature.ix.mtx2]))

X1.train = t(array_reshape(G1.train, c(116*116, n.train)))
X1.test = t(array_reshape(G1.test, c(116*116, n.test)))

X2.train = t(array_reshape(G2.train, c(116*116, n.train)))
X2.test = t(array_reshape(G2.test, c(116*116, n.test)))

X.train = cbind(X1.train, X2.train)
X.test = cbind(X1.test, X2.test)

library(randomForest)
rf = randomForest(X.train, label.train)
sum(predict(rf, X.test) == label.test)/n.test

# try out lambda
logistic.list1 = cv.glmnet(x = X.train, y = label.train, family = "binomial", standardize = F, alpha = 1)
sum(predict(logistic.list1, s = logistic.list1$lambda.min, newx = X.test, type = "class") == label.test)/n.test
logistic.list1$lambda
#

lambda.vec = exp(seq(log(200), log(2), length.out = 100))
lambda.vec = exp(seq(log(0.05), log(0.05*1e-2), length.out = 100))

# separately train cv 
# 0 two lambdas
# 1 one lambda
# 2 one lambda and w
set.seed(10)
logistic.list1 = cv.glmnet(x = X1.train, y = label.train, family = "binomial", standardize = T, alpha = 0, lambda = lambda.vec)
logistic.list2 = cv.glmnet(x = X2.train, y = label.train, family = "binomial", standardize = T, alpha = 0, lambda = lambda.vec)
prob1 = predict(logistic.list1, s = logistic.list1$lambda.min, newx = X1.test, type = "response")
prob2 = predict(logistic.list2, s = logistic.list2$lambda.min, newx = X2.test, type = "response")
prob = (prob1 + prob2)/2
# pred1 = prob2pred(prob1)
# pred2 = prob2pred(prob2)
pred = prob2pred(prob)
# pred.mtx = cbind(pred1, pred2, pred, label.test)
# prob.mtx = cbind(prob1, prob2, prob, label.test)
# row.names(prob.mtx) = seq(1, length(label.test))
acc.ml = sum(pred == label.test)/length(label.test)

set.seed(10)
ml0 = cv.logistic0(X1.train, X2.train, label.train, 10, lambda.vec, 0)
logistic.list1 = glmnet(x = X1.train, y = label.train, family = "binomial", standardize = T, alpha = 0, lambda = lambda.vec)
logistic.list2 = glmnet(x = X2.train, y = label.train, family = "binomial", standardize = T, alpha = 0, lambda = lambda.vec)
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

set.seed(10)
ml1 = cv.logistic1(X1.train, X2.train, label.train, 10, lambda.vec, 0, F)
logistic.list1 = glmnet(x = X1.train, y = label.train, family = "binomial", standardize = F, alpha = 0, lambda = lambda.vec)
logistic.list2 = glmnet(x = X2.train, y = label.train, family = "binomial", standardize = F, alpha = 0, lambda = lambda.vec)
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

set.seed(10)
ml2 = cv.logistic2(X1.train, X2.train, label.train, 10, lambda.vec, 0, T)
logistic.list1 = glmnet(x = X1.train, y = label.train, family = "binomial", standardize = T, alpha = 0, lambda = lambda.vec)
logistic.list2 = glmnet(x = X2.train, y = label.train, family = "binomial", standardize = T, alpha = 0, lambda = lambda.vec)
prob1 = predict(logistic.list1, s = lambda.vec[ml2[[1]]], newx = X1.test, type = "response")
prob2 = predict(logistic.list2, s = lambda.vec[ml2[[1]]], newx = X2.test, type = "response")

prob1 = predict(logistic.list1, s = lambda.vec[100], newx = X1.test, type = "response")
prob2 = predict(logistic.list2, s = lambda.vec[100], newx = X2.test, type = "response")
w = seq(0,1,length.out = 11)
prob = prob1*w[ml2[[2]]] + prob2*(1-w[ml2[[2]]])
pred1 = prob2pred(prob1)
pred2 = prob2pred(prob2)
pred = prob2pred(prob)
pred.mtx = cbind(pred1, pred2, pred, label.test)
prob.mtx = cbind(prob1, prob2, prob, label.test)
row.names(prob.mtx) = seq(1, length(label.test))
acc.ml2 = sum(pred == label.test)/length(label.test)



file.name = 'accuracy.csv'
write.table(t(c(acc.ml0, acc.ml1, acc.ml2, acc.ml)), file = file.name, sep = ',', append = T, col.names = F, row.names = F)










