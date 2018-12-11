library(reticulate)
library(SGL)
library(grplasso)
library(glasso)

scale_data = function(data.list){
  data.list = lapply(data.list, function(x) sapply(1:dim(x)[1], function(i) scale(x[i,,])))
  data.array.list = lapply(data.list, function(x) aperm(array_reshape(x, dim = c(116,137,dim(x)[2])), c(3,2,1)))
  data.concat.list = lapply(data.list, function(x) array_reshape(aperm(array_reshape(x, dim = c(116,137,dim(x)[2])), c(3,2,1)), dim = c(dim(x)[2]*137, 116), order = "C"))
  return(data.concat.list)
}

prob2pred = function(prob){
  pred = rep(1,length(prob))
  pred[prob>0.5] = 2
  return(pred)
}

add.diagonal = function(Beta, col_ix, have.intercept = F){
  if (have.intercept){
    Beta = Beta[-1,]
  }
  p = nrow(Beta)+1
#  print(p)
  n = ncol(Beta)
  newBeta = matrix(0, ncol = n, nrow = p)
  if (col_ix == 1){
    newBeta[col_ix,] = 1
    newBeta[(col_ix+1):p,] = Beta[col_ix:(p-1),]
  }else if (col_ix == p){
    newBeta[col_ix,] = 1
    newBeta[1:(p-1),] = Beta
  }else{
    newBeta[1:(col_ix-1),] = Beta[1:(col_ix-1),]
    newBeta[col_ix,] = 1
    newBeta[(col_ix+1):p,] = Beta[col_ix:(p-1),]
  }
  return(newBeta)
}

mtx2graph = function(mtx){
  p = ncol(mtx)
  graph = sapply(1:p, function(col_ix) add.diagonal(as.matrix(mtx[, col_ix]), col_ix))
  return(graph)
}

graph.features <- function(A){
  A = as.matrix(A)
  return(NetworkToolbox::clustcoeff(A, weighted = T)[[2]])
}

getX = function(data.concat.train.list, data.concat.test.list, l_test, lambda_ix, col_ix){
  fit.list = lapply(1:4, function(l) glmnet(x = data.concat.train.list[[l]][,-col_ix], y = data.concat.train.list[[l]][,col_ix], standardize = F, nlambda = 10, lambda.min.ratio = 0.25, intercept = F))
  res.list = lapply(1:4, function(l) data.concat.train.list[[l]][,col_ix] - predict(fit.list[[l_test]], newx = data.concat.train.list[[l]][,-col_ix], s = fit.list[[l]]$lambda[lambda_ix])) # true residual list
  
  # given residuals refit every subject
  coef.mtx.list = list()
  for (l in 1:4){
    coef.list = list()
    for (subject_ix in 1:n.train.vec[l]){
      res_fit = glmnet(x = data.concat.train.list[[l]][137*(subject_ix-1)+1:137, -col_ix], y = res.list[[l]][137*(subject_ix-1)+1:137], standardize = F, nlambda = 10, lambda.min.ratio = 0.25, intercept = F)
      coef.list[[subject_ix]] = coef(res_fit, s=res_fit$lambda[lambda_ix])
    }
    coef.mtx = do.call(cbind, coef.list)
    coef.mtx.list[[l]] = coef.mtx
  }
  
  X.train = as.matrix(t(do.call(cbind, coef.mtx.list[1:2])))
  
  res4.mtx.list = list()
  for (l_test in 1:4){
    res4.list = lapply(1:4, function(l) data.concat.test.list[[l_test]][,col_ix] - predict(fit.list[[l]], newx = data.concat.test.list[[l_test]][,-col_ix], s = fit.list[[l]]$lambda[lambda_ix])) # candidate residual list
    res4.mtx.list[[l_test]] = do.call(cbind, res4.list)
  }
  
  coef.test.mtx.list = list()
  for (l_test in 1:4){
    coef.test.list = list()
    for (subject_ix in 1:n.test.vec[[l_test]]){
      res_fit = lapply(1:4, function(l) glmnet(x = data.concat.test.list[[l_test]][137*(subject_ix-1)+1:137, -col_ix], y = res4.mtx.list[[l_test]][137*(subject_ix-1)+1:137, l], standardize = F, nlambda = 10, lambda.min.ratio = 0.25, intercept = F))
      coef.test.list[[subject_ix]] = lapply(1:4, function(l) coef(res_fit[[l]], s=res_fit[[l]]$lambda[lambda_ix]))
      coef.test.list[[subject_ix]] = do.call(cbind, coef.test.list[[subject_ix]])
    }
    coef.test.mtx.list[[l_test]] = coef.test.list
  }
  array_test_mtx.list = lapply(1:4, function(l_test) aperm(sapply(coef.test.mtx.list[[l_test]], function(X) as.matrix(X), simplify = "array"), c(2,3,1)))
  X.test = do.call(rbind, lapply(1:2, function(l) array_test_mtx.list[[l]][l_test,,]))
  return(list(X.train, X.test))
}

getX.group = function(data.concat.train.list, data.concat.test.list, l_test, lambda_ix, col_ix){
  Y.train = do.call(c, lapply(data.concat.train.list[1:L], function(x) x[,col_ix]))
  x.list = lapply(data.concat.train.list[1:L], function(x) x[,-col_ix])
  X.train.group = matrix(0, nrow = length(Y.train), ncol = 115*L)
  for (l in 1:L){
    X.train.group[(vec[l]+1):vec[l+1], 115*(l-1)+1:115] = x.list[[l]]
  }
  group_ix = rep(seq(1,115), L)
  lambda_max = lambdamax(x = X.train.group, y = Y.train, index = group_ix, model = LinReg(), center = F, standardize = F)
  lambdas = exp(seq(log(lambda_max), log(lambda_max*0.25), length.out = 10))
  group_fit = grplasso(x = X.train.group, y = Y.train, lambda = lambdas, index = group_ix, model = LinReg(), center = F, standardize = F)
  #  group_fit = grplasso(x = X.train.group, y = Y.train, lambda = lambda_max*0.25, index = group_ix, model = LinReg(), center = F, standardize = F)
  group_coef = array_reshape(group_fit$coefficients, c(115, L, 10), order = "F")
  #  group_coef = array_reshape(group_fit$coefficients, c(115, 4, 1), order = "F")
  res.list = lapply(1:L, function(l) data.concat.train.list[[l]][,col_ix] - data.concat.train.list[[l]][,-col_ix]%*%group_coef[,,lambda_ix][,l_test])
  
  # given residuals refit every subject
  coef.mtx.list = list()
  for (l in 1:L){
    coef.list = list()
    for (subject_ix in 1:n.train.vec[l]){
      res_fit = glmnet(x = data.concat.train.list[[l]][137*(subject_ix-1)+1:137, -col_ix], y = res.list[[l]][137*(subject_ix-1)+1:137], standardize = F, nlambda = 10, lambda.min.ratio = 0.25, intercept = F)
      coef.list[[subject_ix]] = coef(res_fit, s=res_fit$lambda[lambda_ix])
    }
    coef.mtx = do.call(cbind, coef.list)
    coef.mtx.list[[l]] = coef.mtx
  }
  X.train = as.matrix(t(do.call(cbind, coef.mtx.list)))
  
  # testing...
  res4.mtx.list = list()
  for (l_test in 1:L){
    #res4.list = lapply(1:4, function(l) data.concat.test.list[[l_test]][,col_ix] - predict(fit.list[[l]], newx = data.concat.test.list[[l_test]][,-col_ix], s = fit.list[[l]]$lambda[lambda_ix])) # candidate residual list
    res4.list = lapply(1:L, function(l) data.concat.test.list[[l_test]][,col_ix] - data.concat.test.list[[l_test]][,-col_ix]%*%group_coef[,,lambda_ix][,l]) # candidate residual list
    res4.mtx.list[[l_test]] = do.call(cbind, res4.list)
  }
  
  coef.test.mtx.list = list()
  for (l_test in 1:L){
    coef.test.list = list()
    for (subject_ix in 1:n.test.vec[[l_test]]){
      res_fit = lapply(1:L, function(l) glmnet(x = data.concat.test.list[[l_test]][137*(subject_ix-1)+1:137, -col_ix], y = res4.mtx.list[[l_test]][137*(subject_ix-1)+1:137, l], standardize = F, nlambda = 10, lambda.min.ratio = 0.25, intercept = F))
      coef.test.list[[subject_ix]] = lapply(1:L, function(l) coef(res_fit[[l]], s=res_fit[[l]]$lambda[lambda_ix]))
      coef.test.list[[subject_ix]] = do.call(cbind, coef.test.list[[subject_ix]])
    }
    coef.test.mtx.list[[l_test]] = coef.test.list
  }
  array_test_mtx.list = lapply(1:L, function(l_test) aperm(sapply(coef.test.mtx.list[[l_test]], function(X) as.matrix(X), simplify = "array"), c(2,3,1)))
  X.test = do.call(rbind, lapply(1:L, function(l) array_test_mtx.list[[l]][l_test,,]))
  return(list(X.train, X.test))
}

# getX.group.parallel = function(data.concat.train.list, data.concat.test.list, lambda_ix, col_ix, library = "SGL"){
#   L = length(data.concat.train.list)
#   Y.train = do.call(c, lapply(data.concat.train.list, function(x) x[,col_ix]))
#   X.list = lapply(data.concat.train.list, function(x) x[,-col_ix])
#   X.train.group = matrix(0, nrow = length(Y.train), ncol = 115*L)
#   vec = c(0, cumsum(sapply(X.list, nrow)))
#   for (l in 1:L){
#     X.train.group[(vec[l]+1):vec[l+1], 115*(l-1)+1:115] = X.list[[l]]
#   }
#   group_ix = rep(seq(1,115), L)
#   
#   if (library == "SGL"){
#     # use SGL
#     SGLdata = list(x = X.train.group, y = Y.train)
#     SGL.list = SGL(data = SGLdata, index = group_ix, min.frac = 0.25, nlam = 10, alpha = 0.8, standardize = F)
#     #beta.list = lapply(1:10, function(lam_ix) add.diagonal(matrix(SGL.list$beta[,lam_ix], ncol = L), col_ix))
#     group_coef = array(SGL.list$beta, dim = c(115, L, 10))
#   }else if (library == "grplasso"){
#     # use grplasso
#     lambda_max = lambdamax(x = X.train.group, y = Y.train, index = group_ix, model = LinReg(), center = F, standardize = F)
#     lambdas = exp(seq(log(lambda_max), log(lambda_max*0.25), length.out = 10))
#     group_fit = grplasso(x = X.train.group, y = Y.train, lambda = lambdas, index = group_ix, model = LinReg(), center = F, standardize = F)
#     group_coef = array_reshape(group_fit$coefficients, c(115, L, 10), order = "F")
#   }
#   
#   X.list = list()
#   for (l_test in 1:L){
#     res.list = lapply(1:L, function(l) data.concat.train.list[[l]][,col_ix] - data.concat.train.list[[l]][,-col_ix]%*%group_coef[,,lambda_ix][,l_test])
#     # given residuals refit every subject
#     coef.mtx.list = list()
#     for (l in 1:L){
#       coef.list = list()
#       for (subject_ix in 1:n.train.vec[l]){
#         res_fit = glmnet(x = data.concat.train.list[[l]][137*(subject_ix-1)+1:137, -col_ix], y = res.list[[l]][137*(subject_ix-1)+1:137], standardize = F, nlambda = 10, lambda.min.ratio = 0.25, intercept = F)
#         coef.list[[subject_ix]] = coef(res_fit, s=res_fit$lambda[lambda_ix])
#       }
#       coef.mtx = do.call(cbind, coef.list)
#       coef.mtx.list[[l]] = coef.mtx
#     }
#     X.train = as.matrix(t(do.call(cbind, coef.mtx.list)))
#     
#     res.test.list = lapply(1:L, function(l) data.concat.test.list[[l]][,col_ix] - data.concat.test.list[[l]][,-col_ix]%*%group_coef[,,lambda_ix][,l_test])
#     # given residuals refit every subject
#     coef.mtx.test.list = list()
#     for (l in 1:L){
#       coef.list = list()
#       for (subject_ix in 1:n.test.vec[l]){
#         res_fit = glmnet(x = data.concat.test.list[[l]][137*(subject_ix-1)+1:137, -col_ix], y = res.test.list[[l]][137*(subject_ix-1)+1:137], standardize = F, nlambda = 10, lambda.min.ratio = 0.25, intercept = F)
#         coef.list[[subject_ix]] = coef(res_fit, s=res_fit$lambda[lambda_ix])
#       }
#       coef.mtx = do.call(cbind, coef.list)
#       coef.mtx.test.list[[l]] = coef.mtx
#     }
#     X.test = as.matrix(t(do.call(cbind, coef.mtx.test.list)))
#     X.list[[l_test]] = list(X.train = X.train, X.test = X.test)
#   }
#   return(X.list)
# }
# different lambda
cv.logistic0 = function(X1, X2, label, nfolds, lambda.vec, alpha, standardize = T){
  flds = createFolds(label, k = nfolds, list = TRUE, returnTrain = FALSE)
  len_lam = length(lambda.vec)
  acc.mtx.list = list()
  for (k in 1:nfolds){
    X1.train = X1[unlist(flds[-k]), ]
    X1.val = X1[unlist(flds[k]), ]
    X2.train = X2[unlist(flds[-k]), ]
    X2.val = X2[unlist(flds[k]), ]
    label.train = label[unlist(flds[-k])]
    label.val = label[unlist(flds[k])]
    logistic.list1 = glmnet(x = X1.train, y = label.train, family = "binomial", standardize = standardize, alpha = alpha, lambda = lambda.vec)
    logistic.list2 = glmnet(x = X2.train, y = label.train, family = "binomial", standardize = standardize, alpha = alpha, lambda = lambda.vec)
    prob1 = predict(logistic.list1, newx = X1.val, type = "response") # n.val by n.lambda
    prob2 = predict(logistic.list2, newx = X2.val, type = "response")
    prob.list = lapply(1:len_lam, function(i) sapply(1:len_lam, function(j) (prob1[,i]+prob2[,j])/2))
    pred.list = lapply(prob.list, function(prob.mtx) apply(prob.mtx, 2, prob2pred))
    acc.list = lapply(pred.list, function(pred.mtx) apply(pred.mtx, 2, function(pred.vec) sum(pred.vec== label.val)/length(label.val)))
    acc.mtx.list[[k]] = do.call(rbind, acc.list)
  }
  acc.array = array(unlist(acc.mtx.list), dim = c(len_lam, len_lam, nfolds))
  acc.mtx = apply(acc.array, c(1,2), function(x) mean(x))
  sd.mtx = apply(acc.array, c(1,2), function(x) sd(x))
  TF.mtx = acc.mtx == max(acc.mtx)
  sd.vec = sd.mtx[TF.mtx]
  min.sd = min(sd.vec)
  TF.mtx = TF.mtx*(sd.mtx==min.sd)
  id.mtx.max = expand.grid(seq(1,len_lam), seq(1,len_lam))[as.logical(as.vector(TF.mtx)),]
  max.id = id.mtx.max[which.max(apply(id.mtx.max, 1, sum)), ]
  #max.id = which.max(acc.mtx)
  #max.id1 = (max.id - 1)%%len_lam + 1
  #max.id2 = ceiling(max.id/len_lam)
  return(list(max.id[[1]], max.id[[2]], sd.vec, id.mtx.max, acc.mtx, acc.mtx.list))
}

# same lambda
cv.logistic1 = function(X1, X2, label, nfolds, lambda.vec, alpha, standardize = T){
  flds = createFolds(label, k = nfolds, list = TRUE, returnTrain = FALSE)
  len_lam = length(lambda.vec)
  acc.mtx.list = list()
  for (k in 1:nfolds){
    X1.train = X1[unlist(flds[-k]), ]
    X1.val = X1[unlist(flds[k]), ]
    X2.train = X2[unlist(flds[-k]), ]
    X2.val = X2[unlist(flds[k]), ]
    label.train = label[unlist(flds[-k])]
    label.val = label[unlist(flds[k])]
    logistic.list1 = glmnet(x = X1.train, y = label.train, family = "binomial", standardize = standardize, alpha = alpha, lambda = lambda.vec)
    logistic.list2 = glmnet(x = X2.train, y = label.train, family = "binomial", standardize = standardize, alpha = alpha, lambda = lambda.vec)
    prob1 = predict(logistic.list1, newx = X1.val, type = "response") # n.val by n.lambda
    prob2 = predict(logistic.list2, newx = X2.val, type = "response")
    
    prob.mtx = (prob1+prob2)/2
    pred.mtx = apply(prob.mtx, 2, prob2pred)
    acc.vec = apply(pred.mtx, 2, function(pred.vec) sum(pred.vec==label.val)/length(label.val))
    acc.mtx.list[[k]] = acc.vec
    # prob.list = lapply(1:len_lam, function(i) sapply(1:len_lam, function(j) (prob1[,i]+prob2[,j])/2))
    # pred.list = lapply(prob.list, function(prob.mtx) apply(prob.mtx, 2, prob2pred))
    # acc.list = lapply(pred.list, function(pred.mtx) apply(pred.mtx, 2, function(pred.vec) sum(pred.vec== label.val)/length(label.val)))
    # acc.mtx.list[[k]] = do.call(rbind, acc.list)
  }
  # acc.array = array(unlist(acc.mtx.list), dim = c(len_lam, len_lam, nfolds))
  # acc.mtx = apply(acc.array, c(1,2), function(x) mean(x))
  # sd.mtx = apply(acc.array, c(1,2), function(x) sd(x))
  # TF.mtx = acc.mtx == max(acc.mtx) 
  # sd.vec = sd.mtx[TF.mtx]
  # min.sd = min(sd.vec)
  # TF.mtx = TF.mtx*(sd.mtx==min.sd)  
  # id.mtx.max = expand.grid(seq(1,len_lam), seq(1,len_lam))[as.logical(as.vector(TF.mtx)),]
  # max.id = id.mtx.max[which.max(apply(id.mtx.max, 1, sum)), ]
  # #max.id = which.max(acc.mtx)
  # #max.id1 = (max.id - 1)%%len_lam + 1
  # #max.id2 = ceiling(max.id/len_lam)
  # return(list(max.id[[1]], max.id[[2]], sd.vec, id.mtx.max, acc.mtx, acc.mtx.list))
  acc.mtx = do.call(rbind, acc.mtx.list)
  acc.vec = apply(acc.mtx, 2, mean)
  return(list(which.max(acc.vec), acc.mtx, acc.mtx.list))
}

# same lambda different w
cv.logistic2 = function(X1, X2, label, nfolds, lambda.vec, alpha, standardize = T){
  flds = createFolds(label, k = nfolds, list = TRUE, returnTrain = FALSE)
  len_lam = length(lambda.vec)
  acc.mtx.list = list()
  w = seq(0, 1, length.out = 11)
  for (k in 1:nfolds){
    X1.train = X1[unlist(flds[-k]), ]
    X1.val = X1[unlist(flds[k]), ]
    X2.train = X2[unlist(flds[-k]), ]
    X2.val = X2[unlist(flds[k]), ]
    label.train = label[unlist(flds[-k])]
    label.val = label[unlist(flds[k])]
    logistic.list1 = glmnet(x = X1.train, y = label.train, family = "binomial", standardize = standardize, alpha = alpha, lambda = lambda.vec)
    logistic.list2 = glmnet(x = X2.train, y = label.train, family = "binomial", standardize = standardize, alpha = alpha, lambda = lambda.vec)
    prob1 = predict(logistic.list1, newx = X1.val, type = "response") # n.val by n.lambda
    prob2 = predict(logistic.list2, newx = X2.val, type = "response")
    
    prob.mtx.list = lapply(w, function(w_) prob1*w_ + prob2*(1-w_))
    pred.mtx.list = lapply(prob.mtx.list, function(prob.mtx) apply(prob.mtx, 2, prob2pred))
    acc.mtx = sapply(pred.mtx.list, function(pred.mtx) apply(pred.mtx, 2, function(pred.vec) sum(pred.vec==label.val)/length(label.val)))
    acc.mtx.list[[k]] = acc.mtx
    # prob.list = lapply(1:len_lam, function(i) sapply(1:len_lam, function(j) (prob1[,i]+prob2[,j])/2))
    # pred.list = lapply(prob.list, function(prob.mtx) apply(prob.mtx, 2, prob2pred))
    # acc.list = lapply(pred.list, function(pred.mtx) apply(pred.mtx, 2, function(pred.vec) sum(pred.vec== label.val)/length(label.val)))
    # acc.mtx.list[[k]] = do.call(rbind, acc.list)
  }
  acc.array = array(unlist(acc.mtx.list), dim = c(len_lam, 11, nfolds))
  acc.mtx = apply(acc.array, c(1,2), function(x) mean(x))
  # sd.mtx = apply(acc.array, c(1,2), function(x) sd(x))
  # TF.mtx = acc.mtx == max(acc.mtx) 
  # sd.vec = sd.mtx[TF.mtx]
  # min.sd = min(sd.vec)
  # TF.mtx = TF.mtx*(sd.mtx==min.sd)  
  # id.mtx.max = expand.grid(seq(1,len_lam), seq(1,len_lam))[as.logical(as.vector(TF.mtx)),]
  # max.id = id.mtx.max[which.max(apply(id.mtx.max, 1, sum)), ]
  
  max.id = which.max(acc.mtx)
  max.id1 = (max.id - 1)%%len_lam + 1
  max.id2 = ceiling(max.id/len_lam)
  return(list(max.id1, max.id2, acc.mtx, acc.mtx.list))
  
  # acc.mtx = do.call(rbind, acc.mtx.list)
  # acc.vec = apply(acc.mtx, 2, mean)
  # return(list(which.max(acc.vec), acc.mtx, acc.mtx.list))
}




cv.svm = function(X, label, nfolds, cost.vec){
  flds = createFolds(label, k = nfolds, list = TRUE, returnTrain = FALSE)
  len_cost = length(cost.vec)
  acc.mtx.list = list()
  for (k in 1:nfolds){
    X.train = X[unlist(flds[-k]), ]
    X.val = X[unlist(flds[k]), ]
    label.train = label[unlist(flds[-k])]
    label.val = label[unlist(flds[k])]
    
    svm.list = lapply(cost.vec, function(cost) svm(x = X.train, y = label.train, scale = T, kernel = "linear", cost = cost))
    acc.vec = sapply(1:len_cost, function(ix) sum(predict(svm.list[[ix]], X.val) == label.val)/length(label.val))
    acc.mtx.list[[k]] = acc.vec
  }
  acc.mtx = do.call(rbind, acc.mtx.list)
  acc.vec = apply(acc.mtx, 2, mean)
  ix.max = which.max(acc.vec)
  svm.ml = svm(x = X, y = label, scale = T, kernel = "linear", cost = cost.vec[ix.max])
  return(list(ix.max, acc.mtx, acc.mtx.list, svm.ml))
}

getGraph.parallel = function(data.concat, col_ix, lambda_ix = 10, library = 'grplasso'){
  n = dim(data.concat)[1]/137
  Y = data.concat[, col_ix]
  X = matrix(0, nrow = 137*n, ncol = 115*n)
  for (i in 1:n){
    X[137*(i-1)+1:137, 1:115+(i-1)*115] = data.concat[137*(i-1)+1:137, -col_ix]
  }
  index = rep(seq(1, 115), n)
  if (library == 'SGL'){
    SGLdata = list(x = X, y = Y)
    SGL.list = SGL(data = SGLdata, index = index, min.frac = 0.25, nlam = 10, alpha = 0, standardize = F)
    beta.list = lapply(1:10, function(lam_ix) add.diagonal(matrix(SGL.list$beta[,lam_ix], ncol = n), col_ix))
  }else if(library == 'grplasso'){
    lambda_max = lambdamax(x = X, y = Y, index = index, model = LinReg(), center = F, standardize = F)
    lambdas = exp(seq(log(lambda_max), log(lambda_max*0.25), length.out = 10))
    group_fit = grplasso(x = X, y = Y, lambda = lambdas, index = index, model = LinReg(), center = F, standardize = F)
    group_coef = array_reshape(group_fit$coefficients, c(115, n, 10), order = "F")
    beta.list = lapply(1:10, function(lam_ix) add.diagonal(group_coef[,,lam_ix], col_ix))
  }
  return(beta.list[[lambda_ix]])
}

getGraph2.parallel = function(data_array, col_ix, glasso_ix, ix = 10, lambda_ix = 10, library = 'SGL'){
  n = dim(data_array)[1]
  Sigma_list = list()
  for (sub_ix in 1:n){
    Y = data_array[sub_ix,,col_ix]
    X = data_array[sub_ix,,-col_ix]
#    fit = glmnet(x = X, y = Y, standardize = F, nlambda = 10, alpha = 0, intercept = F)
    fit = glmnet(x = X, y = Y, standardize = F, nlambda = 10, lambda.min.ratio = 0.25, alpha = 1, intercept = F)
    Yhat = predict(fit, newx = X)
    E = Yhat - Y
    Sigma = E[,ix]%*%t(E[,ix])
    Sigma_list[[sub_ix]] = Sigma
  }
  
  # compute C
  Sigma_array = sapply(Sigma_list, function(X) X, simplify = "array")
  Sigma = apply(Sigma_array, c(1,2), mean)
  glassofit = glassopath(Sigma)
  S = glassofit$w[,,glasso_ix]
  C = glassofit$wi[,,glasso_ix]
  C_half = expm::sqrtm(C)
  
  X.group = matrix(0, nrow = n*137, ncol = 115*n)
  Y.group = rep(0, n*137)
  for (sub_ix in 1:n){
    Y = C_half%*%data_array[sub_ix,,col_ix]
    X = C_half%*%data_array[sub_ix,,-col_ix]
    Y.group[(137*(sub_ix-1)+1):(137*sub_ix)] = Y
    X.group[(137*(sub_ix-1)+1):(137*sub_ix), (115*(sub_ix-1)+1):(115*sub_ix)] = X
  }
  group.ix = rep(seq(1,115), n)
  if (library == "SGL"){
    SGLdata = list(x = X.group, y = Y.group)
    SGL.list = SGL(data = SGLdata, index = group.ix, min.frac = 0.25, nlam = 10, alpha = 0.8, standardize = F)
    beta.list = lapply(1:10, function(lam_ix) add.diagonal(matrix(SGL.list$beta[,lam_ix], ncol = n), col_ix))
  }else if (library == "grplasso"){
    lambda_max = lambdamax(x = X.group, y = Y.group, index = group.ix, model = LinReg(), center = F, standardize = F)
    lambdas = exp(seq(log(lambda_max), log(lambda_max*0.25), length.out = 10))
    group_fit = grplasso(x = X.group, y = Y.group, lambda = lambdas, index = group.ix, model = LinReg(), center = F, standardize = F)
    group_coef = array_reshape(group_fit$coefficients, c(115, n, 10), order = "F")
    beta.list = lapply(1:10, function(lam_ix) add.diagonal(group_coef[,,lam_ix], col_ix))
  }
  return(beta.list[[lambda_ix]])
}


# get group coef from training data
getX.group = function(data.concat.train.list, col_ix, library = "SGL", alpha = 0.8){
  L = length(data.concat.train.list)
  Y.train = do.call(c, lapply(data.concat.train.list, function(x) x[,col_ix]))
  X.list = lapply(data.concat.train.list, function(x) x[,-col_ix])
  X.train.group = matrix(0, nrow = length(Y.train), ncol = 115*L)
  vec = c(0, cumsum(sapply(X.list, nrow)))
  for (l in 1:L){
    X.train.group[(vec[l]+1):vec[l+1], 115*(l-1)+1:115] = X.list[[l]]
  }
  group_ix = rep(seq(1,115), L)
  
  if (library == "SGL"){
    # use SGL
    SGLdata = list(x = X.train.group, y = Y.train)
    SGL.list = SGL(data = SGLdata, index = group_ix, min.frac = 0.25, nlam = 10, alpha = alpha, standardize = F)
    #beta.list = lapply(1:10, function(lam_ix) add.diagonal(matrix(SGL.list$beta[,lam_ix], ncol = L), col_ix))
    group_coef = array(SGL.list$beta, dim = c(115, L, 10))
  }else if (library == "grplasso"){
    # use grplasso
    lambda_max = lambdamax(x = X.train.group, y = Y.train, index = group_ix, model = LinReg(), center = F, standardize = F)
    lambdas = exp(seq(log(lambda_max), log(lambda_max*0.25), length.out = 10))
    group_fit = grplasso(x = X.train.group, y = Y.train, lambda = lambdas, index = group_ix, model = LinReg(), center = F, standardize = F)
    group_coef = array_reshape(group_fit$coefficients, c(115, L, 10), order = "F")
  }
  return(group_coef)
}

# get res coef given group coef X.group.list
getX.res.parallel = function(data.concat.train.list, data.concat.test.list, X.group.list, lambda_ix, col_ix, library = "SGL", alpha = 0.8){
  L = length(data.concat.train.list)
  group_coef = X.group.list[[col_ix]] # group_coef to be 116 by 115 by L
  X.list = list()
  for (l_test in 1:L){
    res.list = lapply(1:L, function(l) data.concat.train.list[[l]][,col_ix] - data.concat.train.list[[l]][,-col_ix]%*%group_coef[,l_test])
    # given residuals refit every subject
    coef.mtx.list = list()
    for (l in 1:L){
      coef.list = list()
      for (subject_ix in 1:n.train.vec[l]){
        res_fit = glmnet(x = data.concat.train.list[[l]][137*(subject_ix-1)+1:137, -col_ix], y = res.list[[l]][137*(subject_ix-1)+1:137], standardize = F, nlambda = 10, lambda.min.ratio = 0.25, intercept = F)
        coef.list[[subject_ix]] = coef(res_fit, s=res_fit$lambda[lambda_ix])[-1]
      }
      coef.mtx = do.call(cbind, coef.list)
      coef.mtx.list[[l]] = coef.mtx
    }
    X.train = as.matrix(t(do.call(cbind, coef.mtx.list)))
    
    res.test.list = lapply(1:L, function(l) data.concat.test.list[[l]][,col_ix] - data.concat.test.list[[l]][,-col_ix]%*%group_coef[,l_test])
    # given residuals refit every subject
    coef.mtx.test.list = list()
    for (l in 1:L){
      coef.list = list()
      for (subject_ix in 1:n.test.vec[l]){
        res_fit = glmnet(x = data.concat.test.list[[l]][137*(subject_ix-1)+1:137, -col_ix], y = res.test.list[[l]][137*(subject_ix-1)+1:137], standardize = F, nlambda = 10, lambda.min.ratio = 0.25, intercept = F)
        coef.list[[subject_ix]] = coef(res_fit, s=res_fit$lambda[lambda_ix])[-1]
      }
      coef.mtx = do.call(cbind, coef.list)
      coef.mtx.test.list[[l]] = coef.mtx
    }
    X.test = as.matrix(t(do.call(cbind, coef.mtx.test.list)))
    X.list[[l_test]] = list(X.train = X.train, X.test = X.test)
  }
  names(X.list) = c('res1coef', 'res2coef')
  return(X.list) 
}

getX.group.parallel = function(data.concat.train.list, data.concat.test.list, lambda_group_ix, lambda_ix, col_ix, library = "SGL", alpha = 0.8){
  L = length(data.concat.train.list)
  Y.train = do.call(c, lapply(data.concat.train.list, function(x) x[,col_ix]))
  X.list = lapply(data.concat.train.list, function(x) x[,-col_ix])
  X.train.group = matrix(0, nrow = length(Y.train), ncol = 115*L)
  vec = c(0, cumsum(sapply(X.list, nrow)))
  for (l in 1:L){
    X.train.group[(vec[l]+1):vec[l+1], 115*(l-1)+1:115] = X.list[[l]]
  }
  group_ix = rep(seq(1,115), L)
  
  if (library == "SGL"){
    # use SGL
    SGLdata = list(x = X.train.group, y = Y.train)
    SGL.list = SGL(data = SGLdata, index = group_ix, min.frac = 0.25, nlam = 10, alpha = alpha, standardize = F)
    #beta.list = lapply(1:10, function(lam_ix) add.diagonal(matrix(SGL.list$beta[,lam_ix], ncol = L), col_ix))
    group_coef = array(SGL.list$beta, dim = c(115, L, 10))
  }else if (library == "grplasso"){
    # use grplasso
    lambda_max = lambdamax(x = X.train.group, y = Y.train, index = group_ix, model = LinReg(), center = F, standardize = F)
    lambdas = exp(seq(log(lambda_max), log(lambda_max*0.25), length.out = 10))
    group_fit = grplasso(x = X.train.group, y = Y.train, lambda = lambdas, index = group_ix, model = LinReg(), center = F, standardize = F)
    group_coef = array_reshape(group_fit$coefficients, c(115, L, 10), order = "F")
  }
  
  X.list = list()
  for (l_test in 1:L){
    res.list = lapply(1:L, function(l) data.concat.train.list[[l]][,col_ix] - data.concat.train.list[[l]][,-col_ix]%*%group_coef[,,lambda_group_ix][,l_test])
    # given residuals refit every subject
    coef.mtx.list = list()
    for (l in 1:L){
      coef.list = list()
      for (subject_ix in 1:n.train.vec[l]){
        res_fit = glmnet(x = data.concat.train.list[[l]][137*(subject_ix-1)+1:137, -col_ix], y = res.list[[l]][137*(subject_ix-1)+1:137], standardize = F, nlambda = 10, lambda.min.ratio = 0.25, intercept = F)
        coef.list[[subject_ix]] = coef(res_fit, s=res_fit$lambda[lambda_ix])[-1]
      }
      coef.mtx = do.call(cbind, coef.list)
      coef.mtx.list[[l]] = coef.mtx
    }
    X.train = as.matrix(t(do.call(cbind, coef.mtx.list)))
    
    res.test.list = lapply(1:L, function(l) data.concat.test.list[[l]][,col_ix] - data.concat.test.list[[l]][,-col_ix]%*%group_coef[,,lambda_group_ix][,l_test])
    # given residuals refit every subject
    coef.mtx.test.list = list()
    for (l in 1:L){
      coef.list = list()
      for (subject_ix in 1:n.test.vec[l]){
        res_fit = glmnet(x = data.concat.test.list[[l]][137*(subject_ix-1)+1:137, -col_ix], y = res.test.list[[l]][137*(subject_ix-1)+1:137], standardize = F, nlambda = 10, lambda.min.ratio = 0.25, intercept = F)
        coef.list[[subject_ix]] = coef(res_fit, s=res_fit$lambda[lambda_ix])[-1]
      }
      coef.mtx = do.call(cbind, coef.list)
      coef.mtx.test.list[[l]] = coef.mtx
    }
    X.test = as.matrix(t(do.call(cbind, coef.mtx.test.list)))
    X.list[[l_test]] = list(X.train = X.train, X.test = X.test)
  }
  names(X.list) = c('res1coef', 'res2coef')
  return(X.list) 
}

myT.test = function(mean1, sd1, mean2, sd2, n.train.vec){
  n1 = n.train.vec[1]
  n2 = n.train.vec[2]
  Tstat = (mean1 - mean2)/sqrt(sd1^2/n1 +sd2^2/n2)
  df = (sd1^2/n1 + sd2^2/n2)^2/((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))
  if (is.na(Tstat) || is.na(df)){
    return(NaN)
  }else{
    return(((pt(Tstat, df)<0.5)*pt(Tstat, df) + (pt(Tstat, df) >=0.5)*(1-pt(Tstat, df)))*2)
  }
}

mtx.feature.ix = function(G1.train, label.train){
  mtx1.mean = apply(G1.train[,,label.train==1], c(1,2), mean)
  mtx2.mean = apply(G1.train[,,label.train==2], c(1,2), mean)
  
  mtx1.sd = apply(G1.train[,,label.train==1], c(1,2), sd)
  mtx2.sd = apply(G1.train[,,label.train==2], c(1,2), sd)
  
  p.vec = c()
  for (c in 1:116){
    for (r in 1:116){
      p.vec = c(p.vec, myT.test(mtx1.mean[r,c], mtx1.sd[r,c], mtx2.mean[r,c], mtx2.sd[r,c], n.train.vec))
    }
  }
  p.mtx = matrix(p.vec, ncol = 116)
  TF.mtx = p.mtx < 0.05
  TF.mtx[is.na(TF.mtx)] = F
  return(list(TF.mtx = TF.mtx, p.mtx = p.mtx))
}










