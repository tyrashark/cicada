if("doParallel" %in% rownames(installed.packages()) == FALSE)
{install.packages("doParallel")
}
library("doParallel")

if("parallel" %in% rownames(installed.packages()) == FALSE)
{install.packages("parallel")
}
library("parallel")


Mixden <- function(x, p, mu, sigma){
  output = mapply(FUN = function(pk, muk, sigk) pk*dnorm(x, muk, sd = sigk), p, mu, sigma)
  return(apply(output, 1, sum))
}


initial_value <- function(ini, x, n, K, p0=NULL, mu0=NULL, nu0=NULL, dyar, commonvar, rep_kmean = 1000){
  # chose initial values for the optimization
  if(ini == "kmeans"){
    ## generate initial value by K-means clustering
    result = kmeans(x, K, iter.max = 100, nstart = rep_kmean)
    result = kmeans(x, K, iter.max = 100, centers = sort(result$centers))
    mu0 = result$centers
    n_K = as.vector(table(result$cluster))
    p0 = n_K/n
    if(commonvar == T){
      nu0 = (n-K)/sum(result$withinss) 
    }else{
      withinss = result$withinss
      #smoothing variance
      withinss[withinss < 0.001] = mean(withinss)
      nu0 = (n_K-1)/result$withinss
    }
  }else if(ini == "hclust"){
    ## generate initial value by hierarchical clustering
    sort_x = sort(x)
    result = hclust(dist(sort_x), method = "average")
    n_K = as.vector(table(cutree(result, k=K)))
    p0 = n_K/n
    
    mu0 = vector("numeric", K)
    s.var = vector("numeric", K)
    for(k in 1:K){
      mu0[k] = mean(sort_x[cutree(result, k=K) == k])
      if(n_K[k] >= 2){
        s.var[k] = var(sort_x[cutree(result, k=K) == k])
      }else{
        s.var[k] = 0
      }
    }
    # smoothing variance
    s.var[s.var < 0.001] = mean(s.var)
    
    if(commonvar == T){
      nu0 = (n-K)/sum(s.var * (n_K-1))
    }else{
      nu0 = 1/s.var
    }
  }else{
    if(commonvar == T){
      if(!(length(p0) == K & length(mu0) == K & length(nu0)== 1)){
        stop("Warning: mismatched initial value dimension with hypothesis")
      }
    }else{
      if(!(length(p0) == K & length(mu0) == K & length(nu0)== K)){
        stop("Warning: mismatched initial value dimension with hypothesis")
      }
    }
  }
  
  init_value = NULL
  init_value$p0 = p0
  init_value$mu0 = mu0
  init_value$nu0 = nu0
  return(init_value)
}


h_optim <- function(x, n, K, ini_method, p0=NULL, mu0=NULL, nu0=NULL, dyar, commonvar, tol=1e-6, max_it=5e3, conv.deal=F){
  
  # local functions
  ## the mixture density of x
  f_obs_v <- function(x, p, mu, nu){
    sigma = 1/sqrt(nu)
    output = mapply(FUN = function(pk, muk, sigk) pk*dnorm(x, muk, sd = sigk), p, mu, sigma)
    return(apply(output, 1, sum))
  }
  
  ## observed loglikelihood
  l_obs <- function(p, mu, nu, x){
    sum(log(f_obs_v(x, p, mu, nu)))
  }
  ## missing loglikelihood
  l_mis <- function(p, mu, nu, x, z){
    sigma = 1/sqrt(nu)
    output = mapply(FUN = function(zk, pk, muk, sigk) zk * (log(pk) + log(dnorm(x, muk, sd = sigk)) - log(f_obs_v(x, p, mu, nu))), 
                    split(z, col(z)), p, mu, sigma)
    return(sum(output))
  }
  ## complete loglikelihood
  l_comp <- function(p, mu, nu, x, z) {
    sigma = 1/sqrt(nu)
    output = mapply(FUN = function(zk, pk, muk, sigk) zk * (log(pk) + log(dnorm(x, muk, sd = sigk))), 
                    split(z, col(z)), p, mu, sigma)
    return(sum(output))
  }
  
  n_iter = 0
  conv = 0
  
  init_value = initial_value(ini_method, x, n, K, p0=p0, mu0=mu0, nu0=nu0, dyar=dyar, commonvar=commonvar)
  X = matrix(x, nrow=n, ncol=K, byrow = F)
  mu_t = init_value$mu0
  nu_t = init_value$nu0
  p_t = init_value$p0
  if(dyar == T){
    M = cbind(rep(1, K), 0:(K-1))
    mu_t = M %*% solve(t(M)%*%M)%*%t(M)%*% mu_t
  }
  max_research = 0 
  
  while(n_iter < max_it){
    ## E-step
    z_hat = mapply(FUN = function(pk, muk, nuk) pk*dnorm(x, muk, sd=1/sqrt(nuk)) / f_obs_v(x, p_t, mu_t, nu_t), p_t, mu_t, nu_t)
    if(is.nan(sum(z_hat))){
      if(max_research > 10){
        warning("NaN occurs during iteration")
        break
      }
    }
    
    ## M-step
    p_up = apply(z_hat, 2, mean)
    
    if(dyar == T){
      if(commonvar == F){
        ### M-step of ECM with constraints on mean vectors
        # 1) up mu given nu
        B_t = apply(z_hat, 2, sum) * nu_t
        d_t = apply(z_hat*X, 2, sum) * nu_t
        
        beta_up = tryCatch(
          solve(t(M)%*%diag(B_t)%*%M), error = function(e){
            matrix(NaN, nrow = ncol(M), ncol = ncol(M))
          }
        )%*%t(M)%*%d_t
        mu_up = M %*% beta_up
        
        # 2) extra E-step
        z_hat = mapply(FUN = function(pk, muk, nuk) pk*dnorm(x, muk, sd=1/sqrt(nuk)) / f_obs_v(x, p_up, mu_up, nu_t), p_up, mu_up, nu_t)
        
        # 3) update nu given updated mu and updated z_hat
        nu_up = apply(z_hat, 2, sum) / apply(z_hat*(matrix(mu_up, nrow = n, ncol = K, byrow = T) - X)^2, 2, sum)
      }else{
        ### M-step of EM with constraints on mean vectors & common variance
        beta_up = tryCatch(
          solve(t(M)%*%diag(apply(z_hat, 2, sum))%*%M), error = function(e){
            matrix(NaN, nrow = ncol(M), ncol = ncol(M))
          }
        )%*%t(M)%*%apply(z_hat*X, 2, sum)
        mu_up = M %*% beta_up
        nu_up = sum(z_hat) / sum(z_hat*(matrix(mu_up, nrow = n, ncol = K, byrow = T) - X)^2)
      }
    }
    if(dyar == F){
      if(commonvar == F){
        ### M-step of EM with non-constraints
        mu_up = apply(z_hat*X, 2, sum) / apply(z_hat, 2, sum)
        nu_up = apply(z_hat, 2, sum) / apply(z_hat*(matrix(mu_up, nrow = n, ncol = K, byrow = T) - X)^2, 2, sum)
      }else{
        ### M-step of EM with common variance
        mu_up = apply(z_hat*X, 2, sum) / apply(z_hat, 2, sum)
        nu_up = sum(z_hat) / sum(z_hat*(matrix(mu_up, nrow = n, ncol = K, byrow = T) - X)^2)
      }
    }
    
    ## Check tolerance for converging the algorithm
    error = abs(l_obs(p_up, mu_up, nu_up, x) - l_obs(p_t, mu_t, nu_t, x))
    
    p_t = p_up; mu_t = mu_up; nu_t = nu_up;
    n_iter = n_iter + 1
    if(is.nan(error)){
      if(max_research > 10){
        warning("NaN occurs during iteration")
        break
      }
      error = 1
      n_iter = 0
      init_value = initial_value(ini_method, x, n, K, p0, mu0, nu0, dyar, commonvar, rep_kmean = 1)
      mu_t = init_value$mu0
      nu_t = init_value$nu0
      p_t = init_value$p0
      if(dyar == T){
        M = cbind(rep(1, K), 0:(K-1))
        mu_t = M %*% solve(t(M)%*%M)%*%t(M)%*% mu_t
      }
      max_research = max_research + 1
    }
    
    if(error < tol){
      conv = 1
      if(dyar==T & commonvar==F){
        cat(paste("ECM algorithm is converged in",n_iter,"steps \n"))
      }
      else cat(paste("EM algorithm is converged in",n_iter,"steps \n"))
      break
    }
  }
  
  if(n_iter == max_it){
    if(dyar==T & commonvar==F){
      cat(paste("ECM algorithm reachs max number of iterations ", max_it, "steps \n"))
    }
    else cat(paste("EM algorithm reachs max number of iterations ", max_it, "steps \n"))
  }
  
  ## Model selection criteria
  if(dyar == T){
    npara = K + 2 - 1
  }else{
    npara = K * 2 - 1
  }
  if(commonvar == T){
    npara = npara + 1
  }else{
    npara = npara + K
  }
  
  loglikelihood = sum(log(f_obs_v(x, p_t, mu_t, nu_t)))  
  BIC = -2*loglikelihood + npara*log(n)
  LR = loglikelihood
  if(n_iter == max_it & conv.deal){
    LR = NaN
  }
  
  ICL = NA
  if(!is.nan(sum(z_hat))){
    z_max = apply(z_hat, 1, which.max)
    Z_map = matrix(F, nrow = n, ncol = K)
    
    Z_map = t(mapply(FUN = function(z, Z){
      Z[z] = T
      return(Z)
    }, z_max, split(Z_map, row(Z_map))))
    
    ICL = BIC - 2*sum(log(z_hat[Z_map]))
  }
  names(BIC) = "BIC"
  names(ICL) = "ICL"
  names(LR) = "LR"
  
  obj = NULL
  obj$x = x
  obj$n = n

  #obj$initial_value = c(p0, mu0, nu0)
  obj$hat_z = z_hat
  obj$est = list(p = c(p_t), mu = c(mu_t), nu = c(nu_t))
  obj$criteria = c(BIC, ICL, LR)
  #obj$tol = tol
  obj$iterations = n_iter
  obj$conv = (conv == 1)
  obj$log_likelihood = log(f_obs_v(x, p_t, mu_t, nu_t))
  #obj$density = function(x) f_obs_v(x, p_t, mu_t, nu_t)
  obj$cluster = table(factor(apply(z_hat, 1, which.max), levels = 1:K))
  obj$classified = apply(z_hat, 1, which.max)
  obj$research = max_research
  return(obj)
}


ECM <- function(x, K, p0=NULL, mu0=NULL, nu0=NULL, ini_method = c("kmeans", "hclust", "select"), tol = 1e-6, max_it = 5e3, delete = NULL){
  
  # check deleted samples
  if(length(delete) > 0) x_del = x[delete]
  x = x[!((1:length(x)) %in% delete)]
  n = length(x)
 
  # The optimization
  h1 = h_optim(x, n, K, ini_method, p0, mu0, nu0, dyar = F, commonvar = F, tol, max_it, conv.deal = F)
  h2 = h_optim(x, n, K, ini_method, p0, mu0, nu0, dyar = F, commonvar = T, tol, max_it, conv.deal = F)
  h3 = h_optim(x, n, K, ini_method, p0, mu0, nu0, dyar = T, commonvar = F, tol, max_it, conv.deal = F)
  h4 = h_optim(x, n, K, ini_method, p0, mu0, nu0, dyar = T, commonvar = T, tol, max_it, conv.deal = F)
  
  result = NULL
  p_est = cbind(h1$est[["p"]], h2$est[["p"]], h3$est[["p"]], h4$est[["p"]])
  mu_est = cbind(h1$est[["mu"]], h2$est[["mu"]], h3$est[["mu"]], h4$est[["mu"]])
  sigma_est = cbind(1/h1$est[["nu"]], 1/h2$est[["nu"]]*rep(1,K), 1/h3$est[["nu"]], 1/h4$est[["nu"]]*rep(1,K))
  colnames(p_est) = paste0("h", 1:4)
  colnames(mu_est) = paste0("h", 1:4)
  colnames(sigma_est) = paste0("h", 1:4)
  result$est = list(p = p_est, mu = mu_est, sigma_sq = sigma_est)
  
  cluster = rbind(h1$cluster, h2$cluster, h3$cluster, h4$cluster)
  rownames(cluster) = paste0("h", 1:4)
  result$cluster = cluster
  
  criteria = rbind(h1$criteria, h2$criteria, h3$criteria, h4$criteria)
  rownames(criteria) = paste0("h", 1:4)
  result$criteria = criteria
  
  result$log_likelihood = list(h1 = h1$log_likelihood, h2 = h2$log_likelihood, h3 = h3$log_likelihood, h4 = h4$log_likelihood)
  result$hat_z = list(h1 = h1$hat_z, h2 = h2$hat_z, h3 = h3$hat_z, h4 = h4$hat_z)
  result$classified = list(h1 = h1$classified, h2 = h2$classified, h3 = h3$classified, h4 = h4$classified)
  result$num_of_research = list(h1 = h1$research, h2 = h2$research, h3 = h3$research, h4 = h4$research)
  if(length(delete)>0) result$x_del = x_del
  
  return(result)
}


compare_h <- function(x, K, hypothesis = paste0("h", 1:4), p0=NULL, mu0=NULL, nu0=NULL, ini_method = c("kmeans", "hclust", "select"), tol = 1e-6, max_it = 5e3, conv.deal = T){
  
  hypothesis = sort(hypothesis)
  n = length(x)
  
  h1=NULL;h2=NULL;h3=NULL;h4=NULL
  # The optimization
  if("h1" %in% hypothesis) h1 = h_optim(x, n, K, ini_method, p0, mu0, nu0, dyar = F, commonvar = F, tol, max_it, conv.deal=conv.deal)
  if("h2" %in% hypothesis) h2 = h_optim(x, n, K, ini_method, p0, mu0, nu0, dyar = F, commonvar = T, tol, max_it, conv.deal=conv.deal)
  if("h3" %in% hypothesis) h3 = h_optim(x, n, K, ini_method, p0, mu0, nu0, dyar = T, commonvar = F, tol, max_it, conv.deal=conv.deal)
  if("h4" %in% hypothesis) h4 = h_optim(x, n, K, ini_method, p0, mu0, nu0, dyar = T, commonvar = T, tol, max_it, conv.deal=conv.deal)
  
  result = NULL
  p_est = cbind(h1$est[["p"]], h2$est[["p"]], h3$est[["p"]], h4$est[["p"]])
  mu_est = cbind(h1$est[["mu"]], h2$est[["mu"]], h3$est[["mu"]], h4$est[["mu"]])
  sigma_est = cbind(1/h1$est[["nu"]], 1/h2$est[["nu"]]*rep(1,K), 1/h3$est[["nu"]], 1/h4$est[["nu"]]*rep(1,K))
  
  result$est = list(p = p_est, mu = mu_est, sigma_sq = sigma_est)
  
  criteria = rbind(h1$criteria, h2$criteria, h3$criteria, h4$criteria)

  result$criteria = criteria
  
  return(result)
}



### bootstrap

para_boot_LR <- function(x, K, bootnum, delete= NULL, max_it = 5e3, max_it.boot = 1e2, tol = 1e-6, tol.boot = 1e-6, conv.deal = T){
  
  # check deleted samples
  if(length(delete) > 0) x_del = x[delete]
  x = x[!((1:length(x)) %in% delete)]
  bootsize = length(x)
  

  origin = compare_h(x, K=K, hypothesis = paste0("h",1:4), ini_method = "kmeans", max_it=max_it, tol = tol)
  BIC = origin$criteria[,1]
  ICL = origin$criteria[,2]
  LR = origin$criteria[,3]
  
  bootstrap_h3null = lapply(1:bootnum, FUN = function(i){
    z_bnull = rmultinom(bootsize, 1, prob = origin$est$p[,3])
    x_bnull = mapply(FUN = function(mu, sigma) rnorm(bootsize, mean = mu, sd = sqrt(sigma)), origin$est$mu[,3], origin$est$sigma_sq[,3])
    return(apply(x_bnull * t(z_bnull), 1, sum))
  })
  
  bootstrap_h4null = lapply(1:bootnum, FUN = function(i){
    z_bnull = rmultinom(bootsize, 1, prob = origin$est$p[,4])
    x_bnull = mapply(FUN = function(mu, sigma) rnorm(bootsize, mean = mu, sd = sqrt(sigma)), origin$est$mu[,4], origin$est$sigma_sq[,4])
    return(apply(x_bnull * t(z_bnull), 1, sum))
  })
  bootstrap_h2null = lapply(1:bootnum, FUN = function(i){
    z_bnull = rmultinom(bootsize, 1, prob = origin$est$p[,2])
    x_bnull = mapply(FUN = function(mu, sigma) rnorm(bootsize, mean = mu, sd = sqrt(sigma)), origin$est$mu[,2], origin$est$sigma_sq[,2])
    return(apply(x_bnull * t(z_bnull), 1, sum))
  })
  
  result_boot = list()
  for(i in 1:bootnum){
    temp = NULL
    #h3,h1
    temp$h3null = compare_h(bootstrap_h3null[[i]], K=K, hypothesis = c("h3","h1"), ini_method = "kmeans", max_it=max_it.boot, tol=tol.boot,  conv.deal=conv.deal)
    ##h4,h3,h2,h1
    temp$h4null = compare_h(bootstrap_h4null[[i]], K=K, hypothesis = paste0("h",1:4), ini_method = "kmeans", max_it=max_it.boot, tol=tol.boot,  conv.deal= conv.deal)
    #h2,h1
    temp$h2null = compare_h(bootstrap_h2null[[i]], K=K, hypothesis = c("h2","h1"), ini_method = "kmeans", max_it=max_it.boot, tol=tol.boot,  conv.deal= conv.deal)
    
    result_boot[[i]] = temp
  }
  
  LR_h4null = matrix(0, nrow = bootnum, ncol = 4)
  LR_h3null = matrix(0, nrow = bootnum, ncol = 2)
  LR_h2null = matrix(0, nrow = bootnum, ncol = 2)
  
  for(i in 1:bootnum){
    LR_h3null[i, ] = result_boot[[i]]$h3null$criteria[,3]
    LR_h4null[i, ] = result_boot[[i]]$h4null$criteria[,3]
    LR_h2null[i, ] = result_boot[[i]]$h2null$criteria[,3]
  }
  
  LRbooth3h1 = -2*(LR_h3null[1:bootnum,2] - LR_h3null[1:bootnum,1])
  LRbooth4h3 = -2*(LR_h4null[1:bootnum,4] - LR_h4null[1:bootnum,3])
  LRbooth4h2 = -2*(LR_h4null[1:bootnum,4] - LR_h4null[1:bootnum,2])
  LRbooth4h1 = -2*(LR_h4null[1:bootnum,4] - LR_h4null[1:bootnum,1])
  LRbooth2h1 = -2*(LR_h2null[1:bootnum,2] - LR_h2null[1:bootnum,1])
  LRh3h1 = -2*(LR[3] - LR[1])
  LRh4h3 = -2*(LR[4] - LR[3])
  LRh4h2 = -2*(LR[4] - LR[2])
  LRh4h1 = -2*(LR[4] - LR[1])
  LRh2h1 = -2*(LR[2] - LR[1])
  
  nanh3h1 = sum(is.nan(LRbooth3h1))
  nanh4h3 = sum(is.nan(LRbooth4h3))
  nanh4h2 = sum(is.nan(LRbooth4h2))
  nanh4h1 = sum(is.nan(LRbooth4h1))
  nanh2h1 = sum(is.nan(LRbooth2h1))
  
  return(data.frame(h3h1 = sum(LRbooth3h1 > LRh3h1, na.rm = T)/(bootnum+1-nanh3h1), h4h3 = sum(LRbooth4h3 > LRh4h3, na.rm = T)/(bootnum+1-nanh4h3),
                      h4h2 = sum(LRbooth4h2 > LRh4h2, na.rm = T)/(bootnum+1-nanh4h2), h4h1 = sum(LRbooth4h1 > LRh4h1, na.rm = T)/(bootnum+1-nanh4h1), h2h1 = sum(LRbooth2h1 > LRh2h1, na.rm = T)/(bootnum+1-nanh2h1),
                      nanh3h1 = nanh3h1, nanh4h3 = nanh4h3, nanh4h2 = nanh4h2, nanh4h1 = nanh4h1, nanh2h1 = nanh2h1))
}


find_comp <- function(x, K = 3:8, ini_method = "kmeans", max_it = 1e3, tol = 1e-6){
  K_min = min(K)
  K_max = max(K)
  BIC = matrix(0, nrow = K_max-K_min+1, ncol = 4)
  rownames(BIC) = K_min:K_max
  colnames(BIC) = paste0("h", 1:4)
  
  ICL = BIC
  result = list()
  
  for(k in 1:(K_max-K_min+1)){
    result[[k]] = ECM(x, K=k+K_min-1, ini_method = ini_method, delete = NULL, max_it = max_it, tol = tol)
    BIC[k, ] = result[[k]]$criteria[,1]
    ICL[k, ] = result[[k]]$criteria[,2]
  }
  
  op_comp = apply(ICL, 2, which.min)+K_min - 1
  result = data.frame("h1" = op_comp[1], "h2" = op_comp[2], "h3" = op_comp[3], "h4" = op_comp[4], "min_hypo" =  paste0("h", which.min(apply(ICL, 2, min))))
  rownames(result) = ""
  return(result)
}



hypo_test <- function(tab, alpha = 0.05, adj.method = c("bonferroni", "punzo")){
  
  ##Bonferroni correction
  
  if(adj.method == "bonferroni"){
    test = mapply(FUN = function(h4h1, h3h1, h2h1){
      htest = c(F,F,F)
      names(htest) = c("h4", "h3", "h2")
      
      if(h4h1 < alpha){
        htest[1] = T
       if(h3h1 < alpha/2) htest[2] = T
       if(h2h1 < alpha/2) htest[3] = T
      }
      
      if(htest[1] == F){
        result = "h4"
      }else if(htest[2] == F & htest[3] == T){
        result = "h3"
      }else if(htest[2] == T & htest[3] == F){
        result = "h2"
      }else if(htest[2] == T & htest[3] == T){
        result = "h1"
      }else{
        result = "h4"
      }
      return(result)
    }, tab$h4h1, tab$h3h1, tab$h2h1)
  }
  if(adj.method == "punzo"){
    test = mapply(FUN = function(h4h1, h3h1, h2h1){
      htest = c(F,F,F)
      names(htest) = c("h4", "h3", "h2")
      
      if(h4h1 < alpha){
        htest[1] = T
        if(max(h4h1, h3h1) < alpha) htest[2] = T
        if(max(h4h1, h2h1) < alpha) htest[3] = T
      }
      
      if(htest[1] == F){
        result = "h4"
      }else if(htest[2] == F & htest[3] == T){
        result = "h3"
      }else if(htest[2] == T & htest[3] == F){
        result = "h2"
      }else if(htest[2] == T & htest[3] == T){
        result = "h1"
      }else{
        result = "h4"
      }
      return(result)
    }, tab$h4h1, tab$h3h1, tab$h2h1)
  }
  
  result = NULL
  result$hypo = table(factor(test, levels = paste0("h", 1:4)))
  return(result)
}


minimize_ICL <- function(tab, Kmin=3, Kmax=8){
  min_hypo = table(factor(tab[, "min_hypo"], levels = paste0("h", 1:4)))
  ind = as.numeric(substr(tab[, "min_hypo"],2,2))
  selection = c()
  for(i in 1:length(ind)){
    selection = c(selection, tab[i, ind[i]])
  }
  result = NULL
  result$selection = table(factor(selection, levels = Kmin:Kmax))
  result$min_hypo = min_hypo
  return(result)
}





