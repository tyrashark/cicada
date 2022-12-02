
ECM_raw <- function(x, K, p0=NULL, mu0=NULL, nu0=NULL, initial = c("kmeans", "hclust", "select"), dyar = T, commonvar = T, methods = c("Projection", "EMP", "Lasso"), tol = 1e-5, max_it = 5e3, delete = NULL){
  
  # check deleted samples
  if(length(delete) > 0) x_del = x[delete]
  x = x[!((1:length(x)) %in% delete)]
  n = length(x)
  
  # chose initial values for the optimization
  if(initial == "kmeans"){
    ## generate initial value by K-means clustering
    result = kmeans(x, K, iter.max = 100, nstart = 1000)
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
  }else if(initial == "hclust"){
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
  if(dyar == T){
    M = cbind(rep(1, K), 0:(K-1))
    logmu0 = log(mu0)
    mu0 = exp(M %*% solve(t(M)%*%M)%*%t(M)%*% logmu0)
  }
  
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
 
  
  # The optimization
  n_iter = 0
  conv = 0
  z_hat = matrix(0, nrow = n, ncol = K)
  X = matrix(x, nrow=n, ncol=K, byrow = F)
 
  mu_t = mu0
  nu_t = nu0
  p_t = p0
  

  
  # EMP preparation
  if(methods == "EMP"){
    B = matrix(0, nrow = K-2, ncol =K)
    B[1,] = c(1,-2,1, rep(0, K-3))
    for(k in 1:(K-2)){
      B[k,] = c(rep(0, k-1), c(1,-2,1), rep(0, K-(k+2)))
    }
    
    ## derivative of Q by mu
    dmuQ = function(mu, nu, z){
      nu_v = vector("numeric", K) * nu
      apply(z*(X - matrix(mu, nrow = n, ncol = K, byrow = T))%*%diag(nu_v), 2, sum)
    }
  }
  
  while(n_iter < max_it){
    ## E-step
    z_hat = mapply(FUN = function(pk, muk, nuk) pk*dnorm(x, muk, sd=1/sqrt(nuk)) / f_obs_v(x, p_t, mu_t, nu_t), p_t, mu_t, nu_t)
    
    ## M-step
    p_up = apply(z_hat, 2, mean)
    
    if(dyar == T){
      if(commonvar == F){
        ### M-step of EM with non-constraints
        mu_up = apply(z_hat*X, 2, sum) / apply(z_hat, 2, sum)
        if(methods =="projection"){
          mu_up_p = exp(M %*% solve(t(M)%*%M)%*%t(M)%*% log(mu_up))
          mu_up = mu_up_p
        }
        nu_up = apply(z_hat, 2, sum) / apply(z_hat*(matrix(mu_up, nrow = n, ncol = K, byrow = T) - X)^2, 2, sum)
      }else{
        ### M-step of EM with common variance
        mu_up = apply(z_hat*X, 2, sum) / apply(z_hat, 2, sum)
        if(methods == "projection"){
          mu_up_p =  exp(M %*% solve(t(M)%*%M)%*%t(M)%*% log(mu_up))
          mu_up = mu_up_p
        }
        nu_up = sum(z_hat) / sum(z_hat*(matrix(mu_up, nrow = n, ncol = K, byrow = T) - X)^2)
      }
      if(methods == "EMP"){
        ## P-step
        if(K >= 3){
          mu_j = mu_up
          nu_EM = vector("numeric", K) * nu_up
          p_iter = 0
          s = 1
          B_mu = B * matrix(1/mu_j, nrow = K-2, ncol = K, byrow = T)
          ### P-step condition
          Q_incre = crossprod(dmuQ(z_hat, mu_j, nu_EM), t(B_mu)%*%solve(B_mu%*%t(B_mu))%*%B%*%log(mu_j))/crossprod(dmuQ(z_hat, mu_j, nu_EM), (diag(1, K) - P_B)%*%dmuQ(z_hat, mu_j, nu_EM))
          a0 = 0.001
          while(s <= 50){
            a = a0*(1/2)^s  
            if(a < 2*abs(Q_incre)) break
            s = s+1
          }
          a = a* sign(Q_incre)
          if(s > 1){
            while(p_iter < 50){
              P_B = t(B_mu)%*%solve(B_mu%*%t(B_mu))%*%B_mu
              mu_j_up = mu_j - t(B_mu)%*%solve(B_mu%*%t(B_mu))%*%B%*%log(mu_j) + a*(diag(1, K) - P_B)%*%dmuQ(z_hat, mu_j, nu_EM)
              #### P-step convergence
              if(sum(abs(mu_j - mu_j_up))<1e-2) break
              mu_j = mu_j_up
              B_mu = B * matrix(1/mu_j, nrow = K-2, ncol = K, byrow = T)
              p_iter = p_iter + 1
            }
            mu_up = mu_j_up
          }
        }
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
  
  z_max = apply(z_hat, 1, which.max)
  Z_map = matrix(F, nrow = n, ncol = K)
  
  Z_map = t(mapply(FUN = function(z, Z){
    Z[z] = T
    return(Z)
  }, z_max, split(Z_map, row(Z_map))))
  
  ICL = BIC - 2*sum(log(z_hat[Z_map]))
  
  names(BIC) = "BIC"
  names(ICL) = "ICL"
  #  names(loglikelihood) = "log-likelihood"
  #  names(npara) = "df"
  #  names(n) = "sample size"
  
  obj = NULL
  obj$x = x
  obj$n = n
  
  if(length(delete)>0) obj$x_del = x_del
  
  obj$initial_value = c(p0, mu0, nu0)
  obj$hat_z = z_hat
  obj$est = list(p = p_t, mu = mu_t, nu = nu_t)
  obj$criteria = c(BIC, ICL)
  obj$tol = tol
  obj$iterations = n_iter
  obj$conv = (conv == 1)
  obj$log_likelihood = log(f_obs_v(x, p_t, mu_t, nu_t))
  obj$density = function(x) f_obs_v(x, p_t, mu_t, nu_t)
  
  return(obj)
}

