setwd("~/")
source("find_mus.R")
source("Theta.hat_CLIME.R")
libraries = c("mvtnorm", "Matrix", "hdm", "doSNOW", "quantreg", "sandwich", "dplyr" ,"Rmosek", "doMC", "matrixStats")
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

c1 = makeCluster(16) # core numbers
registerDoSNOW(c1)

sim_scene = function(n,p,m,q,tau,df,H,W,rho,beta,idx,gamma){
  Pi = kronecker(matrix(rep(1,q/m),nrow=q/m,ncol=1),diag(m))*(q/m+q/m*0.5^m)^(-1/(q/m))
  u = array(0, dim=c(n,p,m+2))
  innov_Z = eta_Z = M_Z = list(length=p)
  for(i in 1:p){
    u[,i,] = rmvnorm(n, rep(0,m+1), diag(m+1))
    innov_Z[[i]] = numeric(0)
    for(j in 1:q){
      innov_Z[[i]] = cbind(innov_Z[[i]], rt(n+500+1, df)/sqrt(df/(df-2)))
    }
    eta_Z[[i]] = matrix(0,n+500,q)
    for(t in (1+1):(n+1+500)){
      eta_Z[[i]][t-1,] = innov_Z[[i]][t,]*sqrt(0.8*innov_Z[[i]][t-1,]^2+0.2)
    }
    M_Z[[i]] = list()
    for(k in 1:(500+1)){
      M_Z[[i]][[k]] = matrix(rnorm(q*q), nrow=q, ncol=q)
    }
  }
  A_Z = list(length=p)
  for(i in 1:p){
    A_Z[[i]] = list()
    for(j in 1:(500+1)){
      A_Z[[i]][[j]] = (j-1+1)^(-tau-1)*M_Z[[i]][[j]]
    }
  }
  Z = array(0, dim=c(n,p,q))
  for(i in 1:p){
    for(t in 1:(n)){
      for(tt in 1:(500+1)){
        Z[t,i,] = Z[t,i,] + A_Z[[i]][[tt]]%*%eta_Z[[i]][t-tt+500+1,]
      }
    }
    Z[,i,] = scale(Z[,i,])
  }
  X = array(0, dim=c(n,p,m))
  eps = matrix(0,n,p)
  for(i in 1:p){
    eps[,i] = u[,i,1]  
    v = u[,i,2:(m+1)]
    X[,i,] = Z[,i,]%*%Pi + v
  }
  Y = matrix(0, nrow=n, ncol=p)
  invmat = solve(diag(p) - rho*H)
  for(t in 1:n){
    Y[t,] = invmat%*%(eps[t,] + X[t,,]%*%gamma)
  }
  
  ## estimation step 1: RMD for beta
  # preliminary estimation 
  sigmas = matrix(0,p,q)
  ghats = epshats = numeric(0)
  for(i in 1:p){
    XX = cbind(Y%*%W[i,],X[,i,],Y[,-sort(c(i,which(W[i,]!=0)))])
    pp = ncol(XX)
    fit0 = rlasso(XX, Y[,i], intercept = FALSE, penalty = list(homoscedastic = "none", lambda.start = 1.5*sqrt(n)*qnorm(1-0.1/(2*pp))), post=FALSE)
    epshat = fit0$residuals
    for(j in 1:q){
      sigmas[i,j] = sd(Z[,i,j]*epshat)
      ghats = c(ghats, Z[,i,j]*epshat)
    }
    epshats = c(epshats, epshat)
  }
  # arrange the constraints (A*beta = b)
  A = matrix(0, nrow = (p*q), ncol = (1+m+p^2))
  for(i in 1:p){
    for(j in 1:q){
      A[((i-1)*p+j),1] = Z[,i,j]%*%Y%*%W[i,]/n
      A[((i-1)*p+j),2:(m+1)] = Z[,i,j]%*%X[,i,]/n
      A[((i-1)*p+j),(m+2+(p)*(i-1)):(m+1+(p)*i)] = Z[,i,j]%*%Y/n
    }
  }
  b = rep(0, (p*q))
  for(i in 1:p){
    for(j in 1:q)
      b[((i-1)*p+j)] = Z[,i,j]%*%Y[,i]/n
  }
  pp = ncol(A)
  con = 1.1
  lambda = qnorm(1-0.1/(2*q*p))/sqrt(n)*sd(ghats)*con
  con1 = cbind(A, -A)
  con2 = c(1,rep(0,pp-1),-1,rep(0,pp-1))
  C = matrix(0, nrow = p+length(which(W!=0)), ncol = pp) # constraints of diag(Delta) = 0 and Delta[i,j]=0 if W[i,j]!=0
  for(i in 1:p){
    C[i,(m+1+p*(i-1)+i)] = 1
  }
  for(i in (p+1):(p+length(which(W!=0)))){
    C[i,(m+1+which(t(W)!=0)[i-p])] = 1
  }
  con3 = cbind(C, -C)
  A.sparse = rbind(con1, con2, con3) %>% Matrix(sparse=T)
  objective = rep(1, 2*pp)
  blc = c(-lambda + b, -1, rep(0,p+length(which(W!=0))))
  buc = c(lambda + b, 1, rep(0,p+length(which(W!=0))))
  blx = rep(0, 2*pp)
  bux = rep(Inf, 2*pp)
  mosek.prob = list(
    sense="min",
    c=objective,
    A=A.sparse,
    bc=rbind(blc, buc),
    bx=rbind(blx, bux)
  )
  mosek.res = try(mosek(mosek.prob, list(verbose=0)), silent=TRUE)
  beta.hat = mosek.res$sol$bas$xx[1:pp] - mosek.res$sol$bas$xx[(pp+1):(2*pp)]
  rho.hat = beta.hat[1]
  gamma.hat = beta.hat[2:(m+1)]
  Delta.hat = matrix(beta.hat[-c(1:(m+1))],p,byrow=TRUE)
  res = matrix(0,n,p)
  for(i in 1:p){
    res[,i] = Y[,i]-rho.hat*Y%*%W[i,]-X[,i,]%*%gamma.hat-Y%*%Delta.hat[i,]
  }
  g.hat = numeric(0)
  for(i in 1:p){
    g.hat = rbind(g.hat, t(Z[,i,])%*%res[,i]/n)
  }
  
  ## estimation step 2: Debiasing and testing
  Omega.inv = Omega = numeric(0)
  G1 = G2 = numeric(0)
  for(i in 1:p){
    G1 = rbind(G1, -t(Z[,i,])%*%cbind(Y%*%W[i,], X[,i,])/n)
    G2 = bdiag(G2, -t(Z[,i,])%*%Y[,-sort(c(i,which(W[i,]!=0)))]/n)
    omega = t(Z[,i,]*as.vector(res[,i]))%*%(Z[,i,]*as.vector(res[,i]))/n
    Omega = bdiag(Omega, omega)
    omega.inv = try(solve(omega), silent = TRUE)
    if(class(omega.inv)[1]=="try-error"){
      mus = find_mus(omega)*1.2
      omega.inv = inv_CLIME(omega, mus)
      omega.inv = omega.inv*(abs(omega.inv)<= abs(t(omega.inv)))+ t(omega.inv)*(abs(omega.inv)> abs(t(omega.inv)))
    }
    Omega.inv = bdiag(Omega.inv, omega.inv)
  }
  G2 = G2[,-1]
  G = cbind(G1, G2)
  thres = 0.1*sqrt(log(nrow(G))/n)
  G = G*as.numeric(abs(G)>=thres)
  Omega = Omega[,-1]
  Omega.inv = Omega.inv[,-1]
  
  G11 = G[,1]
  G22 = G[,-1]
  D1 = t(G11)%*%Omega.inv%*%G11
  D2 = t(G22)%*%Omega.inv%*%G22
  D1.inv = try(solve(D1), silent = TRUE)
  if(class(D1.inv)[1]=="try-error"){
    mus1 = find_mus(D1)*1.2
    D1.inv = inv_CLIME(D1, mus1)
    D1.inv = D1.inv*(abs(D1.inv)<= abs(t(D1.inv)))+ t(D1.inv)*(abs(D1.inv)> abs(t(D1.inv)))
  }
  D2.inv = try(solve(D2), silent = TRUE)
  if(class(D2.inv)[1]=="try-error"){
    mus2 = find_mus(D2)*1.2
    D2.inv = inv_CLIME(D2, mus2)
    D2.inv = D2.inv*(abs(D2.inv)<= abs(t(D2.inv)))+ t(D2.inv)*(abs(D2.inv)> abs(t(D2.inv)))
  }
  F = -t(G11)%*%Omega.inv%*%G22%*%D2.inv%*%t(G22)%*%Omega.inv%*%G11
  B = D1.inv - D1.inv%*%solve(1 + F%*%D1.inv)%*%F%*%D1.inv
  A = t(G11)%*%Omega.inv%*%(diag(nrow(G)) - G22%*%D2.inv%*%t(G22)%*%Omega.inv)
  beta.hat2 = as.numeric(beta.hat[1] - B%*%A%*%g.hat)
  stds = sqrt(diag(B%*%A%*%Omega%*%t(A)%*%t(B)))
  stat = as.numeric(sqrt(n)*beta.hat2/stds)
  
  G11 = G[,2:(m+1)]
  G22 = G[,-(2:(m+1))]
  D1 = t(G11)%*%Omega.inv%*%G11
  D2 = t(G22)%*%Omega.inv%*%G22
  D1.inv = try(solve(D1), silent = TRUE)
  if(class(D1.inv)[1]=="try-error"){
    mus1 = find_mus(D1)*1.2
    D1.inv = inv_CLIME(D1, mus1)
    D1.inv = D1.inv*(abs(D1.inv)<= abs(t(D1.inv)))+ t(D1.inv)*(abs(D1.inv)> abs(t(D1.inv)))
  }
  D2.inv = try(solve(D2), silent = TRUE)
  if(class(D2.inv)[1]=="try-error"){
    mus2 = find_mus(D2)*1.2
    D2.inv = inv_CLIME(D2, mus2)
    D2.inv = D2.inv*(abs(D2.inv)<= abs(t(D2.inv)))+ t(D2.inv)*(abs(D2.inv)> abs(t(D2.inv)))
  }
  F = -t(G11)%*%Omega.inv%*%G22%*%D2.inv%*%t(G22)%*%Omega.inv%*%G11
  B = D1.inv - D1.inv%*%solve(diag(m) + F%*%D1.inv)%*%F%*%D1.inv
  A = t(G11)%*%Omega.inv%*%(diag(nrow(G)) - G22%*%D2.inv%*%t(G22)%*%Omega.inv)
  beta.hats = as.numeric(beta.hat[2:(m+1)] - B%*%A%*%g.hat)
  beta.hat2 = c(beta.hat2, beta.hats)
  stds = sqrt(diag(B%*%A%*%Omega%*%t(A)%*%t(B)))
  stat = c(stat, as.numeric(sqrt(n)*beta.hats/stds))
  
  G11 = G[,((m+2):ncol(G))[1:50]]
  G22 = G[,-((m+2):ncol(G))[1:50]]
  D1 = t(G11)%*%Omega.inv%*%G11
  D2 = t(G22)%*%Omega.inv%*%G22
  D1.inv = try(solve(D1), silent = TRUE)
  if(class(D1.inv)[1]=="try-error"){
    mus1 = find_mus(D1)*1.2
    D1.inv = inv_CLIME(D1, mus1)
    D1.inv = D1.inv*(abs(D1.inv)<= abs(t(D1.inv)))+ t(D1.inv)*(abs(D1.inv)> abs(t(D1.inv)))
  }
  D2.inv = try(solve(D2), silent = TRUE)
  if(class(D2.inv)[1]=="try-error"){
    mus2 = find_mus(D2)*1.2
    D2.inv = inv_CLIME(D2, mus2)
    D2.inv = D2.inv*(abs(D2.inv)<= abs(t(D2.inv)))+ t(D2.inv)*(abs(D2.inv)> abs(t(D2.inv)))
  }
  F = -t(G11)%*%Omega.inv%*%G22%*%D2.inv%*%t(G22)%*%Omega.inv%*%G11
  B = D1.inv - D1.inv%*%solve(diag(50) + F%*%D1.inv)%*%F%*%D1.inv
  A = t(G11)%*%Omega.inv%*%(diag(nrow(G)) - G22%*%D2.inv%*%t(G22)%*%Omega.inv)
  delta.hat2 = beta.hat[idx[1:50]] - B%*%A%*%g.hat
  stds = sqrt(diag(B%*%A%*%Omega%*%t(A)%*%t(B)))
  stat2 = sqrt(n)*delta.hat2/stds
  
  rho.hat2 = beta.hat2[1]
  gamma.hat2 = beta.hat2[2:(m+1)]
  error.rho = rho - rho.hat
  error.rho2 = rho - rho.hat2
  error.gamma = sqrt(sum((gamma[1:m] - gamma.hat[1:m])^2))
  error.gamma2 = sqrt(sum((gamma[1:m] - gamma.hat2)^2))
  dz.1 = as.numeric(beta.hat!=0)[1:(m+1)]
  db2.1 = as.numeric(abs(stat)>=1.96)
  
  delta.hat = beta.hat[idx[1:50]]
  delta.hat2 = delta.hat2
  error.delta = sqrt(sum((beta[idx[1:50]] - delta.hat)^2))
  error.delta2 = sqrt(sum((beta[idx[1:50]] - delta.hat2)^2))
  dz.2 = as.numeric(beta.hat[idx[1:50]]!=0)
  db2.2 = as.numeric(abs(stat2)>=1.96)
  
  error = c(error.rho, error.rho2, error.gamma, error.gamma2, error.delta, error.delta2)
  
  list(error=error,db2.1=db2.1,dz.1=dz.1,db2.2=db2.2,dz.2=dz.2,beta.hat2=beta.hat2,delta.hat2=delta.hat2,beta.hat=beta.hat)#,Z=Z,u=u)
}

PP = 30 #or 50
probs = 0.8
rhos = c(0.7,0.5)
taus = c(1,0.1)
rep= 500 #100
df = 8 

for(l in 1:length(PP)){
  for(ll in 1:length(rhos)){
    for(lll in 1:length(probs)){
      for(llll in 1:length(taus)){
        p = PP[l]
        rho = rhos[ll]
        prob = probs[lll]
        tau = taus[llll]
        if(p<50){m = 30}else{m = 50}
        q = ceiling((p+m+1)/m)*m  #number of iv for each equation 
        if(q<100){n = 100}else{n = 200} 
        
        H0 = matrix(rbinom(p*p,1,prob=0.5),ncol=p,nrow=p) #actual
        diag(H0) = 0
        W0 = matrix(0,p,p) 
        W0[which(H0!=0)] = rbinom(length(which(H0!=0)),1,prob)   #observed
        H = matrix(0,p,p)
        for(i in 1:p){
          H[,i] = H0[,i]/pmax(rowSums(H0),rowSums(W0))
        }
        W = matrix(0,p,p)
        for(i in 1:p){
          W[,i] = W0[,i]/pmax(rowSums(H0),rowSums(W0))
        }
        h = as.vector(t(H))
        w = as.vector(t(W))
        delta = rho*(h-w)
        Delta = t(matrix(delta,p,p))
        gamma = c(rep(10,5),rep(5,3),rep(1,2),rep(0,m-10)) 
        beta = c(rho, gamma, delta)
        
        idx = numeric(0)
        for(i in 1:p){
          idx = c(idx, (i-1)*p+(1:p)[-sort(c(i,which(W[i,]!=0)))])
        }
        idx = idx+m+1
        
        dgp = list(H=H,W=W,beta=beta,idx=idx)
        save(dgp,file=paste("dgp_tau",tau,"_rho",rho,"_p",p,"_prob",prob,".dat",sep=""))
        
        results = foreach(l=1:rep, .packages=c("quantreg", "mvtnorm", "Matrix", "hdm", "sandwich", "matrixStats", "Rmosek", "dplyr"), .inorder=FALSE) %dopar%{ 
          source("find_mus.R")
          source("Theta.hat_CLIME.R")
          sim_scene(n=n,p=p,m=m,q=q,tau=tau,df=df,H=H,W=W,rho=rho,beta=beta,idx=idx,gamma=gamma)
        }
        
        save(results,file=paste("results_tau",tau,"_rho",rho,"_p",p,"_prob",prob,".dat",sep=""))
      }
    }
  }
}

rep = 500 #100
PP = 30 #or 50
probs = 0.8
rhos = c(0.7,0.5)
taus = c(1,0.1)
p = PP[1]
rho = rhos[1]
prob = probs[1]
tau = taus[1]
if(p<50){m = 30}else{m = 50}
load(file=paste("dgp_tau",tau,"_rho",rho,"_p",p,"_prob",prob,".dat",sep=""))
beta = dgp$beta
h = as.vector(t(dgp$H))
w = as.vector(t(dgp$W))
idx = dgp$idx
length(which(beta[idx]!=0))/length(idx)
length(which(h==0))/length(h)
length(which(w==0))/length(w)
load(paste("results_tau",tau,"_rho",rho,"_p",p,"_prob",prob,".dat",sep=""))
errors = tests.x = tests.y = tests2.x = tests2.y = numeric(0)
for(r in 1:rep){
  errors = rbind(errors, results[[r]]$error)
  tests.x = rbind(tests.x, results[[r]]$dz.1)
  tests.y = rbind(tests.y, results[[r]]$dz.2)
  tests2.x = rbind(tests2.x, results[[r]]$db2.1)
  tests2.y = rbind(tests2.y, results[[r]]$db2.2)
}
round(mean(abs(errors[,2])/abs(errors[,1])),4)
round(median(abs(errors[,2])/abs(errors[,1])),4)
round(mean(abs(errors[,4])/abs(errors[,3])),4)
round(median(abs(errors[,4])/abs(errors[,3])),4)
round(mean(abs(errors[,6])/abs(errors[,5])),4)
round(median(abs(errors[,6])/abs(errors[,5])),4)

round(mean(colMeans(tests2.x[,which(beta[1:(m+1)]==0)])),2)
round(mean(colMeans(tests.x[,which(beta[1:(m+1)]==0)])),2)
round(mean(colMeans(tests2.x[,which(beta[1:(m+1)]!=0)])),2)
round(mean(colMeans(tests.x[,which(beta[1:(m+1)]!=0)])),2)
round(mean(colMeans(tests2.y[,which(beta[idx[1:50]]==0)])),2)
round(mean(colMeans(tests.y[,which(beta[idx[1:50]]==0)])),2)
round(mean(colMeans(tests2.y[,which(beta[idx[1:50]]!=0)])),2)
round(mean(colMeans(tests.y[,which(beta[idx[1:50]]!=0)])),2)
