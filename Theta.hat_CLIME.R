inv_CLIME <- function(Sigma.hat, mus) {
  #################################
  # lpSolve method
  
  # px <- ncol(Sigma.hat)
  # c1 <- lapply(1:px,
  #              function(j) { a <- rep(0,3*px); a[j] <- -1; a[j+px] <- 1; a[j+2*px] <- -1; a } ) %>%
  #   reduce(rbind)
  # c2 <- lapply(1:px,
  #              function(j) { a <- rep(0,3*px); a[j] <- 1;a[j+px] <- -1; a[j+2*px] <- -1; a } ) %>%
  #   reduce(rbind)
  # c3 <- lapply(1:px,
  #              function(j) { c(-1*Sigma.hat[j,], Sigma.hat[j,], rep(0, px)) }) %>%
  #   reduce(rbind)
  # c4 <- lapply(1:px,
  #              function(j) { c(Sigma.hat[j,], -1*Sigma.hat[j,], rep(0, px)) }) %>%
  #   reduce(rbind)
  # A <- rbind(c1, c2, c3, c4)
  # dir <- rep("<=", 4*px)
  # obj <- c(rep(0,px), rep(0,px), rep(1,px))
  #
  # Theta.hat <- lapply(
  #   1:px,
  #   function(j) {
  #     mu <- mus[j]
  #     b1 <- rep(mu, px); b1[j] = mu-1
  #     b2 <- rep(mu, px); b2[j] = mu+1
  #     b <- c(rep(0, 2*px), b1, b2)
  #     res_j <- lp(direction="min", objective.in=obj, const.mat=A, const.dir=dir, const.rhs=b)
  #     theta_j <- res_j$solution[1:px] - res_j$solution[(px+1):(2*px)]
  #     theta_j
  #   }
  # ) %>%
  #   reduce(rbind)
  # Theta.hat
  
  #################################
  # Mosek method
  
  px <- ncol(Sigma.hat)
  Id <- diag(1, px)
  zeros <- matrix(ncol=px,nrow=px); zeros[,] <- 0
  Theta.hat <- matrix(ncol=px,nrow=px)
  
  # constraint matrix
  A1 <- cbind(Id, -1*Id)
  A2 <- cbind(-1*Id, -1*Id)
  A3 <- cbind(Sigma.hat, zeros)
  A4 <- cbind(-1*Sigma.hat, zeros)
  A <- rbind(A1, A2, A3, A4) %>%
    Matrix(sparse=T)
  
  # lower and upper bounds on variables
  blx <- c(rep(-Inf, px), rep(0, px))
  bux <- rep(Inf, 2*px)
  
  objective <- c(rep(0,px), rep(1,px))
  
  for ( j in 1:px ) {
    mu <- mus[j]
    
    # lower and upper bounds on constraints
    blc <- rep(-Inf, 4*px)
    buc3 <- rep(mu, px); buc3[j] <- mu+1
    buc4 <- rep(mu, px); buc4[j] <- mu-1
    buc <- c(rep(0, 2*px), buc3, buc4)
    
    mosek.prob <- list(
      sense="min",
      c=objective,
      A=A,
      bc=rbind(blc, buc),
      bx=rbind(blx, bux)
    )
    mosek.res <- try(mosek(mosek.prob, list(verbose=0)), silent=TRUE)
    Theta.hat[j,] <- mosek.res$sol$bas$xx[1:px]
  }
  Theta.hat
}
