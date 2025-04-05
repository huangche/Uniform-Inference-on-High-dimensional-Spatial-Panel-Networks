find_mus <- function(Sigma.hat) {
  ##########################
  # Using lpSolve
  #
  # px <- ncol(Sigma.hat)
  # Id <- diag(1, px)
  # zeros <- matrix(ncol=px,nrow=px); zeros[,] <- 0
  #
  # c1 <- cbind(Sigma.hat, -1*Sigma.hat, -1*Id, rep(0, px))
  # c2 <- cbind(-1*Sigma.hat, Sigma.hat, -1*Id, rep(0, px))
  # c3 <- cbind(zeros, zeros, diag(1, px), rep(-1, px))
  # A <- rbind(c1, c2, c3)
  #
  # obj <- c(rep(0, 3*px), 1)
  # dir <- rep("<=", 3*px)
  # mus <- lapply(
  #   1:px,
  #   function(j) {
  #     b1.j <- rep(0,px); b1.j[j] <- -1
  #     b2.j <- rep(0,px); b2.j[j] <- 1
  #     b3.j <- rep(0,px)
  #     b.j <- c(b1.j, b2.j, b3.j)
  #     res.j <- lp(direction="min", objective.in=obj, const.mat=A, const.dir=dir, const.rhs=b.j)
  #     res.j$objval
  #   }
  # ) %>% as.numeric
  # mus
  
  ##########################
  # Using Mosek
  
  px <- ncol(Sigma.hat)
  Id <- diag(1, px)
  zeros <- matrix(ncol=px,nrow=px); zeros[,] <- 0
  mus <- numeric(px)
  
  # constraint marix A
  A1 <- cbind(Sigma.hat, -Id, rep(0,px))
  A2 <- cbind(-Sigma.hat, -Id, rep(0,px))
  A3 <- cbind(zeros, Id, rep(-1, px))
  A <- rbind(A1, A2, A3) %>%
    Matrix(sparse=T)
  
  # lower bounds on variables
  blx <- c(rep(-Inf, px), rep(0, px), 0)
  bux <- c(rep(Inf, 2*px), Inf)
  
  objective <- c(rep(0,2*px), 1)
  
  for ( j in 1:px ){
    blc <- rep(-Inf, 3*px)
    buc1 <- rep(0, px); buc1[j] <- 1
    buc2 <- rep(0, px); buc2[j] <- -1
    buc3 <- rep(0, px)
    buc <- c(buc1, buc2, buc3)
    
    mosek.prob <- list(
      sense="min",
      c=objective,
      A=A,
      bc=rbind(blc, buc),
      bx=rbind(blx, bux)
    )
    mosek.res <- try(mosek(mosek.prob, list(verbose=0, soldetail=1)), silent=TRUE)
    mus[j] <- mosek.res$sol$bas$pobjval
  }
  mus
}

