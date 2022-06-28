####################################
# Function to Find Nash Equilibria #
####################################



# libraries ---------------------------------------------------------------

library(tidyverse)
library(foreach)


# best response functions -------------------------------------------------

xi <- function(
  Ri = 1/2,
  sigix = 1/2,
  etax = 1/2,
  etay = 1/2
) {
  
  # Specify intercept and slope parameters
  
  delta0 <- sigix * etay 
  delta1 <- sigix - sigix * etay
  delta2 <- sigix * (etax - etay) - etax
  
  # The best responses
  
  yi <- xi <- 0
  xj <- seq(0, 1 - Ri, by = 0.01)
  for(i in 1:length(xj)) {
    xi[i] <- (delta0 + delta1 * Ri + delta2 * xj[i]) 
    xi[i] <- min(max(xi[i], 0), Ri)
    yi[i] <- Ri - xi[i]
  }
  
  # Return i's best responses
  return(
    tibble(
      recipient = 
        rep(c("Recipient x", "Recipient y"), each = length(xj)),
      aidi = c(xi, yi),
      aidj = c(xj, (1 - Ri) - xj)
    )
  )
}

xj <- function(
  Rj = 1/2,
  sigjx = 1/2,
  etax = 1/2,
  etay = 1/2
) {
  
  # Specify intercept and slope parameters
  
  delta0 <- sigjx * etay 
  delta1 <- sigjx - sigjx * etay
  delta2 <- sigjx * (etax - etay) - etax
  
  # The best responses
  yj <- xj <- 0
  xi <- seq(0, 1 - Ri, by = 0.01)
  for(i in 1:length(xi)) {
    xj[i] <- (delta0 + delta1 * Rj + delta2 * xi[i]) 
    xj[i] <- min(max(xj[i], 0), Rj[i])
    yj[i] <- Rj - xj[i]
  }
  
  # Return j's best responses
  return(
    tibble(
      recipient = 
        rep(c("Recipient x", "Recipient y"), each = length(xi)),
      aidj = c(xj, yj),
      aidi = c(xi, (1 - Rj) - xi)
    )
  )
}

# xi() # test
# xj() # test


# Reaction Paths ----------------------------------------------------------

react_paths <- function(
  sigix = 1/2,
  sigjx = 1/2,
  etax = 1/2,
  etay = 1/2
) {
  
  # Reaction slopes
  deltai2 <- sigix * (etax - etay) - etax
  deltaj2 <- sigjx * (etax - etay) - etax
  
  # Return
  return(
    tibble(
      Donor = c("Donor i", "Donor j"),
      Reaction = c(deltai2, deltaj2)
    )
  )
}

# react_paths() # test


# Nash Equilibria ---------------------------------------------------------

nash_eq <- function(
  Ri = 1/2,
  sigix = 1/2,
  sigjx = 1/2,
  etax = 1/2,
  etay = 1/2
) {
  
  # Parameter identities for Donor i
  deltai0 <- sigix * etay 
  deltai1 <- sigix - sigix * etay
  deltai2 <- sigix * (etax - etay) - etax
  
  # Parameter identities for Donor j
  deltaj0 <- sigjx * etay 
  deltaj1 <- sigjx - sigjx * etay
  deltaj2 <- sigjx * (etax - etay) - etax
  
  # Specify Donor j's share of resources
  Rj <- 1 - Ri
  
  # Equilibrium outcomes
  xistar <-
    ((deltai0 + deltai1 * Ri + deltai2 * (deltaj0 + deltaj1 * Rj)) /
    (1 - deltai2 * deltaj2))
  xistar <- min(max(xistar, 0), Ri)
  yistar <- Ri - xistar
  xjstar <- 
    ((deltaj0 + deltaj1 * Rj + deltaj2 * (deltai0 + deltai1 * Ri)) /
    (1 - deltai2 * deltaj2)) 
  xjstar <- min(max(xjstar, 0), Rj)
  yjstar <- Rj - xjstar
  
  # Return equilibrium outcomes
  return(
    tibble(
      donor = c("i", "j"),
      x = c(xistar, xjstar),
      y = c(yistar, yjstar)
    )
  )
}

# nash_eq() # test
# nash_eq(
#   sigix = 3/4,
#   sigjx = 1/4
# )


# Collective Utility ------------------------------------------------------

collective_utility <- function(
  alpha = 1/2,
  sigix = 1/2,
  sigjx = 1/2,
  etax = 1/2,
  etay = 1/2,
  Ri = 1/2,
  bins = 100
) {
  
  # Specify the parameters
  Rj <- 1 - Ri
  xi <- seq(0, Ri, len = bins)
  xj <- seq(0, Rj, len = bins)
  U <- ui <- uj <- 
    matrix(0, length(xi), length(xj))
  for(i in 1:nrow(U)) {
    for(j in 1:ncol(U)) {
      
      # Realized objectives
      Xi <- xi[i] + etax*xj[j]
      Xj <- xj[j] + etax*xi[i]
      yi <- (Ri - xi[i])
      yj <- (Rj - xj[j])
      Yi <- yi + etay*yj
      Yj <- yj + etay*yi
      
      # Utilities for Donors i and j
      ui[i, j] <- (Xi^sigix)*(Yi^(1-sigix))
      uj[i, j] <- (Xj^sigjx)*(Yj^(1-sigjx))
      
      # Collective utility
      U[i, j] <- alpha*ui[i, j] + (1-alpha)*uj[i, j]
    }
  }
  
  rownames(U) <- rownames(uj) <- rownames(ui) <-
    xi
  colnames(U) <- colnames(uj) <- colnames(ui) <-
    xj

  # Return individual and collective utilities
  return(list(U=U,ui=ui,uj=uj))
}

prep_for_contour <- function(
  X
) {
  X %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    gather(key, value, -rowname) %>%
    mutate(
      key = as.numeric(gsub("V","", key)),
      rowname = as.numeric(rowname)
    ) %>%
    rename(
      donor_i = rowname,
      donor_j = key,
      utility = value
    ) -> X
  return(X)
}

# collective_utility()$U %>%
#   prep_for_contour() %>%
#   ggplot() +
#   geom_contour(
#     aes(donor_i, donor_j, z = utility),
#     bins = 20
#   )

find_optimum <- function(
  Ri = 1/2,
  sigix = 1/2,
  sigjx = 1/2,
  etax = 1/2,
  etay = 1/2
) {
  objective <- function(
    Ri,
    sigix,
    sigjx,
    etax,
    etay,
    pars
  ) {
    Rj <- 1 - Ri
    xi <- pars[1]
    xj <- pars[2]
    yi <- Ri - xi
    yj <- Rj - xj
    
    Xi <- max(xi + etax * xj, 0)
    Xj <- max(xj + etax * xi, 0)
    Yi <- max(yi + etay * yj, 0)
    Yj <- max(yj + etay * yi, 0)
    
    ui <- (Xi^sigix) *
      (Yi^(1 - sigix))
    uj <- (Xj^sigjx) *
      (Yj^(1 - sigjx))
    U <- (ui + uj) +
      -999999999 * any(c(Xi, Xj, Yi, Yj) <= 0)
    return(-U)
  }
  opt <- optim(
    par = c(Ri, 1 - Ri) / 2,
    fn = objective,
    Ri = Ri,
    sigix = sigix,
    sigjx = sigjx,
    etax = etax,
    etay = etay,
    method = "L-BFGS-B",
    lower = 0, upper = c(Ri, 1 - Ri)
  )
  return(opt)
}

# find_optimum()$par
# find_optimum(
#   etax = -1/2
# )$par

# Old Code ----------------------------------------------------------------

# # libraries ---------------------------------------------------------------
# 
# library(tidyverse)
# library(foreach)
# 
# 
# # utility function --------------------------------------------------------
# 
# par_obj <- function(
#   alpha = 1/2,
#   X,
#   sigma,
#   eta,
#   Ri
# ) {
#   Rj <- 1 - Ri
#   xi <- X[1]
#   xj <- X[2]
#   Xi <- xi + eta[1]*xj
#   Xj <- xj + eta[1]*xi
#   yi <- (Ri - xi)
#   yj <- (Rj - xj)
#   Yi <- yi + eta[2]*yj
#   Yj <- yj + eta[2]*yi
#   ui <- (Xi^sigma[1])*(Yi^(1-sigma[1]))
#   uj <- (Xj^sigma[2])*(Yj^(1-sigma[2]))
#   U  <- alpha*ui + (1-alpha)*uj
#   
#   return(tibble(U=U,ui=ui,uj=uj))
# }
# 
# 
# 
# # the best response function ----------------------------------------------
# 
# xi <- function(
#   Ri = 1/2,
#   sigma = 1/2,
#   eta = c(1,1)
# ) {
#   Rj <- 1 - Ri
#   x <- seq(0,Rj,len=1000)
#   xi <- 
#     sigma*(Ri + eta[1]*x + eta[2]*(Rj - x)) - eta[1]*x
#   xi <-
#     ifelse(xi > Ri, Ri, xi)
#   xi <-
#     ifelse(xi < 0, 0, xi)
#   return(
#     tibble(x, xi)
#   )
# }
# 
# 
# # function to find the intersection of reaction paths ---------------------
# 
# my_int <- function (curve1, curve2, empirical = TRUE, domain = NULL,
#                     ties) 
# {
#   if (!empirical & missing(domain)) {
#     stop("'domain' must be provided with non-empirical curves")
#   }
#   if (!empirical & (length(domain) != 2 | !is.numeric(domain))) {
#     stop("'domain' must be a two-value numeric vector, like c(0, 10)")
#   }
#   if (empirical) {
#     curve1_f <- approxfun(curve1$x, curve1$y, rule = 2,
#                           ties = ties)
#     curve2_f <- approxfun(curve2$x, curve2$y, rule = 2,
#                           ties = ties)
#     point_x <- uniroot(function(x) curve1_f(x) - curve2_f(x), 
#                        0:1)$root
#     point_y <- curve2_f(point_x)
#   }
#   else {
#     point_x <- uniroot(function(x) curve1(x) - curve2(x), 
#                        domain)$root
#     point_y <- curve2(point_x)
#   }
#   return(list(x = point_x, y = point_y))
# }
# 
# 
# # function to plot reaction paths and equilibria --------------------------
# 
# plot_eq <- function(
#   Ri, sigma, eta, ties = min, bins = 2
# ) {
#   bind_cols(
#     xi(Ri = Ri, sigma = sigma[1], eta = eta) %>%
#       rename(x1 = x),
#     xi(Ri = 1 - Ri, sigma = sigma[2], eta = eta) %>%
#       rename(xj = xi, x2 = x)
#   ) -> df
#   
#   ggplot(df) +
#     geom_point(
#       aes(xi,x1),col='blue',alpha=.05
#     ) +
#     geom_point(
#       aes(x2,xj),col='red',alpha=.05
#     ) -> g
#   
#   xi <- seq(0,Ri,len=20)
#   xj <- seq(0,1-Ri,len=20)
#   foreach(i = 1:length(xi),.combine = 'rbind') %do% {
#     foreach(j = 1:length(xj),.combine = 'rbind') %do% {
#       par_obj(
#         X = c(xi[i],xj[j]),
#         sigma = sigma,
#         eta = eta,
#         Ri = Ri
#       ) %>%
#         mutate(
#           xi = xi[i],
#           xj = xj[j]
#         )
#     }
#   } -> U
#   
#   
#   
#   # curve1 <- df %>% select(xi,x1) %>%
#   #   rename(x = xi, y = x1)
#   # curve2 <- df %>% select(x2,xj) %>%
#   #   rename(x = x2, y = xj) 
#   # if(is.character(try(my_int(curve1,curve2,ties=ties)))){
#   #   stop("Pick a different method for dealing with ties!")
#   # } 
#   #int <- my_int(curve1,curve2,ties=ties)
#   
#   g +
#     # geom_point(
#     #   aes(int$x,int$y),
#     #   size = 3
#     # ) +
#     labs(
#       x = expression(x[1]),
#       y = expression(x[2])
#     ) +
#     scale_x_continuous(
#       breaks = c(0,sigma[1]*Ri,Ri),
#       labels = c('0',
#                  "DT",
#                  expression(R[i]))
#     ) +
#     scale_y_continuous(
#       breaks = c(0,sigma[2]*(1-Ri),1-Ri),
#       labels = c('0',
#                  "DT",
#                  expression(R[j]))
#     ) +
#     geom_vline(
#       xintercept = c(sigma[1]*Ri,Ri),
#       lty = 3
#     ) +
#     geom_hline(
#       yintercept = c(sigma[2]*(1-Ri), 1 - Ri),
#       lty = 3
#     ) +
#     geom_abline(
#       intercept = 0,
#       slope = ((1-Ri))/(Ri),
#       lty = 2
#     ) +
#     # geom_contour(
#     #   data = U,
#     #   aes(xi,xj,z=ui),
#     #   bins = bins,
#     #   col='blue'
#     # ) +
#     # geom_contour(
#     #   data = U,
#     #   aes(xi,xj,z=uj),
#     #   bins = bins,
#     #   col='red'
#     # ) +
#     geom_contour(
#       data = U,
#       aes(xi,xj,z=U),
#       bins = bins,
#       col='black',
#       alpha = 0.8
#     ) +
#     # geom_segment(
#     #   aes(
#     #     x = int$x,
#     #     xend = int$x,
#     #     y = 0,
#     #     yend = int$y
#     #   )
#     # ) +
#     # geom_segment(
#     #   aes(
#     #     x = 0,
#     #     xend = int$x,
#     #     y = int$y,
#     #     yend = int$y
#     #   )
#     # ) +
#     ggridges::theme_ridges(
#       center_axis_labels = T
#     ) +
#     theme(
#       axis.title.y = element_text(angle = 0, vjust = .5),
#       axis.text = element_text(color='black')
#     )
# }
# 
# # Ri <- 3/4
# # sigma <- c(0.25,0.75)
# # eta <- c(-1/2,1/2)
# # 
# # plot_eq(
# #   Ri = Ri, sigma = sigma, 
# #   eta = eta, ties = 
# # ) +
# #   labs(
# #     title = "Reaction Paths for Objective X"
# #   ) 
# 
# 
# 
# # the generalized model ---------------------------------------------------
# 
# # Number of actors and objectives
# # N <- 2 # actors
# # M <- 2 # objectives
# # 
# # x <- matrix(0,nrow=N,ncol=M)
# # colnames(x) <- paste('dim.',1:M)
# # rownames(x) <- paste('state',1:N)
# # 
# # sig <- matrix(0,nrow=N,ncol=M)
# # sig <- apply(sig,c(1,2),function(x) runif(1,0,1))
# # sig <- sig/rowSums(sig)
# # colnames(sig) <- colnames(x)
# # rownames(sig) <- rownames(x)
# # 
# # Leta <- list()
# # for(i in 1:M) {
# #   eta <- matrix(0,nrow=N,ncol=N)
# #   colnames(eta) <- rownames(x)
# #   rownames(eta) <- rownames(x)
# #   eta[upper.tri(eta)] <- runif(choose(N,2),-1,1)/(N-1)
# #   eta[lower.tri(eta)] <- t(eta)[lower.tri(eta)]
# #   Leta[[i]] <- eta
# #   rm(eta)
# # }
# # 
# # R <- runif(N)
# # R <- R/sum(R)
# # 
# # sig[,1]*(Leta[[1]] - Leta[[M]]) - Leta[[1]]*((1:M)[1]==1)
# # sig[,1]*(Leta[[2]] - Leta[[M]]) - Leta[[1]]*((1:M)[1]==2)
# # LJx <- list()
# # for(n in 1:(M-1)) {
# #   Jx <- list()
# #   for(m in 1:(M-1)) {
# #     Jx[[m]] <- 
# #       sig[,m]*(Leta[[m]] - Leta[[M]]) -
# #       Leta[[m]]*(n==m)
# #   }
# #   LJx[[n]] <- Jx
# # }
# # 
# # Jr <- list()
# # for(m in 1:(M-1)) {
# #   Leta1 <- Leta[[m]]
# #   diag(Leta1) <- 1
# #   Jr[[m]] <- sig[,m]*Leta1
# #   rm(Leta1)
# # }
# # 
# # # check whether convergence is guaranteed:
# # spec <- 0
# # for(m in 1:(M-1)) {
# #   spec[m] <- max(abs(eigen(LJx[[m]][[1]])$values))
# # }
# # spec
# # 
# # x <- sig*R
# # diff <- 0
# # for(i in 1:1000) {
# #   if(i>1) {
# #     if(diff[i-1]==0) {
# #       break
# #     }
# #   }
# #   x_nm1 <- x#/rowSums(x) * R
# #   for(m in 1:(M-1)) {
# #     x[,m] <-
# #       cbind(do.call(cbind,LJx[[m]]),Jr[[m]]) %*%
# #       cbind(c(x_nm1[,-M][is.numeric(x_nm1)],R))
# #   }
# #   x <- apply(x,c(1,2),function(x) max(x,0))
# #   x <- x*(x<R) + R*(x>=R)
# #   x[,M] <- cbind(R) - if(ncol(x)==2) x[,-M] else rowSums(x[,-M])
# #   x[,M] <- x[,M]*(x[,M]>0)
# #   #x_fin <- x
# #   #x <- x/rowSums(x) * R
# #   diff[i] <- round(sum(abs(x - x_nm1)),4)
# # }
# # 
# # plot(1:length(diff),diff,
# #      type = 'b',
# #      pch=19,
# #      xlab = 'iteration',
# #      ylab = expression(
# #        Sigma*abs(x[n] - x[n-1])
# #      ),
# #      main = paste(i - 1, 'iterations to convergence'))
# # 
# # # check that expenditures equal resources
# # sum(x)
# # sum(rowSums(x) - R)
# # diff[length(diff)]
# # max(spec)
# # x
# # 
# # 
# # # the generalized model redux ---------------------------------------------
# # 
# # # demand function
# # x <- function(R,sig,eta,x) {
# #   d <- sig*(R + sum(eta*x)) - 
# #   min(max())
# # }
