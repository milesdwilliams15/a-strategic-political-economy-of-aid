#####################################
# Code to Identify Pareto Solutions #
#####################################


# libraries ---------------------------------------------------------------

library(tidyverse)
library(foreach)


# pareto objective --------------------------------------------------------

par_obj <- function(
  alpha = 1/2,
  X,
  sigma,
  eta,
  Ri,
  c
) {
  Rj <- 1 - Ri
  xi <- X[1]
  xj <- X[2]
  Xi <- xi + eta[1]*xj
  Xj <- xj + eta[1]*xi
  yi <- (Ri - c[1]*xi)/c[2]
  yj <- (Rj - c[3]*xj)/c[4]
  Yi <- yi + eta[2]*yj
  Yj <- yj + eta[2]*yi
  ui <- (Xi^sigma[1])*(Yi^(1-sigma[1]))
  uj <- (Xj^sigma[2])*(Yj^(1-sigma[2]))
  U  <- alpha*ui + (1-alpha)*uj
  
  return(tibble(U=U,ui=ui,uj=uj))
}



# find pareto solution ----------------------------------------------------

Ri <- 1/2
xi <- seq(0,Ri,len=100)
xj <- seq(0,1-Ri,len=100)
foreach(i = 1:length(xi),.combine = 'rbind',.verbose = T) %do% {
  foreach(j = 1:length(xj),.combine = 'rbind') %do% {
    par_obj(
      X = c(xi[i],xj[j]),
      sigma = c(0.25,0.75),
      eta = c(-1/2,1/2),
      Ri = Ri,
      c = c(1,1,1,1)
    ) %>%
      mutate(
        xi = xi[i],
        xj = xj[j]
      )
  }
} -> U

# pareto optimal solution
U %>% filter(U==max(U,na.rm=T))

# best solution for i
U %>% na.omit %>% filter(ui==max(ui)) %>%
  filter(U==max(U))

# best solution for j
U %>% na.omit %>% filter(uj==max(uj)) %>%
  filter(U==max(U))

ggplot(U) +
  geom_contour(
    aes(xi,xj,z=ui),
    bins = 100,
    col='blue'
  ) +
  geom_contour(
    aes(xi,xj,z=uj),
    bins = 100,
    col='red'
  ) +
  geom_contour(
    aes(xi,xj,z=U),
    binwidth = 0.01,
    col='black'
  )

ggplot(U) +
  # geom_point(
  #   aes(ui,uj)
  # ) +
  geom_contour(
    aes(ui,uj,z=xj),
    color='blue'
  )

alphas <- seq(0,1,len=100)
foreach(
  i = 1:length(alphas),
  .combine = 'rbind',
  .verbose = T
) %do% {
  
}

optim(
  par = c(0,0),
  fn = par_obj,
  alpha = 1/2,
  sigma = c(.25,.75),
  eta = c(1/2,1/2),
  Ri = 1/2,
  c = c(1,1,1,1),
  method = "L-BFGS-B",
  lower = 0.0001,
  upper = c(Ri-.0001,1-Ri-.0001)
)


X <- seq(0,1,.01)
Y <- X
Ri <- 1/2
Rj <- 1 - Ri
etax <- 1
etay <- 1
ui <- matrix(nrow=length(X),ncol=length(X))
uj <- ui
# for(i in 1:nrow(ui)) {
#   for(j in 1:nrow(uj)) {
#     ui[i,j] <-
#       0.25*log(min(x[i],Ri) + etax*min(x[j],Rj) + .00001) +
#       0.75*log(Ri - min(x[i],Ri) + etay*(Rj - min(x[j],Rj))+ .00001)
#     uj[i,j] <-
#       0.75*log(min(x[j],Rj) + etax*min(x[i],Ri)+ .00001) +
#       0.25*log(Rj - min(x[j],Rj) + etay*(Ri - min(x[i],Ri))+ .00001)
#     #X[i,j] <- 
#   }
# }
for(i in 1:nrow(ui)) {
  for(j in 1:nrow(uj)) {
    ui[i,j] <-
      0.25*log(X[j]) + 0.75*log(Y[i])
    uj[i,j] <-
      0.75*log(X[j]) + 0.25*log(Y[i])
  }
}
par(
  bty = 'l',
  las = 1
)
contour(ui,col='blue',nlevels=10,lty=2,
        xlab = 'X',ylab='Y')
contour(uj,add=T,col='red',nlevels=10,lty=2)
contour((ui+uj)/2,add=T,nlevels=10,lty=2)
abline(
  a = 0,
  b = 0.25/0.75,
  col = 'blue'
)
abline(
  a = 0,
  b = 0.75/0.25,
  col = 'red'
)
abline(
  a = 0,
  b = 1
)
lines(c(0,1),c(1,0))
points(0.75,0.25,pch=19,col='blue')
points(0.25,0.75,pch=19,col='red')
points(0.5,0.5,pch=19)