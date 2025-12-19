# -------------------------------------------------------------------------------
# R script to calculate nonlinear relationship between precip and wetdays

library(ncdf4)
library(minpack.lm)
library(qpcR)

eps <- .Machine$double.xmin

# -------------------------------------------------------------------------------
# file names

args <- commandArgs(trailingOnly <- TRUE)

mon <- 1 # as.numeric(args[1])

# open the input file

ifid <- nc_open("/Users/jkaplan/nosync/datasets/wetfrac30m_t0.2_202512.nc",write=FALSE)

# get dimension lengths

lon  <- ncvar_get(ifid,"lon")
lat  <- ncvar_get(ifid,"lat")
time <- ncvar_get(ifid,"time")

xlen <- length(lon)
ylen <- length(lat)
tlen <- length(time)

# get input data

pre <- ncvar_get(ifid,varid="pre")
wet <- ncvar_get(ifid,varid="wetf")

# close input files

nc_close(ifid)

# ----------------------
# calculate single pixel

# southern amazon

# x <- 238
# y <- 152

# california

# x <- 144
# y <- 249

# central amazon

# x <- 221
# y <- 174

# spain

# x <- 346
# y <- 255

# pixels with negative slope and significant p-value

# x <- 215
# y <- 129

# x <- 357
# y <- 230

# x <- 384
# y <- 232

x <- 194
y <- 233

# ---

t0 <- 1
t1 <- tlen

# ----------------------

ind_all <- pre[x,y,t0:t1]
dep_all <- wet[x,y,t0:t1]

didrain <- ind_all > 0.2

obs <- sum(didrain, na.rm = TRUE)

ind <- ind_all[didrain]
dep <- dep_all[didrain]

lmfit <- lm(dep ~ log(ind))

# rlmfit <- rlm(dep ~ log(ind))

coefs <- coef(lmfit)

# rng <- max(ind) - min(ind)
# 
# ind_same <- abs(max(ind) - min(ind)) == 0
# dep_same <- abs(max(dep) - min(dep)) == 0
# 
# fit <- try(nlsLM(dep ~ (1-exp(-a*indr))^b,start=list(a=4,b=1),lower=c(eps,eps),upper=c(100,3),control = nls.lm.control(maxiter = 100)),TRUE)
# 
# coefs <- coef(fit)

slope <- coefs[2]
intercept <- coefs[1]

rmse <- RMSE(lmfit)

# ----------------------
# plot fit lines as a diagnostic

# two plots per page

par(mfrow=c(1,2))

plot(dep ~ ind, pch=19, cex=0.5, xlab="precipitation", ylab="wet fraction")

xvals <- seq(0,max(ind),length.out = 51)

lines(xvals,slope * log(xvals) + intercept,col="red")

plot(dep ~ log(ind), pch=19, cex=0.5, xlab="log(precipitation)", ylab="wet fraction")

abline(lmfit,col="red")

# ----------------------
