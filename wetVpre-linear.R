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
# work through grid

missing <- -9999.

obs       <- array(data=-32768,dim=c(xlen,ylen))
slope     <- array(data=missing,dim=c(xlen,ylen))
intercept <- array(data=missing,dim=c(xlen,ylen))
r2        <- array(data=missing,dim=c(xlen,ylen))
rmse      <- array(data=missing,dim=c(xlen,ylen))

# ----------------------
# regression on all data

for (y in seq(1,ylen)) {

  cat(" working on row",y,"of",ylen,"\r")
  
  for (x in seq(1,xlen-1)) {
  
    if (!is.na(pre[x,y,1])) {
  
      ind_all <- pre[x,y,]
      dep_all <- wet[x,y,]
      
      didrain <- ind_all > 0.2
      
      obs[x,y] <- sum(didrain, na.rm = TRUE)
      
      if (obs[x,y] < 10) next

      ind <- ind_all[didrain]
      dep <- dep_all[didrain]
      
      lmfit <- lm(dep ~ log(ind))
      
      pval <- summary(lmfit)$coefficients[2,4]
      
      if (pval > 0.01) next  # do not record stats for fits that don't meet a 99% confidence threshold
      
      coefs <- coef(lmfit)

      slope[x,y]     <- coefs[2]
      intercept[x,y] <- coefs[1]
      
      if (slope[x,y] < 0) cat(x,y,slope[x,y],"\n")
      
      r2[x,y] <- summary(lmfit)$adj.r.squared

      rmse[x,y] <- RMSE(lmfit)

    } 
  }
}

cat("\n")
cat("writing annual\n")

outfile <- "wetVpre-linear-annual.nc"

ofid <- nc_open(outfile,write=TRUE)

# ----
ncvar_put(ofid,"obs",obs)
actual_range <- c(min(obs[obs != -32768]),max(obs[obs != -32768]))
ncatt_put(ofid,"obs","actual_range",actual_range,prec="short")

# ----
ncvar_put(ofid,"slope",slope)
actual_range <- c(min(slope[slope != missing]),max(slope[slope != missing]))
ncatt_put(ofid,"slope","actual_range",actual_range,prec="float")

# ----
ncvar_put(ofid,"intercept",intercept)
actual_range <- c(min(intercept[intercept != missing]),max(intercept[intercept != missing]))
ncatt_put(ofid,"intercept","actual_range",actual_range,prec="float")

# ----
ncvar_put(ofid,"r2",r2)
actual_range <- c(min(r2[r2 != missing]),max(r2[r2 != missing]))
ncatt_put(ofid,"r2","actual_range",actual_range,prec="float")

# ----
ncvar_put(ofid,"rmse",rmse)
actual_range <- c(min(rmse[rmse != missing]),max(rmse[rmse != missing]))
ncatt_put(ofid,"rmse","actual_range",actual_range,prec="float")

nc_close(ofid)






# ----------------------
# regressions per month

# outfile <- "wetVpre-linear.nc"
# 
# ofid <- nc_open(outfile,write=TRUE)
# 
# for (m in seq(1,12)) {
#   for (y in seq(1,ylen)) {
#   
#     cat(" working on row",y,"of",ylen,"\r")
#     
#     for (x in seq(1,xlen-1)) {
#     
#       if (!is.na(pre[x,y,1])) {
#     
#         ind_all <- pre[x,y,seq(m,360,12)]
#         dep_all <- wet[x,y,seq(m,360,12)]
#         
#         didrain <- ind_all > 0.2
#         
#         obs[x,y] <- sum(didrain, na.rm = TRUE)
#         
#         if (obs[x,y] < 10) next
#   
#         ind <- ind_all[didrain]
#         dep <- dep_all[didrain]
#         
#         lmfit <- lm(dep ~ log(ind))
#   
#         coefs <- coef(lmfit)
#   
#         slope[x,y]     <- coefs[2]
#         intercept[x,y] <- coefs[1]
#         
#         r2[x,y] <- summary(lmfit)$adj.r.squared
#   
#         rmse[x,y] <- RMSE(lmfit)
#   
#       } 
#     }
#   }
#   
#   cat("\n")
#   cat("writing monthly",m,"\n")
#   
#   ncvar_put(ofid,"obs",obs,start=c(1,1,m),count=c(xlen,ylen,1))
#   ncvar_put(ofid,"slope",slope,start=c(1,1,m),count=c(xlen,ylen,1))
#   ncvar_put(ofid,"intercept",intercept,start=c(1,1,m),count=c(xlen,ylen,1))
#   ncvar_put(ofid,"r2",r2,start=c(1,1,m),count=c(xlen,ylen,1))
#   ncvar_put(ofid,"rmse",rmse,start=c(1,1,m),count=c(xlen,ylen,1))
# 
# }

# nc_close(ofid)

# ----------------------
