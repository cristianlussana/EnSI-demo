#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(RANN))
suppressPackageStartupMessages(library(sp))
suppressPackageStartupMessages(library(raster))
#options(warn=2)
# 
#..............................................................................
# Functions

#+ manage fatal error
boom <- function( str=NA, code=NA, status=1) {
  cat("Fatal Error ")
  if ( !is.na(code)) {
    if ( code == 1) {
      cat("file not found ")
    } else if ( code == 2) {
      cat("input parameter not found ")
    }
  }
  if ( !is.na(str)) cat( str)
  cat("\n")
  quit( status= status )
}

#+ the end 
rip <- function( str=NA, code=NA, t0=NA, status=0) {
  cat( "the End : ")
  if ( !is.na(code) ) {
    if ( code == 0 ) cat( "normal exit : ")
  }
  if ( !is.na(t0)) {
    t1 <- Sys.time()
    cat( paste( "total time=", round(t1-t0,1), attr(t1-t0,"unit")))
  }
  if ( !is.na(str)) cat( str)
  cat("\n")
  quit( status= status )
}

# + replace elements of a string with date-time elements
replaceDate <- function( string=NULL,
                         date.str=NULL,
                         var.str=NULL,
                         year_string="%Y",
                         month_string="%m",
                         day_string="%d",
                         hour_string="%H",
                         min_string="%M",
                         sec_string="%S",
                         var_string="%V",
                         format="%Y-%m-%d %H:%M:%S") {
#------------------------------------------------------------------------------
  if (is.null(string) | is.null(date.str)) return(NULL)
  Rdate<-as.POSIXlt(str2Rdate(date.str,format=format))
  yyyy<-Rdate$year+1900
  mm<-formatC(Rdate$mon+1,width=2,flag="0")
  dd<-formatC(Rdate$mday,width=2,flag="0")
  hh<-formatC(Rdate$hour,width=2,flag="0")
  MM<-formatC(Rdate$min,width=2,flag="0")
  SS<-formatC(Rdate$sec,width=2,flag="0")
  out<-gsub(year_string,yyyy,string)
  out<-gsub(month_string,formatC(mm,width=2,flag="0"),out)
  out<-gsub(day_string,formatC(dd,width=2,flag="0"),out)
  out<-gsub(hour_string,formatC(hh,width=2,flag="0"),out)
  out<-gsub(min_string,formatC(MM,width=2,flag="0"),out)
  out<-gsub(sec_string,formatC(SS,width=2,flag="0"),out)
  out<-gsub(var_string,var.str,out)
  out
}

str2Rdate <- function(ts,format="%Y-%m-%d %H:%M:%S") {
# ===========================================================================
# converts a string into an R (POSIXt,POSIXct) date object
# date objects can be used with arithmetic operations +/-
# ts is a character or a vector of characters with date information
# in the format specified in format
# Output is a date object

     # the lengthy bunch of testing is necessary because strptime needs
     # explicit specifications for month and day, otherwise it returns NA.
     # this extension allows inputs of the format "%Y-%m" in which case the
     # the first day of the month is taken as a reference.

     #Êcheck if year/month/day is specified
     ysp <- length(c(grep("%Y",format,fixed=TRUE),
                     grep("%y",format,fixed=TRUE)))
     msp <- length(c(grep("%m",format,fixed=TRUE),
                     grep("%b",format,fixed=TRUE),
                     grep("%B",format,fixed=TRUE)))
     jsp <- length(c(grep("%j",format,fixed=TRUE)))
     dsp <- length(c(grep("%d",format,fixed=TRUE)))
     if (ysp > 1) { stop("ERROR: Multiple specification of year in 
                         date format.") }
     if (ysp == 0) { stop("ERROR: No year specification in 
                         date format.") }
     if (msp > 1) { stop("ERROR: Multiple specification of month in 
                         date format.") }
     if (dsp > 1) { stop("ERROR: Multiple specification of day in 
                         date format.") }

     # append month or day if not specified
     tss <- ts
     formati <- format
     if (jsp == 0) {
     if (msp == 0) { 
        tss <- paste(tss,"01",sep="")
        formati <- paste(formati,"%m",sep="")
     }
     if (dsp == 0) { 
        tss <- paste(tss,"01",sep="") 
        formati <- paste(formati,"%d",sep="")
     }
     }

     # this is necessary because strptime() returns NA otherwise
     as.POSIXct(strptime(tss,format=formati),tz="GMT")
}


#+ Box-Cox transformation
boxcox<-function(x,lambda) {
  if (lambda==0) {
    return(log(x))
  } else {
    return((x**lambda-1)/lambda)
  }
}

#+ Started Box-Cox transformation
started_boxcox <- function( x, lambda) {
  ic <- 1 
  if (lambda==0) {
    return(log(x))
  } else {
    res <- x
    ix  <- which(x > ic)
    res[ix] <- ic-(ic**lambda-1)/lambda + (x[ix]**lambda - 1) / lambda
    return(res)
  }
}

#+ Started Box-Cox backtransformation
inv_started_boxcox <- function( x, lambda, brrinf=-100) {
  ic <- 1 
  if (lambda==0) {
    res <- exp(x)
  } else {
    res <- x
    ix  <- which(x > ic)
    res[ix] <- (1 + lambda * (x[ix] - 1))**(1./lambda)
  }
  if (any(x < brrinf)) res[which(x<brrinf)]<-0
  res
}

##+ Correction function 
#correction_function_nobc <- function( ya, yb, eps2, Ya) {
#  Ya <- Ya + eps2 * (ya - yb)
#  Ya
#}

##+ Correction function for Box-Cox Inversion (works for lambda=0.5)
#correction_function_bc <- function( xia, xib, eps2, Xia, lambda, brrinf) {
#  Ya1 <- inv_started_boxcox(Xia, lambda, brrinf) + eps2/4 * (xia - xib)**2
#  Ya1
#}

#+ Correction function for Box-Cox Inversion (Evans)
inv_boxcox_wcorrection <- function( Xia_var, Xib_var, Xia, lambda, brrinf) {
  Ea <- tboxcox( Xia, lambda, brrinf)
  if ( nrow(Ea) == 1) {
    Ea_mu <- array( data = mean(Ea), dim = c(1,1))
  } else {
    Ea_mu <- rowMeans(Ea)
  }
  # if Xia_var == 0 then do not modify Ea
  if ( length( ix <- which(Xia_var != 0 & Xib_var == 0)) > 0) {
    Ea[ix,] <- Ea[ix,] + 0.5 * (1-lambda) * (Ea_mu[ix])**(1-2*lambda) * Xia_var[ix]
  } else if ( length( ix <- which(Xia_var != 0 & Xib_var != 0)) > 0) {
    Ea[ix,] <- Ea[ix,] + 0.5 * (1-lambda) * (Ea_mu[ix])**(1-2*lambda) * Xia_var[ix] * (1 - Xia_var[ix]/Xib_var[ix])
  }
  Ea
}

#+ Started Box-Cox Inversion with bias correction
inv_started_boxcox_wcorrection <- function( Xia, Xil, Xir, eps2, lambda, brrinf) {
  Ea <- inv_started_boxcox( Xia, lambda, brrinf)
  if ( nrow(Ea) == 1) {
    xia <- mean(Xia)
    xil <- mean(Xil)
    xir <- mean(Xir)
  } else {
    xia <- rowMeans(Xia)
    xil <- rowMeans(Xil)
    xir <- rowMeans(Xir)
  }
  xo <- inv_started_boxcox( xir + (1 + eps2) * (xia - xil), lambda, brrinf) 
  xl <- inv_started_boxcox( xil, lambda, brrinf)
  xr <- inv_started_boxcox( xir, lambda, brrinf)
  xa <- xl + 1/(1+eps2) * ( xo - xr)
  Ea <- Ea + xa - inv_started_boxcox( xia, lambda, brrinf)
  Ea
}

#+ Inverse Box-Cox transformation
tboxcox<-function(x,lambda,brrinf=-100) {
  if (lambda==0) {
    res<-exp(x)
  } else {
    res<-(1+lambda*x)**(1./lambda)
  }
  if (any(x<brrinf)) res[which(x<brrinf)]<-0
  res
}

#+
perturb_obs <- function( yo, var, setseed) {
  yop <- yo; yop[] <- NA
  if ( !is.na(setseed)) set.seed(setseed)
  if ( var == "RR1") {
    # The sd are taken from Figure 7 of Lussana et al (2023)
    if ( (nix1 <- length( ix1 <- which( !is.na(yo) & yo < 1))) > 0) {
      yop[ix1] <- yo[ix1] + rnorm( nix1, mean=0, sd=0.50) * yo[ix1]
    }
    if ( (nix2 <- length( ix2 <- which( !is.na(yo) & yo >= 1 & yo < 2))) > 0) {
      yop[ix2] <- yo[ix2] + rnorm( nix2, mean=0, sd=0.40) * yo[ix2]
    }
    if ( (nix3 <- length( ix3 <- which( !is.na(yo) & yo >= 2))) > 0) {
      yop[ix3] <- yo[ix3] + rnorm( nix3, mean=0, sd=0.30) * yo[ix3]
    }
    if (any(yop < 0)) yop[which(yop<0)] <- 0
  } else if ( var == "TA") {
    yop <- yo + rnorm( length(yo), mean=0, sd=1)
  }
  yop
}

#+ Define correlations 
corr1d <- function( values, par, label, values2=NULL, values3=NULL) {
# correlations when input vectors are referred to a specific point (e.g. distances from a point)
#
# par is the vector of parameters (1 for each vector of values to be processed)
# label is the vector of labels identifying the specific elaboration
# values is the vector/matrix of values to be processed (one vector is the column of an array)
#-----------------------------------------------------------------------------
  res <- rep( 1, length(values))

  for (i in 1:length(par)) {
    if (i==2) values <- values2
    if (i==3) values <- values3
    if (label[i] == "gaussian") {
      res <- res * exp( -0.5* (values*values) / (par[i]*par[i]) )

    } else if (label[i] == "soar")  {
      res <- res * (1+values/par[i])*exp(-values/par[i])

    } else if (label[i] == "powerlaw")  {
      res <- res * 1 / (1 + 0.5*(values*values)/(par[i]*par[i]))

    } else if (label[i] == "toar")  {
      res <- res * (1 + values/par[i] + (values*values)/(3*par[i]*par[i])) * exp(-values/par[i])

    } else if (label[i] == "linear")  {
      res <- res * (1 - (1-par[i]) * abs(values))

    } else {
      res <- rep( NA, length(values))

    } # end if

  } # end for

  res
}

#+ Define correlations
corr2d <- function( values, par, label, values_as_globalVar=F, 
                    values2=NULL, values3=NULL) {
# correlations between pair of points
#------------------------------------------------------------------------------
  for (i in 1:length(par)) {
    
    if (values_as_globalVar) {
      if (i==1) { 
        mat2 <- envtmp$dist2 
      } else if (i==2) { 
        mat2 <- envtmp$dist2_z 
      } else if (i==3) { 
        mat2 <- envtmp$dist2_laf 
      }
    } else if (i==1) {
      mat2 <- values 
    } else if (i==2) {
      mat2 <- values2 
    } else if (i==3) {
      mat2 <- values3
    } # end if 

    if (i==1) { res <- mat2; res[] <- 1 }

    if (label[i] == "gaussian") {
      res <- res * exp( -0.5* mat2 / (par[i]*par[i]) )

    } else if (label[i] == "soar")  {
      mat_norm <- sqrt(mat2) / par[i]
      res <- res * (1+mat_norm) * exp(-mat_norm)

    } else if (label[i] == "powerlaw")  {
      res <- res * 1 / (1 + 0.5 * mat2 / (par[i]*par[i]))

    } else if (label[i] == "toar")  {
      mat <- sqrt(mat2)
      res <- res * (1 + mat/par[i] + mat2/(3*par[i]*par[i])) * exp(-mat/par[i])

    } else if (label[i] == "linear")  {
      res <- res * (1 - (1-par[i]) * sqrt(mat2))

    } else {
      res <- rep( NA, length(values[,i]))

    } # end if

  } # end for
  res
}

#+ ensemble optimal interpolation time smoother
enoi <- function( i,
                  corr_dynamic  = "soar",
                  corr_static   = "soar",
                  corrz_dynamic = "gaussian",
                  corrz_static  = "gaussian",
                  corrlaf_dynamic = "gaussian",
                  corrlaf_static  = "gaussian",
                  alpha         = 0.5,
                  beta          = 0.5,
                  dh_static     = 10000,
                  dh_loc        = 10000,
                  dz_static     = 10000,
                  dz_loc        = 10000,
                  dlaf_static   = 0.5,
                  dlaf_loc      = 0.5,
                  k_dim_corr    = 10, 
                  m_dim         = 100000,
                  showdots      = F) {
# returned values: analysis, observation error var, analysis error var
# note: the call ‘t(x) %*% y’ (‘crossprod’) or ‘x %*% t(y)’ (‘tcrossprod’)
#------------------------------------------------------------------------------

  if( showdots & i%%(round(m_dim/10)) == 0) cat(".")

  if ( any( is.na(envtmp$El[i,]))) return( rep( NA, dim(envtmp$El)[2]))

  # select the observations to use
  if ( (p <- length( aux <- which(envtmp$nn2$nn.idx[i,]!=0))) == 0) {

    # no observations, analysis=background
    Ea <- envtmp$El[i,]
  } else {

    # available observations
    ixa  <- envtmp$nn2$nn.idx[i,aux]

    di <- array( data=envtmp$D[ixa,], dim=c(p,envtmp$k_dim))
    if (envtmp$k_dim > k_dim_corr) {
      ixb <- order( colMeans( abs(di)))[1:k_dim_corr]
    } else {
      ixb <- 1:envtmp$k_dim
      k_dim_corr <- envtmp$k_dim
    }

    # define vectors
    dist_xy <- envtmp$nn2$nn.dists[i,aux]
    x <- envtmp$obs_x[ixa]
    y <- envtmp$obs_y[ixa]
    z <- envtmp$obs_z[ixa]
    laf <- envtmp$obs_laf[ixa]
    zdist_xy <- envtmp$nn2$nn.zdists[i,aux]
    lafdist_xy <- envtmp$nn2$nn.lafdists[i,aux]
    Zr <- array( data=envtmp$Zr[ixa,ixb], dim=c(p,k_dim_corr))
    dist2_yy  <- outer(x,x,FUN="-")**2 + outer(y,y,FUN="-")**2
    zdist2_yy <- outer(z,z,FUN="-")**2
    lafdist2_yy <- outer(laf,laf,FUN="-")**2
    # combine static and dynamic correlations
    r_lr_xy <- array( data = 
              ( alpha * corr1d( values  = dist_xy, 
                                values2 = zdist_xy, 
                                values3 = lafdist_xy, 
                par   = c( dh_static, dz_static, dlaf_static), 
                label = c( corr_static, corrz_static, corrlaf_static)) + 
            (1-alpha) * corr1d( values  = dist_xy, 
                                values2 = zdist_xy, 
                                values3 = lafdist_xy, 
                par   = c( dh_loc, dz_loc, dlaf_static), 
                label = c( corr_dynamic, corrz_dynamic, corrlaf_dynamic)) * 
                    tcrossprod( envtmp$Xl[i,ixb], Zr)), dim=c(p,1))
    R_rr_yy <- array( data = (
                alpha * corr2d( values  = dist2_yy, 
                                values2 = zdist2_yy, 
                                values3 = lafdist2_yy, 
                par   = c( dh_static, dz_static, dlaf_static), 
                label = c( corr_static, corrz_static, corrlaf_static)) + 
            (1-alpha) * corr2d( values  = dist2_yy,
                                values2 = zdist2_yy,
                                values3 = lafdist2_yy, 
                par   = c( dh_loc, dz_loc, dlaf_static), 
                label = c( corr_dynamic, corrz_dynamic, corrlaf_dynamic)) * 
                    tcrossprod( Zr, Zr)), dim=c(p,p))

    # analysis
    X3 <- crossprod( chol2inv( chol( (R_rr_yy + diag(x=envtmp$eps2[i],p) ))), di)
    Ea <- envtmp$El[i,ixb] + envtmp$gamma[i] * beta * crossprod( r_lr_xy, X3)
  }
  return( Ea)
}

#+ IDI
idi <- function( i,
                 dh=10000,
                 dz=500,
                 dlaf=0.1, 
                 eps2=0.1, 
                 mode="idi" ) {
#------------------------------------------------------------------------------
  
  if ( is.na(envidi$tval[i])) return(NA)
  
  # select the observations to use
  if ( (p <- length( aux <- which(envidi$nn2$nn.idx[i,]!=0))) == 0) {

    # no observations, analysis=background
    if (mode == "idi") {
      res <- 1 / ( 1 + eps2)
    } else if (mode == "cvidi") {
      res <- 0 
    }
  } else {

    # available observations
    ixa  <- envidi$nn2$nn.idx[i,aux]

    # define vectors
    dist_xy <- envidi$nn2$nn.dists[i,aux]
    x <- envidi$x[ixa]
    y <- envidi$y[ixa]
    z <- envidi$z[ixa]
    laf <- envidi$laf[ixa]
    zdist_xy <- envidi$nn2$nn.zdists[i,aux]
    lafdist_xy <- envidi$nn2$nn.lafdists[i,aux]
    dist2_yy  <- outer(x,x,FUN="-")**2 + outer(y,y,FUN="-")**2
    zdist2_yy <- outer(z,z,FUN="-")**2
    lafdist2_yy <- outer(laf,laf,FUN="-")**2
    # combine static and dynamic correlations
    r_xy <- array( data = corr1d( values  = dist_xy, 
                                  values2 = zdist_xy, 
                                  values3 = lafdist_xy, 
                                  par     = c( dh, dz, dlaf), 
                                  label   = c( "gaussian", "gaussian", "linear")), dim=c(p,1))
    R_yy <- array( data = corr2d( values  = dist2_yy, 
                                  values2 = zdist2_yy, 
                                  values3 = lafdist2_yy, 
                                  par     = c( dh, dz, dlaf), 
                                  label   = c( "gaussian", "gaussian", "linear")), dim=c(p,p))

    SRinv <- chol2inv( chol( (R_yy + diag(x=eps2,p) ) ) )
    # idi
    res <- sum( r_xy * as.vector( rowSums( SRinv ) ) )
    if (mode == "cvidi") {
      ii   <- which( ixa == i)
      Wii  <- sum( r_xy * SRinv[ii,])
      # cv-idi
      res <- (res - Wii) / (1-Wii)
    }
  }
  return(res)
}

#+ Inflation
inflation <- function( Ens, fpos, fneg, idi, var) {
#------------------------------------------------------------------------------
  idi[idi > 1] <- 1
  idi[idi < 0] <- 0
  if ( nrow(Ens) == 1) {
    mu <- mean(Ens)
    mu_aux <- rep( mu, length(Ens))
  } else {
    mu <- rowMeans(Ens)
    mu_aux <- array( data=rep( mu, dim(Ens)[2]), dim=dim(Ens))
  }
  pert <- Ens - mu
  if ( length( ix <- which( pert > 0)) > 0) { 
    fpos_aux <- 1 + idi**2 / (idi**2 + (1-idi)**2) * (fpos - 1)
    fpos_aux[idi==0] <- 1
    fpos_mat <- array( data=rep( fpos_aux, dim(Ens)[2]), dim=dim(Ens))
    Ens[ix] <- mu_aux[ix] + fpos_mat[ix] * pert[ix]
    rm( fpos_aux, fpos_mat)
  }
  if ( length( ix <- which( pert < 0)) > 0) {
    fneg_aux <- 1 + idi**2 / (idi**2 + (1-idi)**2) * (fneg - 1)
    fneg_aux[idi==0] <- 1
    fneg_mat <- array( data=rep( fneg_aux, dim(Ens)[2]), dim=dim(Ens))
    Ens[ix] <- mu_aux[ix] + fneg_mat[ix] * pert[ix] 
    rm( fneg_aux, fneg_mat)
  }
  if (var == "RR1") Ens[Ens < 0] <- 0
  return(Ens)
}

###############################################################################
# Command line arguments
t0 <- Sys.time()
p <- arg_parser("twod")
p <- add_argument(p, "--config_file", help="configuration file", type="character", default=NA)
p <- add_argument(p, "--datel", help="date (%Y-%m-%dT%H)", type="character", default=NA)
p <- add_argument(p, "--varl", help="variable (RR1/TA)", type="character", default=NA)
p <- add_argument(p, "--dater", help="date (%Y-%m-%dT%H)", type="character", default=NA)
p <- add_argument(p, "--varr", help="variable (RR1/TA)", type="character", default=NA)
p <- add_argument(p, "--difftime", help="set dater as difference from datel", type="numeric", default=NA)
p <- add_argument(p, "--difftime_unit", help="time units", type="character", default=NA)
p <- add_argument(p, "--n_adjusts", help="total number of adjustements", type="integer", default=NA)
p <- add_argument(p, "--adjust", help="adjustment to run", type="integer", default=NA)
p <- add_argument(p, "--n_steps", help="total number of steps within one adjustment", type="integer", default=NA)
p <- add_argument(p, "--step", help="step to run", type="integer", default=NA)
p <- add_argument(p, "--cores", help="cores", type="integer", default=NA)
p <- add_argument(p, "--k_dim", help="ensembles", type="integer", default=30)
p <- add_argument(p, "--bvarl", help="variable (RR1/TA). dimensions=(number of varl)", type="character", nargs=Inf, default=NA)
p <- add_argument(p, "--bdater", help="date (%Y-%m-%dT%H)", type="character", nargs=Inf, default=NA)
p <- add_argument(p, "--bvarr", help="variable (RR1/TA). dimensions=(number of varl, number of steps)", type="character", nargs=Inf, default=NA)
p <- add_argument(p, "--bdifftime", help="set dater as difference from datel. dimensions=(number of varl, number of steps)", nargs=Inf, type="numeric", default=NA)
p <- add_argument(p, "--bdifftime_unit", help="time units. dimensions=(number of varl, number of steps)", type="character", default=NA)
p <- add_argument(p, "--ffin_bgl", help="LEFT-background input filename", type="character", default=NA)
p <- add_argument(p, "--ffin_bgl_template", help="LEFT-background input filename template", type="character", default=NA)
p <- add_argument(p, "--ffin_bgr", help="RIGHT-background input filename", type="character", default=NA)
p <- add_argument(p, "--ffin_bgr_template", help="RIGHT-background input filename template", type="character", default=NA)
p <- add_argument(p, "--ffin_l", help="left-updated ensemble input filename", type="character", default=NA)
p <- add_argument(p, "--ffin_l_template", help="left-updated ensemble input filename template", type="character", default=NA)
p <- add_argument(p, "--ffin_r", help="right-updated ensemble input filename", type="character", default=NA)
p <- add_argument(p, "--ffin_r_template", help="right-updated ensemble input filename template", type="character", default=NA)
p <- add_argument(p, "--ffin_obs", help="input observation filename", type="character", default=NA)
p <- add_argument(p, "--ffin_obs_template", help="input observation filename template", type="character", default=NA)
p <- add_argument(p, "--ffin_cvobs", help="input observation filename used for selecting cv-observations", type="character", default=NA)
p <- add_argument(p, "--ffin_cvobs_template", help="input observation filename used for selecting cv-observations", type="character", default=NA)
p <- add_argument(p, "--ffin_thinobs", help="input observation filename used for thinning", type="character", default=NA)
p <- add_argument(p, "--ffin_thinobs_template", help="input observation filename template used for thinning", type="character", default=NA)
p <- add_argument(p, "--ffout", help="output filename", type="character", default=NA)
p <- add_argument(p, "--ffout_template", help="output filename template", type="character", default=NA)
p <- add_argument(p, "--pmax", help="max number of nearby stations", type="numeric", default=NA)
p <- add_argument(p, "--bpmax", help="max number of nearby stations. dimensions=(number of adjustments, number of varl, number of steps)", type="numeric", nargs=Inf, default=NA)
p <- add_argument(p, "--corr_dynamic", help="dynamic correlation function (radial)", type="character", default=NA)
p <- add_argument(p, "--bcorr_dynamic", help="dynamic correlation function (radial). dimensions=(number of adjustments, number of varl, number of steps)", type="character",  nargs=Inf, default=NA)
p <- add_argument(p, "--corr_static", help="static correlation function (radial)", type="character", default=NA)
p <- add_argument(p, "--bcorr_static", help="static correlation function (radial). dimensions=(number of adjustments, number of varl, number of steps)", type="character",  nargs=Inf, default=NA)
p <- add_argument(p, "--corrz_dynamic", help="dynamic correlation function (vertical)", type="character", default=NA)
p <- add_argument(p, "--bcorrz_dynamic", help="dynamic correlation function (vertical). dimensions=(number of adjustments, number of varl, number of steps)", type="character",  nargs=Inf, default=NA)
p <- add_argument(p, "--corrz_static", help="static correlation function (vertical)", type="character", default=NA)
p <- add_argument(p, "--bcorrz_static", help="static correlation function (vertical). dimensions=(number of adjustments, number of varl, number of steps)", type="character", nargs=Inf, default=NA)
p <- add_argument(p, "--corrlaf_dynamic", help="dynamic correlation function (land area fraction)", type="character", default=NA)
p <- add_argument(p, "--bcorrlaf_dynamic", help="dynamic correlation function (land area fraction). dimensions=(number of adjustments, number of varl, number of steps)", type="character",  nargs=Inf, default=NA)
p <- add_argument(p, "--corrlaf_static", help="static correlation function (land area fraction)", type="character", default=NA)
p <- add_argument(p, "--bcorrlaf_static", help="static correlation function (land area fraction). dimensions=(number of adjustments, number of varl, number of steps)", type="character", nargs=Inf, default=NA)
p <- add_argument(p, "--eps2", help="obs-to-backg vars ratio", type="numeric", default=NA)
p <- add_argument(p, "--beps2", help="obs-to-backg vars ratio", type="numeric", nargs=Inf, default=NA)
p <- add_argument(p, "--gamma", help="left-to-right backg stdev ratio", type="numeric", default=NA)
p <- add_argument(p, "--bgamma", help="left-to-right backg stdev ratio", type="numeric", nargs=Inf, default=NA)
p <- add_argument(p, "--dh_static", help="static de-correlation length", type="numeric", default=NA)
p <- add_argument(p, "--bdh_static", help="static de-correlation length", type="numeric", nargs=Inf, default=NA)
p <- add_argument(p, "--dh_loc", help="localization de-correlation length", type="numeric", default=NA)
p <- add_argument(p, "--bdh_loc", help="localization de-correlation length", type="numeric", nargs=Inf, default=NA)
p <- add_argument(p, "--dh_idi", help="localization de-correlation length", type="numeric", default=NA)
p <- add_argument(p, "--bdh_idi", help="localization de-correlation length", type="numeric", nargs=Inf, default=NA)
p <- add_argument(p, "--dz_static", help="static de-correlation length", type="numeric", default=NA)
p <- add_argument(p, "--bdz_static", help="static de-correlation length", type="numeric", nargs=Inf, default=NA)
p <- add_argument(p, "--dz_loc", help="localization de-correlation length", type="numeric", default=NA)
p <- add_argument(p, "--bdz_loc", help="localization de-correlation length", type="numeric", nargs=Inf, default=NA)
p <- add_argument(p, "--dlaf_static", help="static de-correlation length", type="numeric", default=NA)
p <- add_argument(p, "--bdlaf_static", help="static de-correlation length", type="numeric", nargs=Inf, default=NA)
p <- add_argument(p, "--dlaf_loc", help="localization de-correlation length", type="numeric", default=NA)
p <- add_argument(p, "--bdlaf_loc", help="localization de-correlation length", type="numeric", nargs=Inf, default=NA)
p <- add_argument(p, "--dz_idi", help="localization de-correlation length", type="numeric", default=NA)
p <- add_argument(p, "--bdz_idi", help="localization de-correlation length", type="numeric", nargs=Inf, default=NA)
p <- add_argument(p, "--dlaf_idi", help="localization de-correlation length", type="numeric", default=NA)
p <- add_argument(p, "--bdlaf_idi", help="localization de-correlation length", type="numeric", nargs=Inf, default=NA)
p <- add_argument(p, "--alpha", help="tuning: correlation inflation", type="numeric", default=NA)
p <- add_argument(p, "--balpha", help="tuning: correlation inflation", type="numeric", nargs=Inf, default=NA)
p <- add_argument(p, "--beta", help="tuning: de-correlation length inflation based on the time difference", type="numeric", default=NA)
p <- add_argument(p, "--bbeta", help="tuning: de-correlation length inflation based on the time difference", type="numeric", nargs=Inf, default=NA)
p <- add_argument(p, "--pobs", help="Perturb observations", flag=T)
p <- add_argument(p, "--pobs_setseed", help="Seed random number generator", type="numeric", default=NA)
p <- add_argument(p, "--bpobs_setseed", help="Seed random number generator", type="numeric", nargs=Inf, default=NA)
p <- add_argument(p, "--cvobs", help="Cross-validation required", flag=T)
p <- add_argument(p, "--cvobs_setseed", help="Seed random number generator for Cross-validation", type="numeric", default=NA)
p <- add_argument(p, "--bcvobs_setseed", help="Seed random number generator for Cross-validation", type="numeric", nargs=Inf, default=NA)
p <- add_argument(p, "--cvobs_perc", help="Percentage of observations reserved fro Cross-Validation (0-100)", type="numeric", default=NA)
p <- add_argument(p, "--bcvobs_perc", help="Percentage of observations reserved fro Cross-Validation (0-100)", type="numeric", nargs=Inf, default=NA)
p <- add_argument(p, "--thinobs", help="Observation thinning required", flag=T)
p <- add_argument(p, "--thinobs_setseed", help="Seed random number generator for Thinning", type="numeric", default=NA)
p <- add_argument(p, "--bthinobs_setseed", help="Seed random number generator for Thinning", type="numeric", nargs=Inf, default=NA)
p <- add_argument(p, "--thinobs_perc", help="Percentage of observations reserved fro Thinning (0-100)", type="numeric", default=NA)
p <- add_argument(p, "--bthinobs_perc", help="Percentage of observations reserved fro Thinning (0-100)", type="numeric", nargs=Inf, default=NA)
p <- add_argument(p, "--skip", help="Skip the spatial interpolation but pass the information to the next step", flag=T)
p <- add_argument(p, "--lambda", help="Started Box-Cox transformation lambda", type="numeric", default=NA)
p <- add_argument(p, "--blambda", help="Started Box-Cox transformation lambda", type="numeric", nargs=Inf, default=NA)
p <- add_argument(p, "--transf_type", help="Data transformation (bc=Box-Cox; sbc=started Box-Cox)", type="character", default=NA)
p <- add_argument(p, "--btransf_type", help="Data transformation (bc=Box-Cox; sbc=started Box-Cox)", type="character", nargs=Inf, default=NA)
p <- add_argument(p, "--inflate", help="inflation", flag=T)
p <- add_argument(p, "--inflate_fpos", help="inflation factor when anomaly is positive", type="numeric", default=NA)
p <- add_argument(p, "--binflate_fpos", help="inflation factor when anomaly is positive", type="numeric", nargs=Inf, default=NA)
p <- add_argument(p, "--inflate_fneg", help="inflation factor when anomaly is negative", type="numeric", default=NA)
p <- add_argument(p, "--binflate_fneg", help="inflation factor when anomaly is negative", type="numeric", nargs=Inf, default=NA)
p <- add_argument(p, "--mstep", help="step to optimize memory usage during elaboration", type="integer", default=50000)

argv <- parse_args(p)

#-----------------------------------------------------------------------------
# read configuration file
if (!is.na(argv$config_file)) {
  if (file.exists(argv$config_file)) {
    source(argv$config_file)
    argv_tmp<-append(argv,conf)
    names_argv_tmp<-names(argv_tmp)
    argv_def<-list()
    names_argv_def<-integer(0)
    k<-0
    for (i in 1:length(argv_tmp)) {
      if (names_argv_tmp[i] %in% names_argv_def) next
      k<-k+1
      j<-which(names_argv_tmp==names_argv_tmp[i])
      argv_def[[k]]<-argv_tmp[[j[length(j)]]]
      names_argv_def<-c(names_argv_def,names_argv_tmp[i])
    }
    names(argv_def)<-names_argv_def
    rm(argv_tmp,names_argv_tmp,names_argv_def)
    rm(argv)
    argv<-argv_def
    rm(argv_def)
  } else {
    print("WARNING: config file not found")
    print(argv$config_file)
  }
}

#------------------------------------------------------------------------------
# Set parameters from big config files
if ( length( ixvarl <- which( argv$bvarl == argv$varl)) == 0 ) boom( "varl/bvarl", 2)
if ( argv$step > argv$n_steps) boom( "step", 2)
argv$varr <- argv$bvarr[argv$step,ixvarl]
argv$difftime <- argv$bdifftime[argv$step,ixvarl]
argv$difftime_unit <- argv$bdifftime_unit[argv$step,ixvarl]
argv$dater <- format( format = "%Y-%m-%dT%H", tz="UTC",
  seq( strptime( argv$datel, format="%Y-%m-%dT%H", tz="UTC"), 
       length=2, 
       by=paste(argv$difftime,argv$difftime_unit))[2])
if ( is.na(argv$ffin_bgl) ) {
  argv$ffin_bgl <- replaceDate( string   = argv$ffin_bgl_template,
                                date.str = argv$datel,
                                format   = "%Y-%m-%dT%H",
                                var.str  = argv$varl)
}
if ( is.na(argv$ffin_bgr) ) {
  argv$ffin_bgr <- replaceDate( string   = argv$ffin_bgr_template,
                                date.str = argv$dater,
                                format   = "%Y-%m-%dT%H",
                                var.str  = argv$varr)
}
if ( is.na(argv$ffin_l) ) {
  argv$ffin_l <- replaceDate( string   = argv$ffin_l_template,
                              date.str = argv$datel,
                              format   = "%Y-%m-%dT%H",
                              var.str  = argv$varl)
}
if ( is.na(argv$ffin_r) ) {
  argv$ffin_r <- replaceDate( string   = argv$ffin_r_template,
                              date.str = argv$dater,
                              format   = "%Y-%m-%dT%H",
                              var.str  = argv$varr)
}
if ( is.na(argv$ffin_obs) ) {
  if ( argv$varr == "RR1" ) obs_varr <- "rr1"
  if ( argv$varr == "TA"  ) obs_varr <- "ta"
  argv$ffin_obs <- replaceDate( string   = argv$ffin_obs_template,
                                date.str = argv$dater,
                                format   = "%Y-%m-%dT%H",
                                var.str  = obs_varr)
  if ( argv$varl == "RR1" ) obs_varl <- "rr1"
  if ( argv$varl == "TA"  ) obs_varl <- "ta"
  ffin_obs_l <- replaceDate( string   = argv$ffin_obs_template,
                             date.str = argv$datel,
                             format   = "%Y-%m-%dT%H",
                             var.str  = obs_varl)
}
if ( is.na(argv$ffout) ) {
  argv$ffout <- replaceDate( string   = argv$ffout_template,
                             date.str = argv$datel,
                             format   = "%Y-%m-%dT%H",
                             var.str  = argv$varl)
}
if (argv$cvobs & is.na(argv$ffin_cvobs)) {
  if ( argv$varl == "RR1" ) obs_varl <- "rr1"
  if ( argv$varl == "TA"  ) obs_varl <- "ta"
  argv$ffin_cvobs <- replaceDate( string   = argv$ffin_cvobs_template,
                                  date.str = argv$datel,
                                  format   = "%Y-%m-%dT%H",
                                  var.str  = obs_varl)
}
if (argv$thinobs & is.na(argv$ffin_thinobs)) {
  if ( argv$varl == "RR1" ) obs_varl <- "rr1"
  if ( argv$varl == "TA"  ) obs_varl <- "ta"
  argv$ffin_thinobs <- replaceDate( string   = argv$ffin_thinobs_template,
                                    date.str = argv$datel,
                                    format   = "%Y-%m-%dT%H",
                                    var.str  = obs_varl)
}
if ( is.na( argv$pmax))  
  argv$pmax <- argv$bpmax[argv$step,ixvarl,argv$adjust]
if ( is.na( argv$corr_dynamic))  
  argv$corr_dynamic <- argv$bcorr_dynamic[argv$step,ixvarl,argv$adjust]
if ( is.na( argv$corr_static))  
  argv$corr_static <- argv$bcorr_static[argv$step,ixvarl,argv$adjust]
if ( is.na( argv$corrz_dynamic))  
  argv$corrz_dynamic <- argv$bcorrz_dynamic[argv$step,ixvarl,argv$adjust]
if ( is.na( argv$corrz_static))  
  argv$corrz_static <- argv$bcorrz_static[argv$step,ixvarl,argv$adjust]
if ( is.na( argv$corrlaf_dynamic))  
  argv$corrlaf_dynamic <- argv$bcorrlaf_dynamic[argv$step,ixvarl,argv$adjust]
if ( is.na( argv$corrlaf_static))  
  argv$corrlaf_static <- argv$bcorrlaf_static[argv$step,ixvarl,argv$adjust]
if ( is.na( argv$eps2))  
  argv$eps2 <- argv$beps2[argv$step,ixvarl,argv$adjust]
if ( is.na( argv$gamma))  
  argv$gamma <- argv$bgamma[argv$step,ixvarl,argv$adjust]
if ( is.na( argv$dh_static))  
  argv$dh_static <- argv$bdh_static[argv$step,ixvarl,argv$adjust]
if ( is.na( argv$dh_loc))  
  argv$dh_loc <- argv$bdh_loc[argv$step,ixvarl,argv$adjust]
if ( is.na( argv$dh_idi))  
  argv$dh_idi <- argv$bdh_idi[argv$step,ixvarl,argv$adjust]
if ( is.na( argv$dz_static))  
  argv$dz_static <- argv$bdz_static[argv$step,ixvarl,argv$adjust]
if ( is.na( argv$dz_loc))  
  argv$dz_loc <- argv$bdz_loc[argv$step,ixvarl,argv$adjust]
if ( is.na( argv$dlaf_static))  
  argv$dlaf_static <- argv$bdlaf_static[argv$step,ixvarl,argv$adjust]
if ( is.na( argv$dlaf_loc))  
  argv$dlaf_loc <- argv$bdlaf_loc[argv$step,ixvarl,argv$adjust]
if ( is.na( argv$dlaf_idi))  
  argv$dlaf_idi <- argv$bdlaf_idi[argv$step,ixvarl,argv$adjust]
if ( is.na( argv$alpha))  
  argv$alpha <- argv$balpha[argv$step,ixvarl,argv$adjust]
if ( is.na( argv$beta))  
  argv$beta <- argv$bbeta[argv$step,ixvarl,argv$adjust]
if ( argv$pobs & is.na( argv$pobs_setseed))  
  argv$pobs_setseed <- argv$bpobs_setseed[argv$step,ixvarl,argv$adjust]
if ( argv$cvobs & is.na( argv$cvobs_setseed))  
  argv$cvobs_setseed <- argv$bcvobs_setseed[argv$step,ixvarl,argv$adjust]
if ( argv$cvobs & is.na( argv$cvobs_perc))  
  argv$cvobs_perc <- argv$bcvobs_perc[argv$step,ixvarl,argv$adjust]
if (is.na( argv$lambda) & any(!is.na(argv$blambda)))
  argv$lambda <- argv$blambda[argv$step,ixvarl,argv$adjust]
if (!is.na( argv$lambda) & is.na(argv$transf_type))
  argv$transf_type <- argv$btransf_type[argv$step,ixvarl,argv$adjust]
if ( argv$inflate & is.na( argv$inflate_fpos))  
  argv$inflate_fpos <- argv$binflate_fpos[argv$step,ixvarl,argv$adjust]
if ( argv$inflate & is.na( argv$inflate_fneg))  
  argv$inflate_fneg <- argv$binflate_fneg[argv$step,ixvarl,argv$adjust]

#------------------------------------------------------------------------------
# Background
#cat("Background:")
if ( !file.exists(argv$ffin_bgl)) boom( str=argv$ffin_bgl, code=1)
load(argv$ffin_bgl)
Ebl <- in2d_grid_val
if ( length( ixg <- grep(  "in2d", ls())) > 0 ) rm( list=ls()[ixg])
if ( length( ixg <- grep( "saved", ls())) > 0 ) rm( list=ls()[ixg])

if ( !file.exists(argv$ffin_bgr)) boom( str=argv$ffin_bgr, code=1)
load(argv$ffin_bgr)
Ebr <- in2d_grid_val
if ( length( ixg <- grep(  "in2d", ls())) > 0 ) rm( list=ls()[ixg])
if ( length( ixg <- grep( "saved", ls())) > 0 ) rm( list=ls()[ixg])

if ( !file.exists(argv$ffin_l)) boom( str=argv$ffin_l, code=1)
load(argv$ffin_l)
El <- in2d_grid_val
if ( length( ixg <- grep(  "in2d", ls())) > 0 ) rm( list=ls()[ixg])
if ( length( ixg <- grep( "saved", ls())) > 0 ) rm( list=ls()[ixg])

if ( !file.exists(argv$ffin_r)) boom( str=argv$ffin_r, code=1)
load(argv$ffin_r)
Er <- in2d_grid_val
#if ( length( ixg <- grep(  "in2d", ls())) > 0 ) rm( list=ls()[ixg])
if ( length( ixg <- grep( "saved", ls())) > 0 ) rm( list=ls()[ixg])

print( paste(" argv$ffin_obs:", argv$ffin_obs))
print( paste(" argv$ffin_bgl:", argv$ffin_bgl))
print( paste(" argv$ffin_bgr:", argv$ffin_bgr))
print( paste("   argv$ffin_l:", argv$ffin_l))
print( paste("   argv$ffin_r:", argv$ffin_r))
#cat( paste("size of El=", round(object.size(El)/(1024*1024)), "Mb\n"))
#cat( paste("size of Er=", round(object.size(Er)/(1024*1024)), "Mb\n"))

#------------------------------------------------------------------------------
# Observations

#cat("Observations:")
if ( ffin_obs_l != argv$ffin_obs) {
  load(ffin_obs_l)
  in2d_obs_val_l <- in2d_obs_val
  in2d_obs_y_l   <- in2d_obs_y
  in2d_obs_x_l   <- in2d_obs_x
  in2d_obs_z_l   <- in2d_obs_z
  in2d_obs_laf_l <- in2d_obs_laf
  if ( length( ixg <- grep( "saved", ls())) > 0 ) rm( list=ls()[ixg])
  rm( in2d_obs_val, in2d_obs_y, in2d_obs_x, in2d_obs_z, in2d_obs_laf)
}
if (!file.exists(argv$ffin_obs)) boom(str=argv$ffin_obs, code=1)
load(argv$ffin_obs)
npoints <- length(in2d_obs_val)
yo_orig <- array(data=NA, dim=c(npoints, argv$k_dim))
for (e in 1:argv$k_dim) yo_orig[,e] <- in2d_obs_val
# perturb observations (if needed)
# NOTE: if we use dynamical correlations, we should not perturb observations
#       use perturbed observations only with full static correlations
yo <- yo_orig
if (argv$pobs) yo <- perturb_obs( yo=yo, var=argv$varr, setseed=argv$pobs_setseed)
cat("\n")
#cat(paste("size of yo=", round(object.size(yo)/(1024*1024)), "Mb\n"))
k_dim <- dim(yo)[2]
if ( ffin_obs_l == argv$ffin_obs) {
  in2d_obs_val_l <- in2d_obs_val
  in2d_obs_y_l   <- in2d_obs_y
  in2d_obs_x_l   <- in2d_obs_x
  in2d_obs_z_l   <- in2d_obs_z
  in2d_obs_laf_l <- in2d_obs_laf
}
print( paste(" number of grid points:", dim(Er)[1]))
print( paste("   number of ensembles:", dim(Er)[2]))
print( paste("number of observations:", dim(yo)[1]))
#------------------------------------------------------------------------------
# Create environment
envtmp <- new.env( parent = emptyenv())
envtmp$y       <- in2d_grid_y
envtmp$m_dim   <- length(envtmp$y)
envtmp$x       <- in2d_grid_x
envtmp$z       <- getValues(in2d_grid_z)
envtmp$laf     <- getValues(in2d_grid_laf)
envtmp$k_dim   <- k_dim
envtmp$p_dim   <- npoints
envtmp$obs_y   <- in2d_obs_y
envtmp$obs_x   <- in2d_obs_x
envtmp$obs_z   <- in2d_obs_z
envtmp$obs_laf <- in2d_obs_laf
rmaster <- in2d_grid_z

#------------------------------------------------------------------------------
# Transform data (if needed)

if (!is.na(argv$lambda)) {
  if (argv$transf_type == "sbc" | argv$transf_type == "sbc0") {
    envtmp$yo <- started_boxcox( yo, lambda=argv$lambda)
    envtmp$El_tot <- started_boxcox( El, lambda=argv$lambda)
    envtmp$Er <- started_boxcox( Er, lambda=argv$lambda)
  } else if (argv$transf_type == "bc") {
    envtmp$yo <- boxcox( yo, lambda=argv$lambda)
    envtmp$El_tot <- boxcox( El, lambda=argv$lambda)
    envtmp$Er <- boxcox( Er, lambda=argv$lambda)
  }
} else {
  envtmp$yo <- yo
  envtmp$El_tot <- El
  envtmp$Er <- Er
}

#------------------------------------------------------------------------------
# Compute idi

if (argv$inflate) {
  envidi <- new.env( parent = emptyenv())
  # observation vectors
  envidi$p_dim <- length(obsop)
  envidi$y   <- in2d_obs_y
  envidi$x   <- in2d_obs_x
  envidi$z   <- in2d_obs_z
  envidi$laf <- in2d_obs_laf
  # target vectors
  envidi$tval <- envtmp$El_tot[,1]
  envidi$ty   <- envtmp$y
  envidi$tx   <- envtmp$x
  envidi$tz   <- envtmp$z
  envidi$tlaf<- envtmp$laf
  envidi$m_dim <- length(envidi$tx)
  ii <- 1
  jj <- argv$mstep
  res_idi <- vector( mode="numeric", length=envidi$m_dim)
  res_idi[] <- NA 
  if (ii == 1 & jj > envidi$m_dim) jj <-  envidi$m_dim
  while (jj <= envidi$m_dim) {
    m_dim_tmp <- jj - ii + 1
    envidi$nn2 <- nn2( cbind( envidi$x, envidi$y), 
                       query = cbind( envidi$tx[ii:jj], envidi$ty[ii:jj]), 
                       k = argv$pmax, searchtype = "radius", 
                       radius = (7*argv$dh_idi))
    idx <- envidi$nn2$nn.idx
    if (any(idx == 0)) idx[which(idx==0)] <- NA
    envidi$nn2$nn.zdists   <- array( data = abs(envidi$tz[ii:jj] - envidi$z[idx]), dim=dim(idx))
    envidi$nn2$nn.lafdists <- array( data = abs(envidi$tlaf[ii:jj] - envidi$laf[idx]), dim=dim(idx))
    rm(idx)
    if (!is.na(argv$cores)) {
      res_tmp <- t( mcmapply( idi,
                              1:m_dim_tmp,
                              mc.cores=argv$cores,
                              SIMPLIFY=T,
                              MoreArgs=list( dh   = argv$dh_idi,
                                             dz   = argv$dz_idi,
                                             dlaf = argv$dlaf_idi,
                                             mode = "idi")))
    # no-multicores
    } else {
      res_tmp <- t( mapply( idi,
                            1:m_dim_tmp,
                            SIMPLIFY=T,
                            MoreArgs=list( dh   = argv$dh_idi,
                                           dz   = argv$dz_idi,
                                           dlaf = argv$dlaf_idi,
                                           mode = "idi")))
    }
    res_idi[ii:jj,] <- res_tmp
    ii <- ii + argv$mstep
    jj <- jj + argv$mstep
    if ( ii <= envidi$m_dim & jj > envidi$m_dim) jj <- envidi$m_dim
  }
  grid_idi  <- envtmp$El_tot[,1]
  grid_idi[] <- NA
  grid_idi[] <- res_idi
  rm(envidi, res_idi, res_tmp, obsop)
} else {
  grid_idi <- numeric(0) 
}
#------------------------------------------------------------------------------
# Skip this observation time (used for testing)
if (argv$skip) {
  in2d_obs_val  <- in2d_obs_val
  in2d_grid_val <- envtmp$El_tot
  in2d_grid_idi <- grid_idi
  argv_saved <- argv
  yo_saved <- yo
  envtmp_saved <- envtmp
  El_saved <- envtmp$El_tot 
  Er_saved <- Er
  Ea_saved <- envtmp$El_tot
  Ebr_saved <- Ebr
  Ebl_saved <- Ebl
  save(file=argv$ffout, envtmp_saved, El_saved, Er_saved, Ebl_saved, Ebr_saved, Ea_saved, yo_saved, argv_saved, in2d_obs_val, in2d_grid_val, in2d_grid_y, in2d_grid_x, in2d_grid_z, in2d_grid_idi, in2d_grid_laf, in2d_obs_y, in2d_obs_x, in2d_obs_z, in2d_obs_laf)
  q()
}

#------------------------------------------------------------------------------
# -- EnSI --
envtmp$HEr <- envtmp$yo; envtmp$HEr[] <- NA
HEbr <- envtmp$HEr
for (e in 1:envtmp$k_dim) {
  rmaster[] <- envtmp$Er[,e]
  envtmp$HEr[,e] <- extract(rmaster, cbind(envtmp$obs_x,envtmp$obs_y))
  rmaster[] <- Ebr[,e]
  HEbr[,e] <- extract(rmaster, cbind(envtmp$obs_x,envtmp$obs_y))
}
envtmp$eps2  <- rep(argv$eps2, envtmp$m_dim)
envtmp$D    <- array(data=NA,dim=c(envtmp$p_dim,envtmp$k_dim))
envtmp$D <- envtmp$yo - envtmp$HEr
Ebl_mean <- rowMeans(Ebl)
Ebl_sd   <- apply(Ebl, FUN=function(x){sd(x)}, MAR=1)
HEbr_mean <- rowMeans(HEbr)
HEbr_sd <- apply( HEbr, FUN=function(x){sd(x)}, MAR=1)
Ebr_sd   <- apply( Ebr, FUN=function(x){sd(x)}, MAR=1)
if (argv$gamma <= 0) argv$gamma <- mean(Ebl_sd,na.rm=T) / mean(Ebr_sd,na.rm=T)
#print(mean(Ebl_sd,na.rm=T))
#print(mean(Ebr_sd,na.rm=T))
#print(mean(Ebl_sd,na.rm=T)/mean(Ebr_sd,na.rm=T))
envtmp$Zr    <- array( data=NA, dim=c(envtmp$p_dim,envtmp$k_dim))
envtmp$Zr[,] <- 1/sqrt(envtmp$k_dim-1) * (HEbr - HEbr_mean) / HEbr_sd
envtmp$Zr[!is.finite(envtmp$Zr)] <- 0
rm(HEbr, HEbr_mean, HEbr_sd)
ii <- 1
jj <- argv$mstep
res <- array( data=NA, dim=c(envtmp$m_dim, envtmp$k_dim))
if ( jj > envtmp$m_dim ) jj <-  envtmp$m_dim
while (jj <= envtmp$m_dim) {
  m_dim_tmp <- jj - ii + 1
  envtmp$eps2  <- rep(argv$eps2, m_dim_tmp) 
  envtmp$gamma <- rep(argv$gamma, m_dim_tmp)
  envtmp$El <- envtmp$El_tot[ii:jj,]
  envtmp$Xl    <- array( data=NA, dim=c(m_dim_tmp,envtmp$k_dim))
  envtmp$Xl[,] <- 1/sqrt(envtmp$k_dim-1) * (Ebl[ii:jj,] - Ebl_mean[ii:jj]) / Ebl_sd[ii:jj]
  envtmp$Xl[!is.finite(envtmp$Xl)] <- 0 
  envtmp$nn2 <- nn2( cbind( envtmp$obs_x, envtmp$obs_y), 
                     query = cbind( envtmp$x[ii:jj], envtmp$y[ii:jj]), 
                     k = argv$pmax, searchtype = "radius", 
                     radius = (7*argv$dh_loc))
  idx <- envtmp$nn2$nn.idx
  if (any(idx == 0)) idx[which(idx==0)] <- NA
  envtmp$nn2$nn.zdists   <- array( data = abs(  envtmp$z[ii:jj] - envtmp$obs_z[idx]), dim=dim(idx))
  envtmp$nn2$nn.lafdists <- array( data = abs(envtmp$laf[ii:jj] - envtmp$obs_laf[idx]), dim=dim(idx))
  rm(idx)
  if (!is.na(argv$cores)) {
    res_tmp <- t( mcmapply( enoi,
                            1:m_dim_tmp,
                            mc.cores=argv$cores,
                            SIMPLIFY=T,
                            MoreArgs=list( corr_dynamic  = argv$corr_dynamic,
                                           corr_static   = argv$corr_static,
                                           corrz_dynamic = argv$corrz_dynamic,
                                           corrz_static  = argv$corrz_static,
                                           corrlaf_dynamic = argv$corrlaf_dynamic,
                                           corrlaf_static  = argv$corrlaf_static,
                                           alpha         = argv$alpha,
                                           beta          = argv$beta,
                                           dh_static     = argv$dh_static,
                                           dh_loc        = argv$dh_loc,
                                           dz_static     = argv$dz_static,
                                           dz_loc        = argv$dz_loc,
                                           dlaf_static   = argv$dlaf_static,
                                           dlaf_loc      = argv$dlaf_loc,
                                           k_dim_corr    = envtmp$k_dim,
                                           showdots = F)))
  # no-multicores
  } else {
    res_tmp <- t( mapply( enoi,
                          1:m_dim_tmp,
                          SIMPLIFY=T,
                          MoreArgs=list( corr_dynamic  = argv$corr_dynamic,
                                         corr_static   = argv$corr_static,
                                         corrz_dynamic = argv$corrz_dynamic,
                                         corrz_static  = argv$corrz_static,
                                         corrlaf_dynamic = argv$corrlaf_dynamic,
                                         corrlaf_static  = argv$corrlaf_static,
                                         alpha         = argv$alpha,
                                         beta          = argv$beta,
                                         dh_static     = argv$dh_static,
                                         dh_loc        = argv$dh_loc,
                                         dz_static     = argv$dz_static,
                                         dz_loc        = argv$dz_loc,
                                         dlaf_static   = argv$dlaf_static,
                                         dlaf_loc      = argv$dlaf_loc,
                                         k_dim_corr    = envtmp$k_dim,
                                         showdots = F)))
  }
  res[ii:jj,] <- res_tmp
  ii <- ii + argv$mstep
  jj <- jj + argv$mstep
  if ( ii <= envtmp$m_dim & jj > envtmp$m_dim) jj <- envtmp$m_dim
}
#------------------------------------------------------------------------------
# Back-transform data (if needed)

if (!is.na(argv$lambda)) {
  if (argv$transf_type == "sbc") {
    if (argv$varl == "RR1") { min <- 0 } else { min <- (-9999)}
    Ea <- inv_started_boxcox_wcorrection( Xia=res, Xil=envtmp$El_tot, Xir=envtmp$Er, eps2=argv$eps2, lambda=argv$lambda, brrinf=min) 
  } else if (argv$transf_type == "sbc0") {
    if (argv$varl == "RR1") { min <- 0 } else { min <- (-9999)}
    Ea <- inv_started_boxcox( res, argv$lambda, min) 
  } else if (argv$transf_type == "bc") {
    if (argv$varl == "RR1") { min <- 0 } else { min <- (-9999)}
    avar <- apply( res, FUN=function(x){var(x)}, MAR=1)
    bvar <- apply( envtmp$El_tot, FUN=function(x){var(x)}, MAR=1)
    Ea <- inv_boxcox_wcorrection( Xia=res, Xia_var=avar, Xib_var=bvar , lambda=argv$lambda, brrinf=min) 
  }
} else {
  Ea <- res
}
if (argv$varl == "RR1") Ea[Ea<0] <- 0
#------------------------------------------------------------------------------
# Inflation
if (argv$inflate) 
  Ea <- inflation( Ens  = Ea,
                   fpos = argv$inflate_fpos, 
                   fneg = argv$inflate_fneg, 
                   idi  = grid_idi, 
                   var  = argv$varl)

#------------------------------------------------------------------------------
# Write output

npoints_l <- length(in2d_obs_val_l)
ya <- array(data=NA, dim=c(npoints_l, argv$k_dim))
yb <- array(data=NA, dim=c(npoints_l, argv$k_dim))
for (e in 1:envtmp$k_dim) {
  rmaster[] <- Ea[,e]
  ya[,e] <- extract(rmaster, cbind(in2d_obs_x_l,in2d_obs_y_l))
  rmaster[] <- envtmp$El_tot[,e]
  yb[,e] <- extract(rmaster, cbind(in2d_obs_x_l,in2d_obs_y_l))
}
rmse_b <- sqrt( mean( (yb[,1] - in2d_obs_val_l)**2 ))
rmse_a <- sqrt( mean( (ya[,1] - in2d_obs_val_l)**2 ))
mae_b <- mean( abs( yb[,1] - in2d_obs_val_l ))
mae_a <- mean( abs( ya[,1] - in2d_obs_val_l ))
cat("\n")
in2d_obs_val  <- in2d_obs_val
in2d_grid_val <- Ea
in2d_grid_idi <- grid_idi
argv_saved <- argv
yo_saved <- yo
envtmp_saved <- envtmp
El_saved <- envtmp$El_tot 
Er_saved <- Er
Ea_saved <- Ea
Ebr_saved <- Ebr
Ebl_saved <- Ebl
print( paste("          gamma:", round( argv$gamma, 4) ))
print( paste("   corr_dynamic:", argv$corr_dynamic))
print( paste("    corr_static:", argv$corr_static))
print( paste("  corrz_dynamic:", argv$corrz_dynamic))
print( paste("   corrz_static:", argv$corrz_static))
print( paste("corrlaf_dynamic:", argv$corrlaf_dynamic))
print( paste(" corrlaf_static:", argv$corrlaf_static))
print( paste("          alpha:", round( argv$alpha,4) ))
print( paste("           beta:", round( argv$beta,4) ))
print( paste("      dh_static:", round( argv$dh_static,4) ))
print( paste("         dh_loc:", round( argv$dh_loc,4) ))
print( paste("      dz_static:", round( argv$dz_static,4) ))
print( paste("         dz_loc:", round( argv$dz_loc,4) ))
print( paste("    dlaf_static:", round( argv$dlaf_static,4) ))
print( paste("       dlaf_loc:", round( argv$dlaf_loc,4) ))
print( paste("          k_dim:", argv$k_dim))
print("")
print( paste( "yo:", round( min( in2d_obs_val_l, na.rm=T), 2),
                     round( max( in2d_obs_val_l, na.rm=T), 2) ))
print( paste( "yb:", round( min( yb[,1], na.rm=T), 2),
                     round( max( yb[,1], na.rm=T), 2) ))
print( paste( "ya:", round( min( ya[,1], na.rm=T), 2),
                     round( max( ya[,1], na.rm=T), 2) ))
print( paste( "El:", round( min( El_saved[,1], na.rm=T), 2),
                     round( max( El_saved[,1], na.rm=T), 2) ))
print( paste( "Ea:", round( min( Ea_saved[,1], na.rm=T), 2),
                     round( max( Ea_saved[,1], na.rm=T), 2) ))
print( paste( " mae_b:", round( mae_b, 2) ))
print( paste( " mae_a:", round( mae_a, 2) ))
print( paste( "rmse_b:", round( rmse_b, 2) ))
print( paste( "rmse_a:", round( rmse_a, 2) ))
save(file=argv$ffout, envtmp_saved, El_saved, Er_saved, Ebl_saved, Ebr_saved, Ea_saved, yo_saved, argv_saved, in2d_obs_val, in2d_grid_val, in2d_grid_y, in2d_grid_x, in2d_grid_z, in2d_grid_idi, in2d_grid_laf, in2d_obs_y, in2d_obs_x, in2d_obs_z, in2d_obs_laf)
cat(paste("written file:",argv$ffout,"\n"))
t1 <- Sys.time()
cat( paste( "total time=", round(t1-t0,1), attr(t1-t0,"unit"),"\n"))

#------------------------------------------------------------------------------
# Normal exit
q()
#load("/home/cristianl/data/ensi-metno/output/ensi-metno_2d/case_20230807/RR1_k01/analysis_RR1_adj01_20230807T03Z_RR1_20230807T05Z_step5.rdata")
# rmaster[]<-Ea_saved[,1]
# breaks<-c(0,seq(0.1,28,length=9));col<-c("gray",rev(rainbow(8))); image(rmaster1,breaks=breaks,col=col)
# rmaster-rmaster1
