#.RcoclustEnv <- new.env(parent = emptyenv());
#.onLoad <- function(libname, pkgname){}

.noGenerics <- TRUE

.onUnload <- function(libpath)
  library.dynam.unload("Rcoclust", libpath)

# -------------------------------------------------------------

.check_data <- function(envrdata) {
  if (exists("envrdata"))
    if (!is.null("envrdata"))
      if (!is.null("envrdata$n") && !is.null("envrdata$d") && 
          !is.null("envrdata$Ais_byrow") && !is.null("envrdata$Ajs_byrow") && 
          !is.null("envrdata$Ajs_bycol") && !is.null("envrdata$Axs_bycol") && 
          !is.null("envrdata$Axs_byrow") && !is.null("envrdata$Ais_bycol") &&
          !is.null("envrdata$nnzi") && !is.null("envrdata$nnzj") &&
          !is.null("envrdata$nnz") &&  !is.null("envrdata$name"))
        if (length(envrdata$Ais_byrow)==envrdata$nnz && length(envrdata$Ajs_byrow)==envrdata$nnz && 
            length(envrdata$Ajs_bycol)==envrdata$nnz && length(envrdata$Axs_bycol)==envrdata$nnz && 
            length(envrdata$Axs_byrow)==envrdata$nnz && length(envrdata$Ais_bycol)==envrdata$nnz)
          if (max(envrdata$Ais_byrow)<envrdata$n && min(envrdata$Ais_byrow)>=0 &&
              max(envrdata$Ais_bycol)<envrdata$n && min(envrdata$Ais_bycol)>=0 &&
              max(envrdata$Ajs_byrow)<envrdata$d && min(envrdata$Ajs_byrow)>=0 &&
              max(envrdata$Ajs_bycol)<envrdata$d && min(envrdata$Ajs_bycol)>=0 &&
              min(envrdata$Axs_byrow)>=0 && min(envrdata$Axs_bycol)>=0)
              if ((sum(envrdata$nnzi)==envrdata$nnz) && (sum(envrdata$nnzj)==envrdata$nnz))
               {return (1);} else { return (0); }
}

# -------------------------------------------------------------

randp <- function(envrdata,dimr,sgr,vect_Ar,transfrm,debug) {
  if (.check_data(envrdata) && dimr>=1 && sgr>0 &&  
      transfrm%in%0:3 && debug%in%0:1 && !is.null(vect_Ar) && 
      length(as.numeric(vect_Ar))==envrdata$n*dimr) { 
    ans1 <- .Call("RPMAT", 
        as.integer(envrdata$Ais_byrow), as.integer(envrdata$Ajs_byrow),
        as.numeric(envrdata$Axs_byrow), as.integer(envrdata$nnzi), as.integer(envrdata$nnzj), 
        as.numeric(vect_Ar),
        as.numeric(c(envrdata$n,envrdata$d,envrdata$nnz,dimr,sgr,0,transfrm,debug)));
    return (1);
  }
  return (0);
}

# -------------------------------------------------------------

cc_coclus <- function(g,envrdata,zi_coclus,wj_coclus,transfrm,maxiter,debug) {
  if (.check_data(envrdata) && g%in%2:envrdata$n &&
      transfrm%in%0:3 && maxiter%in%1:300 && debug%in%0:1 &&   
      !is.null(zi_coclus) && length(as.numeric(zi_coclus))==envrdata$n &&
      !is.null(wj_coclus) &&  length(as.numeric(wj_coclus))==envrdata$d) {  
    ans1 <- .Call("COCLUS",
                  as.integer(envrdata$Ais_byrow), as.integer(envrdata$Ajs_byrow),
                  as.numeric(envrdata$Axs_byrow), as.integer(envrdata$nnzi), 
                  as.integer(envrdata$Ais_bycol), as.integer(envrdata$Ajs_bycol),
                  as.numeric(envrdata$Axs_bycol), as.integer(envrdata$nnzj), 
                  as.integer(zi_coclus), as.integer(wj_coclus), 
                  as.numeric(c(envrdata$n, envrdata$d, g, envrdata$nnz, transfrm, maxiter, debug,1E-7)));
    
    return (list(obj=ans1, zi=zi_coclus, wj=wj_coclus));
  }
  return (NULL);
}

# -------------------------------------------------------------

cc_fddkm <- function(g,envrdata,zi_fddkm,wj_fddkm,vect_cik_fddkm,vect_wjk_fddkm,alpha,beta,delta,transfrm,maxiter,debug) {
  if (.check_data(envrdata) && g%in%2:envrdata$n && alpha>1 && beta>1 && 
      (delta==-1 || (delta>=min(envrdata$Axs_byrow) && delta<=max(envrdata$Axs_byrow))) &&
      transfrm%in%0:3 && maxiter%in%1:300 && debug%in%0:1 &&   
      !is.null(zi_fddkm) && length(as.numeric(zi_fddkm))==envrdata$n &&
      !is.null(wj_fddkm) &&  length(as.numeric(wj_fddkm))==envrdata$d &&
      !is.null(vect_cik_fddkm) &&  length(as.numeric(vect_cik_fddkm))==envrdata$n*g &&
      !is.null(vect_wjk_fddkm) &&  length(as.numeric(vect_wjk_fddkm))==envrdata$d*g ) {
      ans1 <- .Call("FDDKM",
                 as.integer(envrdata$Ais_byrow), as.integer(envrdata$Ajs_byrow), 
                 as.numeric(envrdata$Axs_byrow), as.integer(envrdata$nnzi), 
                 as.integer(envrdata$Ais_bycol), as.integer(envrdata$Ajs_bycol), 
                 as.numeric(envrdata$Axs_bycol), as.integer(envrdata$nnzj), 
                 as.numeric(vect_cik_fddkm), as.numeric(vect_wjk_fddkm),
                 as.integer(zi_fddkm), as.integer(wj_fddkm),
                 as.numeric(c(envrdata$n, envrdata$d, g, envrdata$nnz, alpha, beta, delta, transfrm, maxiter, debug,1E-7)));
      
      return (list(obj=ans1, zi=zi_fddkm, wj=wj_fddkm));
  }
  return (NULL);
}

# -------------------------------------------------------------

cc_ddkm <- function(g,envrdata,zi_ddkm,wj_ddkm,delta,transfrm,maxiter,debug) {
  if (.check_data(envrdata) && g%in%2:envrdata$n &&
      (delta==-1 || (delta>=min(envrdata$Axs_byrow) && delta<=max(envrdata$Axs_byrow))) &&
      transfrm%in%0:3 && maxiter%in%1:300 && debug%in%0:1 &&   
      !is.null(zi_ddkm) && length(as.numeric(zi_ddkm))==envrdata$n &&
      !is.null(wj_ddkm) &&  length(as.numeric(wj_ddkm))==envrdata$d) {
    ans1 <- .Call("DDKM",
                  as.integer(envrdata$Ais_byrow), as.integer(envrdata$Ajs_byrow),
                  as.numeric(envrdata$Axs_byrow), as.integer(envrdata$nnzi), 
                  as.integer(envrdata$Ais_bycol), as.integer(envrdata$Ajs_bycol),
                  as.numeric(envrdata$Axs_bycol), as.integer(envrdata$nnzj), 
                  as.integer(zi_ddkm), as.integer(wj_ddkm),
                  as.numeric(c(envrdata$n, envrdata$d, g, envrdata$nnz, delta, transfrm, maxiter, debug,1E-7)));
    
    return (list(obj=ans1, zi=zi_ddkm, wj=wj_ddkm));
  }
  return (NULL);
}

# -------------------------------------------------------------
