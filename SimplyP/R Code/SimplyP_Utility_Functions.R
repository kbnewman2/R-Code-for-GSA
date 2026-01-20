
# File name: 0_SimplyP_Utility_Functions.R
# 30 May 2020
# 10 August 2020: refreshing memory, "clearly" local R code
# 25 May 2021: removed some now in KenUtilities (mytitle);
#    added SimplyP_Projections
# 15 Feb 2024: moved to new Folder for LSM 1.2 work;
#  Also made SimplyP_Run_General SimplyP_Run
#

# 2. CV
# 3. gamma.pars: given mean and CV find shape and rate
# 4. IQR
# 5. exceedance.frac
# 6. zero.one.scale:  transform x on [a,b] to x* on [0,1]
# 7. rev.zero.one: backtransform u from [0,1] to x on [a,b]
# 8. SimplyP_Run: make multiple runs w/ SimplyP and return output
# 9. exp.cov:   exponential covariance function ---
# 10. rel.absPE
# 11. SimplyP_projections
# 12. multiv.chull.sample

cat("brought in   CV, gamma.pars, IQR, exceedance.frac,
    zero.one.scale, rev.zero.one, SimplyP_Run, exp.cov, rel.abs.PE,
    SimplyP_projections \n")


# 2 ----------------------------------------------------
CV      <- function(x) return(sd(x)/mean(x))

# 3 ----------------------------------------------------
gamma.pars <- function(ybar,yCV) {
  # given mean and CV find shape and rate
  a <- 1/yCV^2
  b <- a/ybar
  out <- list(alpha=a,beta=b)
  return(out)
}

# 4 ----------------------------------------------------
IQR <- function(x) {
  out <- quantile(x,probs=0.75)-quantile(x,probs=0.25)
  return(out)
}

# 5 ----------------------------------------------------
exceedance.frac <- function(x,exceed=0.03) {
  # Return fraction of days with x > exceed
  ok <- !is.na(x)
  x <- x[ok]
  n <- length(x)
  out <- sum(x>exceed)/n
  return(out)
}

# 6 ----------------------------------------------------
# transform x on [a,b] to x* on [0,1]
zero.one.scale <- function(x,a,b) (x-a)/(b-a)

# 7 ----------------------------------------------------
rev.zero.one   <- function(u,a,b) u*(b-a)+a


# 8 ----------------------------------------------------
SimplyP_Run <- function(num.param.comb,
                        data.path,
                        # the following 4 are now passed values
                        param.file,
                        input.file,
                        new.start.date,
                        new.end.date,

                        # Parameter information
                        num.par,
                        full.par.names,
                        par.class,
                        par.sample,

                        # Output information
                        num.output.vars,
                        var.names,
                        var.indices,
                        output.names.short=NULL,
                        verbose=FALSE,

                        # Precipitation change
                        original.precip=NULL,
                        original.airtemp=NULL,
                        input_data_start_date=NULL,
                        precip.multipliers=NULL
) {

  #num.param.comb = number of parameter combinations, thus runs w/ SimplyP
  #num.par = number of parameters that will have changed values
  #full.par.names = names of those parameters
  #par.class = the class of those parameters
  #par.sample = matrix of param values w/ dim = num.param.comb by num.par

  #num.output.vars = number of output variables that will be extracted
  #var.names = names of those output variables that SimplyP recognizes
  #var.indices = indicates a catchment (e.g., Coull), or land use type (e.g., arable)
  #output.names.short = customized names for the output variables

  # checks
  if(num.par != length(full.par.names)) stop("Inconsistent no. of params")
  if(num.output.vars != length(var.names)) stop("inconsistent no. of outputs")

  # Initialize
  mobius_setup_from_parameter_and_input_file(
    ParameterFileName=param.file, InputFileName=input.file)

  # The following overrides the start date and number of days
  none <- vector(mode='character', length=0)
  date.seq      <- seq(ymd(new.start.date), ymd(new.end.date),by=1)
  num.timesteps <- length(date.seq)
  mobius_set_parameter_time(Name     = 'Start date',
                            IndexesIn= none,
                            Value    = new.start.date)
  #mobius_set_parameter_uint('End date', none, num.timesteps)

  # Create a 3-D array for output from SimplyP
  output.array <- array(data=NA,
          dim=c(num.output.vars,num.param.comb,num.timesteps),
          dimnames=list(output.names.short,NULL,NULL))

  do.run <- TRUE
  if(do.run) {
  # Loop over the sets of input parameter combinations
  for(i in 1:num.param.comb) {

   if(verbose) { if((i %% 100) == 0) cat("Combination",i,"\n") }

    # If the precipitation gets an across-the-board multiplication
    if(!is.null(precip.multipliers)) {
      new.precip <- original.precip*precip.multipliers[i]
      mobius_setup_from_parameter_file_and_input_series(
        ParameterFileName = param.file,    #TarlandParameters_v0-3.dat
        InputDataStartDate= input_data_start_date,
        AirTemperature    = original.airtemp,
        Precipitation     = new.precip)
      # Need to override the model run start date and number of days to simulate
      mobius_set_parameter_time(Name     = 'Start date',
                                IndexesIn= none,
                                Value    = new.start.date)
      mobius_set_parameter_time(Name     = 'End date',
                                IndexesIn= none,
                                Value    = new.end.date)
      # mobius_set_parameter_uint('End date', none, num.timesteps)
    }

    # Set the parameter values
    for(j in 1:num.par) {
      if(par.class[j] != "") {  #Not a "none" class
        mobius_set_parameter_double(Name=full.par.names[j],
                                    IndexesIn=par.class[j],
                                    Value=par.sample[i,j])
      } else {
        mobius_set_parameter_double(Name=full.par.names[j],
                                    IndexesIn= none,
                                    Value = par.sample[i,j])
      }
    }

    # ----------------- Run SimplyP --------------------------------------
    mobius_run_model()

    # Extract the Outputs
    for(j in 1:num.output.vars) {
      temp <- mobius_get_result_series(
        Name=var.names[j], IndexesIn=var.indices[[j]])
      #output.array[j,i,] <- temp
      if(j==1) cat("output length=",length(temp),"\n")
    }

  }  # End of Loop over the sets of input parameter combinations
  } # end of do.run block

  # re-initialize: this will set to 485 days again
  mobius_setup_from_parameter_and_input_file(
    ParameterFileName=param.file, InputFileName=input.file)

  output <- list(output.array=output.array)
  return(output)

}

# 9 ----------------------------------------------------
exp.cov <- function(theta1,theta2,sigma2,omega.vec) {
  out <- sigma2*exp(- sum((theta1-theta2)^2/omega.vec^2))
  return(out)
}

# 10  ----------------------------------------------------
rel.absPE <- function(est,true) {
  #---- relative MSPE----
  out <- mean(abs((est-true)/true))
  return(out)
}

# 11  ----------------------------------------------------
SimplyP_projections <- function(eta,y,t1,t2,var.name=NULL,plot.it=TRUE,
                                my.cex=1,
                                verbose=FALSE) {
  #eta is a n x T1 matrix of SimplyP output on one time series
  # with column names being character representations of Dates
  #y is a length T2 time series of observations, where obs'ns
  #  have names being character representations of Dates
  # 7 Jan 2022: modifying so that obs'ns, y, need not be sequential

  num.sims <- nrow(eta)
  if(is.character(t1)) {
    s <- ymd(t1)
    e <- ymd(t2)
    date.seq <- seq(s,e,by=1)
    date.seq.char <- as.character(date.seq)
  }
  if(is.Date(t1)) {
    date.seq <- seq(t1,t2,by=1)
    date.seq.char <- as.character(date.seq)
  }
  eta.sub <- eta[,date.seq.char]

  #ok obs'ns
  temp <- as.Date(names(y))
  ok   <- temp >= as.Date(date.seq[1]) & temp <= as.Date(tail(date.seq,n=1))
  y.sub <- y[ok]
  y.dates <- as.Date(names(y.sub))

  #y.sub   <- y[date.seq.char]
  if(verbose) {
    cat("#obs'ns=",length(y.sub),"\n")
    print(y.sub)
  }
  if(plot.it) {
    my.ylim <- range(c(eta.sub,y.sub),na.rm=TRUE)

    # plot(date.seq,y.sub,type="l",xlab="Date",ylab="",
    #      ylim=my.ylim,col="red",main=var.name)
    plot(date.seq,eta.sub[1,],type="n",xlab="Date",ylab="",
         ylim=my.ylim,main=var.name)
    for(i in 1:num.sims) {
      lines(date.seq,eta.sub[i,])
    }
    #lines(y.dates,y.sub,col="red",lwd=2,lty=2)
    points(y.dates,y.sub,col="red",lwd=2,lty=2,pch=3,cex=my.cex)
    # lines(date.seq,y.sub,col="red",lwd=2,lty=2)
  } else {
    out <- list(eta.sub=eta.sub,y.sub=y.sub)
    return(out)
  }

}


# 12  ----------------------------------------------------
# ************************* Function for next Wave Sampling *********************
multiv.chull.sample <- function(X,N=2,inflater=1000,verbose=FALSE) {
  #takes an n by p matrix X, n=#obs'ns, p=#vars
  #and generates a N by p matrix of points in the
  #p dimensional convex hull (...if that makes sense)

  n <- nrow(X)
  p <- ncol(X)
  var.names <- colnames(X)

  num.pairs <- p*(p-1)/2
  if(verbose) cat("number of pairs=",num.pairs,"\n")

  #Find the lower and upper bounds per variable
  L.var <- apply(X,2,min)
  U.var <- apply(X,2,max)

  #Create a list of all the convex hulls
  convex.hull.list <- list()
  tally <- 0
  for(i in 1:(p-1)) {
    for(j in (i+1):p) {
      tally <- tally+1
      convex.hull.list[[tally]] <- convhulln(p=X[,c(i,j)],output.options = TRUE)
    }
  }

  oversample <- p*N*inflater
  undone <- TRUE
  results <- matrix(NA,nrow=1,ncol=p)
  number.iter <- 1

  while(undone) {
    LHS <-  randomLHS(n=oversample,k=p)
    rescaled <- LHS
    dimnames(rescaled)[[2]] <- var.names

    # LHS returns values on (0,1), the following rescales to allowed ranges
    for(i in 1:p) {
      var.range <- U.var[i] - L.var[i]
      rescaled[,i] <- LHS[,i]*var.range+L.var[i]
    }

    #find which values are in the joint convex hull
    ok.set <- matrix(data=FALSE,nrow=oversample,ncol=num.pairs)
    tally <- 0
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        tally <- tally+1
        ch.case <- convex.hull.list[[tally]]
        in_ch <- inhulln(ch.case, rescaled[,c(i,j)])
        ok.set[,tally] <- in_ch
      }
    }

    all.ok <- apply(ok.set,1,all)
    # print(ok.set)
    # cat("\n")
    if(verbose) {
      cat("number ok=",sum(all.ok),"out of",oversample,
          "or",round(100*sum(all.ok)/oversample),"% \n")
    }

    kept.combos <- rescaled[all.ok,]
    results <- rbind(results,kept.combos)
    if(nrow(results >=N)) {
      undone <- FALSE
    } else{
      number.iter <- number.iter+1
      if(verbose) {
        cat("another iteration=",number.iter,"\n")
      }
    }

  }
  out <- list(convex.hull.list=convex.hull.list,num.pairs=num.pairs,
              total.sample=rescaled,output=results[-1,])
  return(out)
}

