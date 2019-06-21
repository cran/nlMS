model_selection<-function(data=NULL, name_dir=NULL, name_results=NULL){ #data=X|Y

  if(is.null(data)) {
    print("No data available. Stop running")
    return()}
  if(is.null(name_dir))
    name_dir=tempdir()
  if(is.null(name_results))
    name_results="nlMS_results"

  o<-order(data[,1])
  data<-data[o,]
  ## include var names
  var_names<-colnames(data)[-1] #exclude first column, supposed to keep X variable
  X<-as.numeric(data[,1])

  n<-length(var_names) #number of variables

  index<-c("AIC", "RSE", "R2", "residuals")
  N<-length(X) #sample size

  models<-c("linear", "quadratic", "cubic", "logistic", "exponential", "power",
            "monod", "haldane", "logit")
  formulas<-c(Y~a+b*X, #linear
              Y~a+b*X+c*X2, #quadratic
              Y~a+b*X+c*X2+d*X3, #cubic
              Y~a/(1+exp(-(X-X0)*b)), #logistic
              Y~a*exp(b*X), #exponential growth
              Y~a*X**b+c, #power
              Y~a*X/(b+X), #monod
              Y~a*X/(b+c*X+X2), #haldane
              Y~X0+1/b*log((X+1)/(a-(X+1)))) #logit

  formulas<-lapply(formulas, stats::as.formula)

  #file name with global results of each model
  file_name_global<-paste0(name_dir, "/", name_results, ".csv")
  #writes first line in output file
  cat("variable", "index", models, "best fit", "\n", file=file_name_global,
      sep=";", append=FALSE)
  table<-matrix(numeric(length(models)*length(index)), nrow=length(index))
  row.names(table)<-c("AIC", "RSE", "R2", "res.norm")


  for(i in 1:n) {#for each variable
    #print(c("variable ", var_names[i]))
    j<-which(colnames(data)==var_names[i]) #column in file corresponding to var i
    Y<-as.numeric(data[,j]) #define variable Y
    X2<-X**2
    X3<-X**3
    starts<-list(list(a=summary(stats::lm(Y~X))$coefficients[1], b=summary(stats::lm(Y~X))$coefficients[2]), #linear
                 list(a=summary(stats::lm(Y~X+X2))$coefficients[1], b=summary(stats::lm(Y~X+X2))$coefficients[2], c=summary(stats::lm(Y~X+X2))$coefficients[3]), #quadratic
                 list(a=summary(stats::lm(Y~X+X2+X3))$coefficients[1], b=summary(stats::lm(Y~X+X2+X3))$coefficients[2], c=summary(stats::lm(Y~X+X2+X3))$coefficients[3], d=summary(stats::lm(Y~X+X2+X3))$coefficients[4]), #cubic
                 list(a=mean(Y[which(X==max(X))]), b=1/stats::sd(X), X0=mean(X)), #logistic
                 list(a=max(0.001, abs(min(Y))), b=log(max(abs(Y[which(X==max(X))]))/max(abs(min(Y)), 0.001))/max(X)), #exponential growth
                 list(a=max(0.0001, exp(log(max(Y))-log(max(Y)/max(0.0001, min(Y)))*log(X[N]))), b=log(max(Y)/max(0.0001, min(Y))), c=max(0.0001, min(Y))), #power
                 list(a=sign(Y[which.max(abs(Y))])*max(abs(Y)), b=mean(X)), #monod
                 list(a=sign(Y[which.max(abs(Y))])*mean(abs(Y[which(X==X[1])])), b=mean(X[1])**2, c=mean(Y)), #haldane
                 list(a=X[N]+2, b=sign(mean(Y[which(X==X[N])])-mean(Y[which(X==X[1])]))*1/stats::sd(Y)), X0=mean(Y)) #logit

    if(is.na(stats::var(Y)) | stats::var(Y)<=1e-40) cat(var_names[i], "var(Y)=0", "\n", file=file_name_global, sep=";", append=TRUE)
    else{
      SST<-sum((Y-mean(Y))**2) #total sum of squares
      table=table*0

      file_name<-paste0(name_dir, "/", var_names[i], ".txt") #file name for the results of each variable
      plot_name<-paste0(name_dir, "/", var_names[i], ".png") #file name for the figures of the models

      #from now, all plots are saved in this file
      grDevices::png(plot_name, width=1000, height = 1000, res=100) #it can also be pdf, jpeg, tiff or bmp
      graphics::par(mfrow=c(3,3),mar=c(3,4,2,0), las=1, cex=1, cex.lab=1, font.lab=1, cex.axis=0.8, family="serif")

      #writes first line in output file
      cat(" variable ", var_names[i], "\n\n", file = file_name)
      for(m in 1:length(models)){
        mod=NULL
        mes<-utils::capture.output(mod<-tryCatch(nlme::gnls(formulas[[m]], start=starts[[m]]), error=function(e) NULL ))
        if(!is.null(mod)) {
          k=length(starts[[m]])
          out<-func(mod, X, Y, models[m], k, SST, N, file_name, formulas[[m]])
          table[1,m]<-out$AIC_val
          table[2,m]<-out$RSE_val
          table[3,m]<-out$R2_val
          table[4,m]<-out$norm
        }
      }
      grDevices::dev.off()

      #write results to file
      for(line in 1:4){
        if(line<3){#
          min_mod=min(subset(table[line,], table[line,]!=0))
          min_mod=models[which(table[line,]==min_mod)]
        }
        #detect nans and infinities
        if(length(which(is.nan(table)))>0)
          print(paste(which(is.nan(table)), "is nan; for ", name_results, "and variable", var_names[i], "\n"))
        if(length(which(is.infinite(table)))>0)
          print(paste(which(is.infinite(table)), " is infinite; for ", name_results, " and variable ", var_names[i], "\n"))
        #write indices
        if(line!=4) cat(var_names[i], index[line], table[line,], "\t", file=file_name_global, sep=";", append = TRUE)
        if(line<3) cat(min_mod, file=file_name_global, sep=";", append = TRUE)
        if(line==4) cat(var_names[i], index[line], replace(replace(x <- table[line,], x==1, "norm"), x==0, "no norm"), "\t", file=file_name_global, sep=";", append = TRUE)
        cat("\n", file=file_name_global, sep=";", append = TRUE)
      }
      cat("\n", file=file_name_global, sep=";", append = TRUE)
    }#end var(Y)>0
  }#end variable
}#end function

func<-function(mod, X, Y, name_mod, k=2, SST, N, file_name=NULL, form=NULL, make_plot=TRUE){
  a=NULL; b=NULL; c=NULL; d=NULL; X0=NULL;
  fitted=mod$fitted
  residuals=Y-fitted

  #fitness index
  AIC_val<-stats::AIC(mod)
  SSE<-sum(residuals**2, na.rm=TRUE) #sum of squared errors
  RSE_val<-sqrt(SSE/(N-k)) # root (mean) squared errors
  R2_val=0

  #writing result in file
  if(!is.null(file_name)){
    cat("\n \n model: ", name_mod," \n", file=file_name, sep='', append=TRUE)
    cmd <- utils::tail(as.character(form),1)
    exp0 <- parse(text=cmd)
    cat(" equation: ", cmd," \n", file=file_name, sep='', append=TRUE)

    a=mod$coefficients[1]
    b=mod$coefficients[2]
    cat(" parameters: a = ", a, "; b = ", b, file=file_name, sep='', append=TRUE)

    c=mod$coefficients[3]
    d=mod$coefficients[4]
    if(name_mod=="haldane"){
      a=a/c
      b=b/c
    }
    if(name_mod=="logistic") X0=mod$coefficients[3]
    if(!is.null(c)) if(!is.na(c)) cat(" ; c = ", c, file=file_name, sep='', append=TRUE)
    if(!is.null(d)) if(!is.na(d)) cat(" ; d = ", d, file=file_name, sep='', append=TRUE)
    if(name_mod=="logistic" || name_mod=="logit") X0=c
    if(!is.null(X0)) if(!is.na(X0)) cat(" ; X0 = ", X0, file=file_name, sep='', append=TRUE)

    if(make_plot==TRUE) {
      graphics::plot(X, Y, pch=20, main=name_mod, ylab="", xlab="", xlim=c(min(0, min(X)), max(0,max(X))), ylim=c(min(0, min(Y)), max(0,max(Y))))
      params<-list(a=a, b=b)
      if(!is.null(c))
        params<-list(a=a, b=b, c=c)
      if(!is.null(d))
        params<-list(a=a, b=b, c=c, d=d)
      if(!is.null(X0))
        params<-list(a=a, b=b, X0=X0)
      z<-seq(min(X), max(X), (max(X)-min(X))/100)
      params$X<-z
      if(name_mod=="quadratic" || name_mod=="haldane") params$X2<-z**2
      if(name_mod=="cubic") {params$X2<-z**2; params$X3<-z**3}
      graphics::lines(z, eval(exp0, params))#f(list(X=z, params)))
    }
    cat("\n AIC = ", AIC_val, "; RSE = ", RSE_val, file=file_name, sep='', append=TRUE)
    if(name_mod=="linear" || name_mod=="quadratic" || name_mod=="cubic")
      cat("; R^2 = ", 1-SSE/SST, file=file_name, sep='', append=TRUE) #R2 value
    a="no-test"
    p=-1
    if(!is.na(stats::var(residuals)) & stats::var(residuals)>0){
      p=stats::shapiro.test(residuals)$p.value
      if(p>=0.05) {a="norm";b=1}
      if(p<0.05) {a="no-norm";b=0}
    }
    cat("\n residuals distr: ", a, " p-value: ", p, " \n", file=file_name, sep='', append=TRUE)
    list(AIC_val=AIC_val, RSE_val=RSE_val, R2_val=R2_val, norm=b)
  }
}

