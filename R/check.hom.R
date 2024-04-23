# Function for testing the homogeneity of variance
# Date: 24/11/2023. Works only with lm objects
check.hom <- function(obj, var = NULL, alternative = c("varIdent",
                                                       "varPower",
                                                       "varExp")){
  # Exclude objects != lm
   if(!inherits(obj, "lm"))
     stop("This method works only with 'lm' object")

  # sets the alternative
  alternative <- match.arg(alternative)

  # Retrieve the dataset from model call
  if(!is.null(eval(obj$call$data)))
    dframe <- eval(obj$call$data)

  # var: uses data column or external vector
  tmp1 <- try(is.vector(var), silent = T)
  if(!is(tmp1, "try-error")){
    ## print(var)
    # Other code to execute
    if(is.null(var)) var <- fitted(obj)
    } else {
    var <- deparse(substitute(var))
    tmp2 <- try(dframe[, var], silent = T)
    if(is(tmp2, "try-error")) {
      stop("Classification variable not found in data or in workspace")
      } else {
        # print("Data")
        var <- tmp2
      }
    }

  # finding the type of variable
  # print(var)
  # print(inherits(var, "numeric"))
  if(!inherits(var, "numeric") & alternative != "varIdent")
    stop("varPower and varExp can only be used with numeric covariates or no covariates")

  # Fitting gls and comparing
  dframe2 <- obj$model
  dframe2$claVar <- factor(var)
  mod1 <- gls(formula(obj$call$formula),
              data = dframe2)
  if(alternative == "varIdent"){
    mod2 <- gls(formula(obj$call$formula),
              data = dframe2,
              weights = varIdent(form = ~ 1|claVar))
    } else if(alternative == "varPower"){
    # print(var)
    dframe2$claVar <- var
    mod2 <- gls(formula(obj$call$formula),
              data = dframe2,
              weights = varPower())
    } else if(alternative == "varExp"){
    dframe2$claVar <- var
    mod2 <- gls(formula(obj$call$formula),
              data = dframe2,
              weights = varExp())
    }
  LRT <- abs(as.numeric(2 * (logLik(mod1) - logLik(mod2))))
  df1 <- attr(logLik(mod1), "df")
  df2 <- attr(logLik(mod2), "df")
  dfLRT <- abs(df1 - df2)
  plev <- pchisq(as.numeric(LRT), dfLRT, lower.tail = F)
  LRTobj <- anova(mod1, mod2)
  returnList <- list("LRT" = LRT, "P-level" = plev, "aovtable" = LRTobj,
                     modHet = mod2)
  invisible(returnList)
}
