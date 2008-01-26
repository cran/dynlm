dynlm <- function(formula, data, subset, weights, na.action,
		method = "qr", model = TRUE, x = FALSE, y = FALSE,
		qr = TRUE, singular.ok = TRUE, contrasts = NULL,
		offset, start = NULL, end = NULL, ...)
{

    ## setup environment in which some auxiliary functions are defined
    Zenv <- new.env(parent = environment(formula)) ##FIXME## better parent?
    ## turn formula into dynformula
    assign("dynformula", function(x) structure(x, class = unique(c("dynformula", oldClass(x)))), envir = Zenv)
    ## define convenience versions of diff() and lag() and season()
    assign("L", function(x, k = 1) {      
      if(length(k) > 1) {
        rval <- lapply(k, function(i) lag(x, k = -i))
	rval <- if(inherits(x, "ts")) do.call("ts.intersect", rval)
	  else do.call("merge", c(rval, list(all = FALSE)))
        colnames(rval) <- k
      } else {
        rval <- lag(x, k = -k)
      }
      return(rval)      
    }, envir = Zenv)
    assign("d", function(x, lag = 1) diff(x, lag = lag), envir = Zenv)
    assign("season", function(x, ref = NULL) 
    {
        freq <- frequency(x)
	stopifnot(freq > 1 && identical(all.equal(freq, round(freq)), TRUE))
	freq <- ofreq <- round(freq)
	freq <- if(freq == 12) month.abb else
	        if(freq == 4) paste("Q", 1:4, sep = "") else
		1:freq	
        rval <- factor(coredata(cycle(x)), labels = freq)
	if(!is.null(ref)) rval <- relevel(rval, ref = ref)
        rval <- zoo(rval, index(x), ofreq)
	return(rval)
    }, envir = Zenv)
    ## model.frame.dynformula relies on merge.zoo
    assign("model.frame.dynformula", function (formula, data = NULL, subset = NULL, 
    na.action = na.omit, drop.unused.levels = FALSE, xlev = NULL, ...) 
    {
	if (is.null(data)) data <- parent.frame()
	if (!is.list(data)) data <- as.list(data)
	args <- as.list(attr(terms(formula), "variables"))[-1]
	args$retclass <- "list"
	args$all <- FALSE
	formula <- terms(formula)
	attr(formula, "predvars") <- as.call(append(merge.zoo, args))
	attr(formula, "predvars")[[1]] <- as.name("merge.zoo")
	NextMethod("model.frame", formula = formula)
    }, envir = Zenv)

    ## original class of the dependent variable
    if(missing(data)) data <- Zenv
    orig.class <- if(is.data.frame(data) || is.environment(data))
                    class(eval(attr(terms(formula), "variables")[[2]], data, Zenv))
		  else class(data)

    ## prepare model.frame call
    ret.x <- x
    ret.y <- y
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    
    ## formula processing for 2SLS (if necessary)    
    if(length(formula[[3]]) > 1 && identical(formula[[3]][[1]], as.name("|"))) {
      twostage <- TRUE
      ff <- formula
      mf$formula[[3]][1] <- call("+")      
      ff1 <- . ~ .
      ff2 <- ~ .
      ff1[[2]] <- ff[[2]]
      ff1[[3]] <- ff[[3]][[2]]
      ff2[[3]] <- ff[[3]][[3]]
      ff2[[2]] <- NULL
    } else {
      twostage <- FALSE
    }
    
    ## call model.frame
    mf[[1]] <- as.name("model.frame")
    mf[[2]] <- as.call(list(as.name("dynformula"), mf[[2]]))
    ## and evaluate in Zenv such that L() and d() are known
    mf <- eval(mf, envir = Zenv)

    ## preserve zoo index properties and zoo factors
    mfna <- attr(mf, "na.action")
    if(length(index(mf[,1])) > nrow(mf)) {
      for(i in 1:NCOL(mf)) attr(mf[,i], "index") <- attr(mf[,i], "index")[-as.vector(mfna)]
    }
    is.zoofactor <- function(x) !is.null(attr(x, "oclass")) && attr(x, "oclass") == "factor"
    for(i in 1:NCOL(mf)) if(is.zoofactor(mf[,i])) mf[,i] <- coredata(mf[,i])
    
    ## process start and end
    mf1 <- mf[,1] ## reference series
    start <- if(is.null(start)) 1
               else {
	         if(length(start) > 1) start <- start[1] + (start[2] - 1)/frequency(mf1)
                 start <- min(which(index(mf1) >= start))
               }
    end <- if(is.null(end)) length(mf1)
             else {
	       if(length(end) > 1) end <- end[1] + (end[2] - 1)/frequency(mf1)
               end <- max(which(index(mf1) <= end))
             }
    if(end < start) {
      warning("empty model frame specified")
      mf <- head(mf, 0)
      mf1 <- head(mf1, 0)
    } else {
      mf <- mf[start:end,,drop=FALSE]
      mf1 <- mf1[start:end]
      if(!is.null(mfna))
        attr(mf, "na.action") <- mfna[as.vector(mfna) >= start & as.vector(mfna) <= end]
    }
        
    ## convert back to "ts" or "numeric"
    if("ts" %in% orig.class && is.regular(mf1, strict = TRUE)) {
      for(i in 1:ncol(mf)) if(!is.factor(mf[,i])) mf[,i] <- zoo:::as.ts.zoo(mf[,i])
    }
    if(all(orig.class == "numeric")) {
      for(i in 1:ncol(mf)) if(!is.factor(mf[,i])) mf[,i] <- as.vector(mf[,i])
    }

    ## assign rownames
    rownames(mf) <- index2char(index(mf1), frequency(mf1))

    if (method == "model.frame")
	return(mf)
    else if (method != "qr")
	warning("method = ", method, " is not supported. Using \"qr\".")
    mt <- attr(mf, "terms")         # allow model.frame to update it
    attr(mt, "predvars") <- NULL    #FIXME#
    attr(mt, "dataClasses") <- NULL #FIXME#
    y <- model.response(mf, "numeric")
    w <- model.weights(mf)
    offset <- model.offset(mf)
    if(!is.null(offset) && length(offset) != NROW(y)) stop("Number of offsets is ",
      length(offset), ", should equal ", NROW(y), " (number of observations)")

    if (is.empty.model(mt)) {
      x <- NULL
      rval <- list(coefficients = numeric(0), residuals = y, fitted.values = 0 * y,
	weights = w, rank = 0, df.residual = length(y))
      if(!is.null(offset)) rval$fitted.values <- offset
    }
    else {
      ## auxiliary regression
      if(twostage) {
        mt <- terms(ff1)
	y1 <- y
        x1 <- model.matrix(mt, mf)
        z <- model.matrix(terms(ff2), mf)
	xnam <- colnames(x1)
	ynam <- colnames(y1)
	auxreg <- if(is.null(w)) lm.fit(z, cbind(y1, x1), offset = offset, singular.ok=singular.ok, ...)
          else lm.wfit(z, cbind(y1, x1), w, offset = offset, singular.ok=singular.ok, ...)
        x <- auxreg$fitted[, -(1:NCOL(y1)), drop = FALSE]
        y <- auxreg$fitted[, 1:NCOL(y1)]
	colnames(x) <- xnam
	colnames(y) <- ynam
      } else {
        x <- model.matrix(mt, mf, contrasts)      
      }
      
      ## main regression
      rval <- if(is.null(w)) lm.fit(x, y, offset = offset, singular.ok=singular.ok, ...)
	else lm.wfit(x, y, w, offset = offset, singular.ok=singular.ok, ...)
    }
    class(rval) <- c(if(is.matrix(y)) "mlm", "lm")
    rval$na.action <- attr(mf, "na.action")
    rval$offset <- offset
    rval$contrasts <- attr(x, "contrasts")
    rval$xlevels <- .getXlevels(mt, mf)
    rval$call <- cl
    rval$terms <- mt
    if(model) rval$model <- mf
    if(ret.x) rval$x <- x
    if(ret.y) rval$y <- y
    rval$index <- index(mf1)
    rval$frequency <- frequency(mf1)
    rval$twostage <- twostage
    if(twostage) {
      rval$formula <- ff
      rval$residuals <- y1 - predict(rval, newdata = as.data.frame(x1))
      rval$fitted <- y1 - rval$residuals
    }

    class(rval) <- c("dynlm", class(rval))
    return(rval)
}

index.dynlm <- time.dynlm <- function(x, ...) {
  x$index
}

start.dynlm <- function(x, ...) {
  start(x$residuals)
}

end.dynlm <- function(x, ...) {
  end(x$residuals)
}

print.dynlm <- function(x, ...) {
  rx <- residuals(x)
  cat(paste("\nTime series regression with \"", class(rx)[1], "\" data:\n", sep = ""))
  cat(paste("Start = ", index2char(index(rx)[1], x$frequency),
            ", End = ",   index2char(index(rx)[length(rx)], x$frequency), "\n", sep = ""))
  NextMethod()
}

summary.dynlm <- function(object, vcov. = NULL, df = NULL, ...) {
  rval <- NextMethod()
  rval$frequency <- object$frequency
  if(any(c(!is.null(vcov.), !is.null(df)))) {
    coefmat <- coeftest(object, vcov. = vcov., df = df)
    attr(coefmat, "method") <- NULL
    class(coefmat) <- "matrix"
    rval$coefficients <- coefmat
    Rmat <- if(attr(object$terms, "intercept"))
      cbind(0, diag(length(coef(object))-1)) else diag(length(coef(object)))
    rval$fstatistic[1] <- linear.hypothesis(object, Rmat, vcov. = vcov.)[2,"F"]
    ## FIXME: make this default for linear.hypothesis.default?
    ## FIXME: seem to need to set r = 0 for hypothesis printing?
    if(!is.null(df)) rval$fstatistic[3] <- df
  }
  class(rval) <- c("summary.dynlm", class(rval))
  return(rval)
}

print.summary.dynlm <- function(x, ...) {
  rx <- residuals(x)
  x$residuals <- coredata(x$residuals)
  cat(paste("\nTime series regression with \"", class(rx)[1], "\" data:\n", sep = ""))
  cat(paste("Start = ", index2char(index(rx)[1], x$frequency),
            ", End = ",   index2char(index(rx)[length(rx)], x$frequency), "\n", sep = ""))
  NextMethod()
}

recresid.dynlm <- function(x, ...) {
  stopifnot(require("strucchange"))
  rval <- NextMethod()
  res <- residuals(x)
  if(inherits(res, "zoo")) {
    rval <- zoo(rval, time(res)[-(1:length(coef(x)))], frequency = attr(res, "frequency"))
  } else {
    rval <- ts(rval, end = end(res), frequency = frequency(res))
  }
  return(rval)
}
