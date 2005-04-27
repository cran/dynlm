dynlm <- function(formula, data, subset, weights, na.action,
		method = "qr", model = TRUE, x = FALSE, y = FALSE,
		qr = TRUE, singular.ok = TRUE, contrasts = NULL,
		offset, start = NULL, end = NULL, ...)
{
    ## setup environment in which some auxiliary functions are defined
    Zenv <- new.env(parent = environment(formula))
    ## turn formula into dynformula
    assign("dynformula", function(x) structure(x, class = unique(c("dynformula", oldClass(x)))), envir = Zenv)
    ## define convenience versions of diff() and lag() and season()
    assign("L", function(x, k = 1) lag(x, k = -k), envir = Zenv)
    assign("d", function(x, k = 1) diff(x, k = k), envir = Zenv)
    assign("season", function(x) 
    {
        freq <- frequency(x)
	stopifnot(freq > 1 && identical(all.equal(freq, round(freq)), TRUE))
	freq <- ofreq <- round(freq)
	freq <- if(freq == 12) month.abb else
	        if(freq == 4) paste("Q", 1:4, sep = "") else
		1:freq	
        zoo(factor(coredata(cycle(x)), labels = freq), index(x), ofreq)
    }, envir = Zenv)
    ## model.frame.dynformula relies on merge.zoo
    assign("model.frame.dynformula", function (formula, data = NULL, subset = NULL, 
    na.action = na.omit, drop.unused.levels = FALSE, xlev = NULL, ...) 
    {
	if (is.null(data)) data <- parent.frame()
	if (!is.list(data)) data <- as.list(data)
	args <- as.list(attr(terms(formula), "variables"))[-1]
	args$retclass <- "data.frame"
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
    
    ## call model.frame
    mf[[1]] <- as.name("model.frame")
    mf[[2]] <- as.call(list(as.name("dynformula"), mf[[2]]))
    ## and evaluate in Zenv such that l() and d() are known
    mf <- eval(mf, envir = Zenv)

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
    }    
        
    ## convert back to "ts" or "numeric"
    if("ts" %in% orig.class && is.regular(mf1, strict = TRUE)) {
      for(i in 1:ncol(mf)) if(!is.factor(mf[,i])) mf[,i] <- as.ts.zoo(mf[,i])
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
    if(!is.null(offset) && length(offset) != NROW(y))
	stop("Number of offsets is ", length(offset),
             ", should equal ", NROW(y), " (number of observations)")

    if (is.empty.model(mt)) {
	x <- NULL
	rval <- list(coefficients = numeric(0), residuals = y,
		  fitted.values = 0 * y, weights = w, rank = 0,
		  df.residual = length(y))
        if(!is.null(offset)) rval$fitted.values <- offset
    }
    else {
	x <- model.matrix(mt, mf, contrasts)
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
    if (model)
	rval$model <- mf
    if (ret.x)
	rval$x <- x
    if (ret.y)
	rval$y <- y
    rval$index <- index(mf[,1])

    class(rval) <- c("dynlm", class(rval))
    return(rval)
}
