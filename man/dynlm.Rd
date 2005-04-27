\name{dynlm}
\alias{dynlm}
\title{Dynamic Linear Models and Time-Series Regression}
\description{
Interface to \code{\link{lm.wfit}} for fitting dynamic linear models
and time-series regression relationships.
}
\usage{dynlm(formula, data, subset, weights, na.action, method = "qr",
  model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE,
  contrasts = NULL, offset, start = NULL, end = NULL, ...)}
}

\arguments{
  \item{formula}{a \code{"formula"} describing the linear model to be fit.
    For details see below and \code{\link{lm}}.}
  \item{data}{an optional \code{"data.frame"} or time-series object (e.g.,
    \code{"ts"} or \code{"zoo"}), containing the variables
    in the model.  If not found in \code{data}, the variables are taken
    from \code{environment(formula)}, typically the environment from which
    \code{lm} is called.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}
  \item{weights}{an optional vector of weights to be used
    in the fitting process. If specified, weighted least squares is used
    with weights \code{weights} (that is, minimizing \code{sum(w*e^2)});
    otherwise ordinary least squares is used.}
  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset.  The \dQuote{factory-fresh}
    default is \code{\link{na.omit}}. Another possible value is
    \code{NULL}, no action. Note, that for time-series regression
    special methods like \code{\link{na.contiguous}}, \code{\link[zoo]{na.locf}}
    and \code{\link[zoo]{na.approx}} are available.}
  \item{method}{the method to be used; for fitting, currently only
    \code{method = "qr"} is supported; \code{method = "model.frame"} returns
    the model frame (the same as with \code{model = TRUE}, see below).}
  \item{model, x, y, qr}{logicals.  If \code{TRUE} the corresponding
    components of the fit (the model frame, the model matrix, the
    response, the QR decomposition) are returned.}
  \item{singular.ok}{logical. If \code{FALSE} (the default in S but
    not in \R) a singular fit is an error.}
  \item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{model.matrix.default}.}
  \item{offset}{this can be used to specify an \emph{a priori}
    known component to be included in the linear predictor
    during fitting.  An \code{\link{offset}} term can be included in the
    formula instead or as well, and if both are specified their sum is used.}
  \item{start}{start of the time period which should be used for fitting the model.}
  \item{end}{end of the time period which should be used for fitting the model.}
  \item{\dots}{additional arguments to be passed to the low level
    regression fitting functions.}
}

\details{
The interface and internals of \code{dynlm} are very similar to \code{\link{lm}},
but currently \code{dynlm} offers two advantages over the direct use of
\code{lm}: 1. extended formula processing, 2. preservation of time-series
attributes.

For specifying the \code{formula} of the model to be fitted, there are
additional functions available which facilitate the specification of 
dynamic models. An example would be \code{d(y) ~ L(y, 2)}, where
\code{d(x, k)} is \code{diff(x, lag = k)} and \code{L(x, k)} is
\code{lag(x, lag = -k)}, note the difference in sign. The default
for \code{k} is in both cases \code{1}.

The specification of dynamic relationships only makes sense if there is an
underlying ordering of the observations. Currently, \code{lm} offers only limited
support for such data, hence a major aim of \code{dynlm} is to preserve 
time-series properties of the data. Explicit support is currently available 
for \code{"ts"} and \code{"zoo"} series. Internally, the data is kept as a \code{"zoo"}
series and coerced back to \code{"ts"} if the original dependent variable was of
that class (and not internal \code{NA}s were created by the \code{na.action}).
}

\seealso{\code{\link[zoo]{zoo}}, \code{\link[zoo]{merge.zoo}}}

\examples{
## multiplicative SARIMA(1,0,0)(1,0,0)_12 model fitted
## to UK seatbelt data
uk <- log10(UKDriverDeaths)
dfm <- dynlm(uk ~ L(uk, 1) + L(uk, 12))
dfm
## explicitely set start and end
dfm <- dynlm(uk ~ L(uk, 1) + L(uk, 12), start = c(1975, 1), end = c(1982, 12))
dfm

## remove lag 12
dfm0 <- update(dfm, . ~ . - L(uk, 12))
anova(dfm0, dfm)

## add season term
dfm1 <- dynlm(uk ~ 1, start = c(1975, 1), end = c(1982, 12))
dfm2 <- dynlm(uk ~ season(uk), start = c(1975, 1), end = c(1982, 12))
anova(dfm1, dfm2)

plot(uk)
lines(fitted(dfm0), col = 2)
lines(fitted(dfm2), col = 4)

## Examples 7.11/7.12 from Greene (1993)
if(require(lmtest)) {
data(USDistLag)
dfm1 <- dynlm(consumption ~ gnp + L(consumption), data = USDistLag)
dfm2 <- dynlm(consumption ~ gnp + L(gnp), data = USDistLag)
plot(USDistLag[, "consumption"])
lines(fitted(dfm1), col = 2)
lines(fitted(dfm2), col = 4)
encomptest(dfm1, dfm2)
}
}

\keyword{regression}
