\name{M1Germany}
\alias{M1Germany}
\encoding{latin1}
\title{German M1 Money Demand}
\description{
German M1 money demand.
}

\usage{data(M1Germany)}

\format{
\code{M1Germany} is a \code{"zoo"} series containing 4 quarterly
time series  from 1960(1) to 1996(3).
 \describe{
    \item{logm1}{logarithm of real M1 per capita,}
    \item{logprice}{logarithm of a price index,}
    \item{loggnp}{logarithm of real per capita gross national product,}
    \item{interest}{long-run interest rate,}
}
}

\details{
This is essentially the same data set as \code{\link[strucchange]{GermanM1}}, 
the important difference is that it is stored as a \code{\link[zoo]{zoo}} series 
and not as a data frame. It does not contain differenced and lagged versions
of the variables (as \code{GermanM1}) does, because these do not have to be 
computed explicitely before applying \code{dynlm}.

The (short) story behind the data is the following (for more detailed information
see \code{\link[strucchange]{GermanM1}}):
Lütkepohl et al. (1999) investigate the linearity and
stability of German M1 money demand: they find a stable regression relation
for the time before the monetary union on 1990-06-01 but a clear structural
instability afterwards. Zeileis et al. (2005) re-analyze this data set
in a monitoring situation.
}

\source{The data is provided by the German central bank and is
available online in the data archive of the Journal of Applied
Econometrics
\url{http://qed.econ.queensu.ca/jae/1999-v14.5/lutkepohl-terasvirta-wolters/}.}

\references{
Lütkepohl H., Teräsvirta T., Wolters J. (1999), Investigating
Stability and Linearity of a German M1 Money Demand Function,
\emph{Journal of Applied Econometrics}, \bold{14}, 511--525.

Zeileis A., Leisch F., Kleiber C., Hornik K. (2005), Monitoring
Structural Change in Dynamic Econometric Models,
\emph{Journal of Applied Econometrics}, \bold{20}, 99--121.
}

\seealso{\code{\link[strucchange]{GermanM1}}}

\examples{
data("M1Germany")
## fit the model of Lütkepohl et al. (1999) on the history period
## before the monetary unification 
histfm <- dynlm(d(logm1) ~ d(L(loggnp, 2)) + d(interest) + d(L(interest)) + d(logprice) +
                           L(logm1) + L(loggnp) + L(interest) +
                           season(logm1, ref = 4),
                data = M1Germany, start = c(1961, 1), end = c(1990, 2))

## fit on extended sample period
fm <- update(histfm, end = c(1995, 4))

if(require("strucchange")) {
  scus <- gefp(fm, fit = NULL)
  plot(scus, functional = supLM(0.1))
}
}

\keyword{datasets}
