\name{BayesianBEKK}
\alias{BayesianBEKK}

\title{Bayesian Estimation of Bivariate Volatility Model}
\description{The Multivariate Generalized Autoregressive Conditional Heteroskedasticity (MGARCH) models are used for modelling the volatile multivariate data sets. In this package a variant of MGARCH called BEKK (Baba, Engle, Kraft, Kroner) proposed by Engle and Kroner (1995) <http://www.jstor.org/stable/3532933> has been used to estimate the bivariate time series data using Bayesian technique.}

\usage{BayesianBEKK(X,sd,iter,burnIn)}
%-
\arguments{
  \item{X}{
a bivarite time series data

}
\item{sd}{
sd is the standard deviation of the Normal priors assigned to each parameter. Default value is taken as 10.

}
\item{iter}{
iter is the number of iterations for the Markov chain Monte Carlo (MCMC) chain. Default value is 100.

}
\item{burnIn}{
burnIn is the number of iterations to be excluded from the estimate of the chain. Default value is 40.

}

}
\details{Engle and Kroner (1995) introduced the BEKK model which is the direct generalization of the univariate Generalized Autoregressive Conditional Heteroskedasticity (GARCH) model. The resulting variance is dependent on the amount of currently available information.In this package we have estimated this model using Bayesian technique. Metropilis-Hasting algorithm is used to generate the posterior density of the model parametrs. Normal distribution is used as priors for each parameter of the model following Fioruci et al., (2014).

}
\value{It returns the Bayesian estimates of the BEKK model.}
\references{
%%
Ardia, D. and Hoogerheide, L.F.(2010)<https://doi:10.32614/RJ-2010-014>

Bauwens et al.(2006)<https://doi.org/10.1002/jae.842>

Bauwens, L. and Rombouts, J.V.K.(2007)<https://doi.org/10.1080/07474930701220576>

Engle, R.F. and Kroner, K.F.(1995)<http://www.jstor.org/stable/3532933>

Fioruci et al.(2014)<https://doi.org/10.1080/02664763.2013.839635>

Metropolis et al.(1953)<https://doi.org/10.1063/1.1699114>

}
\author{
Achal Lama,Girish K Jha, K N Singh and Bishal Gurung
}
\note{This package cab be used to analyse bivariate series only}

\seealso{BEKK11}
\examples{
set.seed(6)
rt1=rnorm(20,30,5)
rt2=rnorm(20,20,1)
data_trial=cbind(rt1,rt2)
BayesianBEKK(data_trial)
}

\keyword{Bayesian}
\keyword{BEKK}
