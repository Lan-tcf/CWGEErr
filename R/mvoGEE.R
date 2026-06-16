#' Unweighted GEE for multiple correlated binary outcomes in cross-sectional data with informative cluster size.
#'
#' Solves the generalized estimating equations for correlated binary
#' responses in clustered data assuming
#' using the method of quasi-least squares.
#'
#' The \code{data} must be provided in case level or equivalently in `long' format.
#'
#' @param formula a formula expression as for other regression models.
#' @param data an optional data frame containing the variables provided in
#' \code{formula}, \code{id}, \code{cluster.var} and \code{time.var}.
#' @param cluster a vector that identifies the clusters.
#' @param unit a vector that identifies the unit within a cluster.
#' @param resp.ind a vector that indicates the responses.
#' @param corr.str a character string that indicates the working correlation structure among the correlated responses.
#' Options include \code{"ind"} for independence, \code{"unstr"} for unstructured, and \code{"exch"} for exchangeable.
#' @param common.slope a character string indicating which variables in the model will
#' have a common slope for each of the responses.
#'
#' @return Returns an object of the class \code{"cwgee"}. This has components:
#' \item{call}{the matched call.}
#' \item{coefficients}{the estimated regression parameter vector of the marginal model.}
#' \item{coef.names}{the variable name of the coefficients.}
#' \item{robust.variance}{the estimated "robust" covariance matrix.}
#' \item{robust.se}{the estimated "robust" standard errors.}
#' \item{wald.chisq}{the Wald Chi-square test statistic for coefficient estimates.}
#' \item{p.value}{the p-value based on a Wald Chi-square test statistic that no covariates are statistically significant.}
#' \item{corr.matrix}{the estimated correlation matrix.}
#' \item{niter}{the number of iterations the model took to converge.}
#' \item{corr.str}{the working correlation structure assumed for the model.}
#' @author Aya Mitani
#' @examples
#' data(perio_base)
#' fitmod <- mvoGEE(formula = y ~ smoking + age + edu, data = perio_base,
#' cluster = subject, resp.ind = outcome, unit = tooth,
#' common.slope = c("smoking", "edu"), corr.str = "exch")
#' summary(fitmod)
#' @export



mvoGEE <- function(formula, data, cluster, resp.ind, unit, corr.str, common.slope = NULL){

  call <- match.call()
  mcall <- match.call(expand.dots = FALSE)
  mf <- match(c("formula", "data", "cluster", "resp.ind", "unit"), names(mcall), 0L)
  m <- mcall[c(1L, mf)]
  if (is.null(m$cluster))
    m$cluster <- as.name("id")
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Terms <- attr(m, "terms")
  y <- model.response(m)

  corr.strs <- c("ind", "exch", "unstr")
  strcheck <- match(corr.str, corr.strs, -1)
  if (strcheck < 1){
    stop("unknown working correlation structure")
  }

  cluster <- model.extract(m, "cluster")
  resp.ind <- model.extract(m, "resp.ind")
  unit <- model.extract(m, "unit")

  nresponse <- max(resp.ind)

  mterms <- attr(m, "terms")
  xvars <- as.character(attr(mterms, "variables"))[-c(1:2)]
  names(m) <- c("y", xvars, "cluster", "resp.ind", "unit")

  if (identical(xvars, common.slope)){
    xvars.expand <- NULL
    xvars.common <- xvars
  }else if (is.null(common.slope)){
    xvars.expand <- xvars
    xvars.common <- NULL
  }else{
    xvars.expand <- setdiff(xvars, common.slope)
    xvars.common <- intersect(xvars, common.slope)
  }

  dummyout <- make.dummy(data=m, resp.ind, nresponse, xvars.expand, xvars.common)[,-1]
  newxvars <- names(dummyout)

  mdat <- as.data.frame(cbind(y, dummyout, resp.ind, cluster, unit))
  mformula <- as.formula(paste("y", "~", paste(c("-1", paste(newxvars, collapse="+")), collapse="+")))
  nbeta <- length(newxvars)
  X_mat <- model.matrix(mformula, data = mdat)
  mframe <- model.frame(mformula, data = mdat)
  Y_var <- model.response(mframe)
  id <- unique.cluster <- unique(cluster)

  # use GLM to get initial beta estimates
  beta0 <- as.vector(coef(glm(mformula, data = mdat, family = poisson(link = "log"))))
  Z <- glm(mformula, data = mdat, family = poisson(link = "log"))$residuals
  #list(Y, nresponse, m, common.slope, unit, cluster, xvars, xvars.expand, xvars.common)

  cwgeeout <- mvoGEEgen(mdat, beta0, nbeta, id, Y_var, nresponse, X_mat, corr.str, Z)
  #list(mdat, beta0, nbeta, id, Y_var, nresponse, X_mat, corr.str, Z)

  coef.names <- newxvars

  result <- list()
  result$call <- match.call()
  result$coefficients <- cwgeeout$coefficients
  result$coef.names <- coef.names
  result$robust.variance <- cwgeeout$robust.variance
  result$robust.se <- cwgeeout$robust.se
  result$wald.chisq <- cwgeeout$wald.chisq
  result$p.value <- cwgeeout$p.value
  result$corr.str <- corr.str
  result$corr.matrix <- cwgeeout$R
  result$niter <- cwgeeout$niter
  class(result) <- "cwgee"
  result

}



