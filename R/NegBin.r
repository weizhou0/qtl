
# from MASS

NegBin <- function (theta = 5, link = "log")
{
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  if (linktemp %in% c("log", "identity", "sqrt"))
    stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  }
  else {
    if (inherits(link, "link-glm")) {
      stats <- link
      if (!is.null(stats$name))
        linktemp <- stats$name
    }
    else stop(gettextf("\"%s\" link not available for negative binomial family; available links are \"identity\", \"log\" and \"sqrt\"",
                       linktemp))
  }
  .Theta <- theta
  env <- new.env(parent = .GlobalEnv)
  assign(".Theta", theta, envir = env)
  variance <- function(mu) mu + mu^2/.Theta
  validmu <- function(mu) all(mu > 0)
  dev.resids <- function(y, mu, wt) 2 * wt * (y * log(pmax(1, y)/mu) - (y + .Theta) * log((y + .Theta)/(mu + .Theta)))
  aic <- function(y, n, mu, wt, dev) {
    term <- (y + .Theta) * log(mu + .Theta) - y * log(mu) +
      lgamma(y + 1) - .Theta * log(.Theta) + lgamma(.Theta) -
      lgamma(.Theta + y)
    2 * sum(term * wt)
  }
  initialize <- expression({
    if (any(y < 0)) stop("negative values not allowed for the negative binomial family")
    n <- rep(1, nobs)
    mustart <- y + (y == 0)/6
  })
  simfun <- function(object, nsim) {
    ftd <- fitted(object)
    rnegbin(nsim * length(ftd), ftd, .Theta)
  }
  environment(variance) <- environment(validmu) <- environment(dev.resids) <- environment(aic) <- environment(simfun) <- env
  famname <- paste("NegBin(", format(round(theta, 4)), ")", sep = "")
  structure(list(family = famname, link = linktemp, linkfun = stats$linkfun,
                 linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids,
                 aic = aic, mu.eta = stats$mu.eta, initialize = initialize,
                 validmu = validmu, valideta = stats$valideta, simulate = simfun, theta = .Theta),
                 class = "family")
}

#*******************************************************************************



