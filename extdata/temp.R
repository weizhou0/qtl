  mu <- fam$linkinv(eta)
      mu.eta.val <- fam$mu.eta(eta)
      mu.eta.val <- ifelse(mu.eta.val == 0, 1e-04, mu.eta.val)
      varmu <- fam$variance(mu)
      varmu <- ifelse(varmu == 0, 1e-04, varmu)

