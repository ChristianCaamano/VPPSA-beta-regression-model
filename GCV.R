### GCV Model Betareg package
GCV_BR<- function(formula, data, p, k = 10, seed = 123) {
  set.seed(seed)
  folds <- createFolds(data[[as.character(formula[[2]])]], k = k)
  mse_folds <- numeric(k)
  
  for (i in seq_along(folds)) {
    train_data <- data[-folds[[i]], ]
    test_data <- data[folds[[i]], ]
    
    model <- betareg(formula, data = train_data)
    pred <- predict(model, newdata = test_data, type = "response")
    n_train <- nrow(train_data)
    mse_folds[i] <- mean((test_data[[as.character(formula[[2]])]] - pred)^2)
  }
  
  mse_total <- mean(mse_folds)
  n_total <- nrow(data)
  gcv <- mse_total / (1 - p/n_total)^2  
  return(gcv)
}

### GCV Model gamlss package

gcv_gamlss <- function(data, k = 20, p = 4, seed = 123) {
    set.seed(seed)
  folds <- createFolds(data$Y, k)
  mse_folds <- numeric(k)
  for (i in seq_along(folds)) {
    train_data <- data[-folds[[i]], ]
    test_data <- data[folds[[i]], ]

    vars_needed <- c("readscr", "avginc", "mealpct", "calwpct")
    has_variation <- all(sapply(train_data[vars_needed], function(x) length(unique(x)) > 1))
    
    if (!has_variation) {
      warning(sprintf("Fold %d omitted due to lack of variation in some variable", i))
      mse_folds[i] <- NA
      next
    }

    model <- try(gamlss(
      formula =  Y ~ readscr + cs(avginc),
      sigma.formula = ~ mealpct + cs(calwpct),
      family = BE(mu.link = "log", sigma.link = "log"),
      data = train_data
    ), silent = TRUE)

    if (inherits(model, "try-error")) {
      warning(sprintf("Model failed to fold %d", i))
      mse_folds[i] <- NA
      next
    }

    pred <- try(predict(model, newdata = test_data, type = "response", what = "mu"), silent = TRUE)
    
    if (inherits(pred, "try-error")) {
      warning(sprintf("Prediction failed on the fold %d", i))
      mse_folds[i] <- NA
    } else {
      mse_folds[i] <- mean((test_data$Y - pred)^2)
    }
  }

  n_total <- nrow(data)
  mse_mean <- mean(mse_folds, na.rm = TRUE)
  gcv <- mse_mean / (1 - p / n_total)^2
  return(gcv)
}



### pseudo R^2 gamlss
R2gamlss<- function(model) {
  n <- length(model$y)
  null_model <-  gamlss(formula = Y ~ 1,
                  sigma.formula = ~1,
                  family=BE(mu.link = "log", sigma.link = "log"))
  
  D_null <- deviance(null_model)
  D_model <- deviance(model)
  
  R2 <- (1 - exp((D_model - D_null)/n))
  
  return(R2)
}

