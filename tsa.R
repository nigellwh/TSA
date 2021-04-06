#install.packages('tseriesChaos')
#install.packages('matrixStats')
#install.packages('nonlinearTseries')
#install.packages('forecast')
require(tseriesChaos)
require(matrixStats)
require(nonlinearTseries)
require(forecast)

#set.seed(4401) #random seed for reproducibility
N = 100 #number of datasets to simulate

phi = (sqrt(5)-1)/2 #golden ratio, for cubic map

##functions to generate clean datasets

#param t in each system's function is number of points in dataset

#Lorenz system; usual params: sigma=10, rho=30, beta=8/3
lorenz <- function(sigma, rho, beta, t){
  lorenz_sys <- function(t, y, parms){
    f <- function(x, y, z, sigma, rho, beta){
      return(sigma*(y-x))
    }
    g <- function(x, y, z, sigma, rho, beta){
      return(x*(rho - z) - y)
    }
    h <- function(x, y, z, sigma, rho, beta){
      return(x*y - beta*z)
    }
    out <- list(c(f(y[1],y[2],y[3],parms[1],parms[2],parms[3]), g(y[1],y[2],y[3],parms[1],parms[2],parms[3]), 
             h(y[1],y[2],y[3],parms[1],parms[2],parms[3])))
    return(out)
  }
  inits <- runif(3, min = 0, max = 1)
  observed <- function(x){
    return(x[1] + x[2])
  }
  timestep = 0.01
  ts <- sim.cont(lorenz_sys, start=0, end=1500+t, dt=timestep, start.x=inits, parms=c(sigma,rho,beta),
                 obs.fun=observed)
  ts = tail(ts, -150001)
  L = length(ts)
  ts <- ts[seq(1, L, L/t)]
  return(ts)
}

#Cubic map; usual params: 
#[periodic] f=0, Q=0, A=1; 
#[chaotic] f=-0.8, Q=0, A=1.5;
#[HH; SNA through the Heagy-Hammel route] f=0.7, Q=0, A=1.88697;
#[S3; SNA through type-3 intermittency] f=0.35, Q=0, A=2.14;
#[2T] f=-0.18, Q=0, A=1.1
cubic <- function(f, Q, A, t){
  x <- vector(mode = "numeric", length = 100+t)
  theta <- vector(mode = "numeric", length = 100+t)
  x[1] <- abs(rnorm(1,0,1))
  theta[1] <- abs(rnorm(1,0,1))
  for (i in 2:(100+t)){
    x[i] <- f*cos(2*pi*theta[i-1]) - A*x[i-1] + (x[i-1])^3
    theta[i] <- theta[i-1] + phi - floor(theta[i-1]+phi)
  }
  y = x/6 + theta/10
  ts <- ts(c(y))
  ts = tail(ts, -100)
  if (is.finite(ts[t])){
    return(ts)
  } else {
    return(cubic(f, Q, A, t))
  }
}

#Logistic map; usual params: [periodic] r=3.5; [chaotic] r=4
logistic <- function(r, t){
  x <- vector(mode = "numeric", length = 1000+t)
  x[1] <- runif(1, min = 0, max = 1)
  for (i in 2:(1000+t)){
    x[i] <- r*x[i-1]*(1 - x[i-1])
  }
  x = tail(x, -1000)
  return(x)
}

#Trended random walk; param d gives slope of linear trend (use d = 0.3)
rw_trend <- function(d, t){
  x <- vector(mode = "numeric", length = 100+t)
  x[1] <- rnorm(1,0,1)
  for (i in 2:(100+t)){
    x[i] <- x[i-1]+d+rnorm(1,0,1)
  }
  x = tail(x, -100)
  return(x)
}

##function to create a matrix of N fn-timeseries (columns), using an array of params (including t)
make_mat <- function(fn, params){
  ts_mat <- matrix(nrow = tail(params, n=1), ncol = N)
  for (column in 1:N){
    ts_mat[, column] <- do.call(fn, as.list(params))
  }
  return(ts_mat)
}

##Add Gaussian noise in proportion to SD of each timeseries (column) in a matrix
noisify <- function(ts_mat, noise_level){
  sds <- colSds(ts_mat)
  new_mat <- ts_mat
  t <- nrow(new_mat)
  for (column in 1:N){
    sd = sds[column]*noise_level
    new_mat[,column] = new_mat[,column] + rnorm(t,0,sd)
  }
  return(new_mat)
}

#for every timeseries in a ts_mat: apply timeLag and estimateEmbeddingDim to find 
#the embedding params, then use nonLinearPrediction and compute prediction errors

#ave RMS error (in SD of current time series) from nonLinearPrediction with estimated optimal time 
#lag and embedding dimension (using timeLag and estimateEmbeddingDim); choose test set length 
#and prediction horizon (test_len and pred_hor)
predErrors <- function(ts_mat, test_len, pred_hor){
  N <- ncol(ts_mat)
  totalErrorinSD <- 0
  for (column in 1:N){
    ts <- ts_mat[, column]
    ts_len <- length(ts)
    error <- 0
    for (i in 0:(test_len - pred_hor)){
      actual_val <- ts[ts_len - test_len + pred_hor + i]
      train_ts <- ts[1:(ts_len - test_len + i)]
      init_radius <- sd(train_ts)/4
      time_lag <- timeLag(train_ts, "ami", "first.value", 0.3)
      emb_dim <- estimateEmbeddingDim(time.series = train_ts, time.lag = time_lag)
      pred_val <- nonLinearPrediction(train_ts, emb_dim, time_lag, pred_hor, init_radius, init_radius/8)
      error <- error + (actual_val - pred_val)**2
    }
    error <- sqrt(error / (test_len - pred_hor + 1))
    errorinSD <- error / sd(ts)
    totalErrorinSD <- totalErrorinSD + errorinSD
  }
  return(totalErrorinSD/N)
}

#Run simulations and NLTSA predictions for all 5x5 conditions for systems

#Generate ave RMS errors for given sys, noise level and dataset length with array of params 
#(excluding t), with proportion of dataset length for testing and prediction horizon
gen_err <- function(sys, params, noise_level, dataset_length, test_len_proportion=0.1, pred_hor=3){
  test_len <- dataset_length * test_len_proportion
  params_with_t <- append(params, dataset_length)
  ts_mat <- make_mat(sys, params_with_t)
  ts_mat <- noisify(ts_mat, noise_level)
  pred_error <- predErrors(ts_mat, test_len, pred_hor)
  return(pred_error)
}

#Generate 5x5 matrix of ave RMS errors from NLTSA by noise levels x dataset lengths, 
#for sys with array of params (excluding t), with proportion of dataset length for testing and 
#prediction horizon
gen_mat <- function(sys, params, test_len_proportion=0.1, pred_hor=3){
  error_mat <- matrix(nrow = 5, ncol = 5)
  noise_levels <- c(0, 0.1, 0.2, 0.3, 0.4)
  dataset_lengths <- c(50, 100, 200, 300, 500)
  for (i in 1:5){
    for (j in 1:5){
      t <- dataset_lengths[j]
      noise_level <- noise_levels[i]
      error_mat[i,j] <- gen_err(sys, params, noise_level, t, test_len_proportion, pred_hor)
    }
  }
  return(error_mat)
}

#for every timeseries in a ts_mat: use auto.arima to generate ARIMA models and predictions, and
#compute prediction errors
predARIMAErrors <- function(ts_mat, test_len, pred_hor){
  N <- ncol(ts_mat)
  totalErrorinSD <- 0
  for (column in 1:N){
    ts <- ts_mat[, column]
    ts_len <- length(ts)
    error <- 0
    for (i in 0:(test_len - pred_hor)){
      actual_val <- ts[ts_len - test_len + pred_hor + i]
      train_ts <- train_ts <- ts[1:(ts_len - test_len + i)]
      model <- auto.arima(train_ts)
      pred_val <- forecast(model, h = pred_hor)$mean[pred_hor]
      error <- error + (actual_val - pred_val)**2
    }
    error <- sqrt(error / (test_len - pred_hor + 1))
    errorinSD <- error / sd(ts)
    totalErrorinSD <- totalErrorinSD + errorinSD
  }
  return(totalErrorinSD/N)
}

#Run simulations and ARIMA predictions for all 5x5 conditions for systems

#Generate ave RMS errors for given sys, noise level and dataset length with array of params 
#(excluding t), with proportion of dataset length for testing and prediction horizon
gen_arima_err <- function(sys, params, noise_level, dataset_length, test_len_proportion=0.1, pred_hor=3){
  test_len <- dataset_length * test_len_proportion
  params_with_t <- append(params, dataset_length)
  ts_mat <- make_mat(sys, params_with_t)
  ts_mat <- noisify(ts_mat, noise_level)
  pred_error <- predARIMAErrors(ts_mat, test_len, pred_hor)
  return(pred_error)
}

#Generate 5x5 matrix of ave RMS errors from ARIMA by noise levels x dataset lengths, 
#for sys with array of params (excluding t), with proportion of dataset length for testing and 
#prediction horizon
gen_arima_mat <- function(sys, params, test_len_proportion=0.1, pred_hor=3){
  error_mat <- matrix(nrow = 5, ncol = 5)
  noise_levels <- c(0, 0.1, 0.2, 0.3, 0.4)
  dataset_lengths <- c(50, 100, 200, 300, 500)
  for (i in 1:5){
    for (j in 1:5){
      t <- dataset_lengths[j]
      noise_level <- noise_levels[i]
      error_mat[i,j] <- gen_arima_err(sys, params, noise_level, t, test_len_proportion, pred_hor)
    }
  }
  return(error_mat)
}

#Calculate prediction errors for one time series with given test length and
#prediction horizon (used for empirical dataset)

#NLTSA function
predtsErrors <- function(ts, test_len, pred_hor){
  ts_len <- length(ts)
  error <- 0
  for (i in 0:(test_len - pred_hor)){
    actual_val <- ts[ts_len - test_len + pred_hor + i]
    train_ts <- ts[1:(ts_len - test_len + i)]
    init_radius <- sd(train_ts)/4
    time_lag <- timeLag(train_ts, "ami", "first.value", 0.2)
    #emb_dim <- estimateEmbeddingDim(time.series = train_ts, time.lag = time_lag)
    pred_val <- nonLinearPrediction(train_ts, 3, time_lag, pred_hor, init_radius, init_radius/8)
    error <- error + (actual_val - pred_val)**2
  }
  error <- sqrt(error / (test_len - pred_hor + 1))
  errorinSD <- error / sd(ts)
  return(errorinSD)
}

#ARIMA function
predARIMAtsErrors <- function(ts, test_len, pred_hor){
  ts_len <- length(ts)
  error <- 0
  for (i in 0:(test_len - pred_hor)){
    actual_val <- ts[ts_len - test_len + pred_hor + i]
    train_ts <- train_ts <- ts[1:(ts_len - test_len + i)]
    model <- auto.arima(train_ts)
    pred_val <- forecast(model, h = pred_hor)$mean[pred_hor]
    error <- error + (actual_val - pred_val)**2
  }
  error <- sqrt(error / (test_len - pred_hor + 1))
  errorinSD <- error / sd(ts)
  return(errorinSD)
}