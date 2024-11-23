filter_gam_time <- function(W, U, alpha = 0.1, offset = 1,
                       df_list = 5, reveal_prop = 0.1, 
                       mute = TRUE){
  
  ## Check the input format
  if(is.numeric(W)){
    W <- as.vector(W)
  }else{
    stop("W should be a numeric vector.")
  }
  
  if(is.numeric(U)){
    U <- as.matrix(U)
  }else{
    stop("U should be numeric.")
  }
  
  if(!is.numeric(reveal_prop)) stop("reveal_prop should be a numeric.")
  if(reveal_prop > 1) stop("reveal_prop should be a numeric between 0 and 1.")
  if(reveal_prop < 0) stop("reveal_prop should be a numeric between 0 and 1.")
  
  
  ## Extract dimensionality
  p <- length(W)
  
  ## check if U is in the correct form
  if(dim(U)[1] != p){
    if(dim(U)[2] == p){
      U <- t(U)
    }
    else{
      stop("Please check the dimensionality of the side information!")
    }
  }
  
  ## Prepare the output
  rejs <- list()
  nrejs <- rep(0,length(alpha))
  rej.path <- c()
  
  ## Initilization
  ordered_alpha <- sort(alpha, decreasing = TRUE)
  W_abs <- abs(W)
  W_sign <- (sign(W) + 1) / 2       # code the sign of W_j's as 0, 1, 1/2
  revealed_sign <- rep(1 / 2, p)
  
  all_id <- 1:p
  revealed_id <- c() 
  unrevealed_id <- all_id
  df_int <- df_list[1]              # set a intial value for df_int in case the model selection fails
  
  cutoff_val <- quantile(W_abs[W != 0], reveal_prop)
  cutoff <- p - sum(W_abs <= cutoff_val)
  count <- 0
  
  ## Update the list of revealed and unrevealed_id
  ind.reveal <- which(W_abs <= cutoff_val)
  ind.reveal <- unrevealed_id[ind.reveal]
  revealed_id <- c(revealed_id, ind.reveal)
  unrevealed_id <- all_id[-revealed_id]
  rej.path <- c(rej.path, ind.reveal)
  revealed_sign[ind.reveal] <- W_sign[ind.reveal]
  
  df_int <- df_list
  
  ## Iteratively reveal hypotheses;
  ## the order is determined by gam()
  start.time = Sys.time()
  fdphat <- compute_fdphat(W, revealed_id, offset = offset)
  
  count <- count + 1
  
  if(dim(U)[2] == 1){
    
    mdl <- suppressWarnings(gam(revealed_sign ~ ns(U, df_int) + W_abs, family = binomial()))
    
  }else{
    
    mdl <- suppressWarnings(gam(revealed_sign ~ U + W_abs, family = binomial()))
    
  }
  
  fitted.pval <- mdl$fitted.values[unrevealed_id]
  
  ## Reveal the W_j with smallest probability of being a positive
  ind.min <- suppressWarnings(which(fitted.pval == min(fitted.pval)))
  
  if(length(ind.min) == 1){
    ind.reveal <- ind.min
  }else{
    ind.reveal <- ind.min[which.min(W_abs[ind.min])]
  }
  
  ind.reveal <- unrevealed_id[ind.reveal]
  revealed_id <- c(revealed_id, ind.reveal)
  unrevealed_id <- all_id[-revealed_id]
  rej.path <- c(rej.path, ind.reveal)
  revealed_sign[ind.reveal] <- W_sign[ind.reveal]
  
  end.time = Sys.time()
  total.time = difftime(end.time, start.time, units='mins')
  return(total.time)
}

## Select the optimal degree of freedom via BIC

select_gam_mdl <- function(df_list, signw, U, W_abs){
  
  val <- c()
  
  for(df in df_list){
    
    gam_mdl <- suppressWarnings(gam(signw  ~ ns(U, df) + W_abs, family = binomial()))
    val <- c(val, BIC(gam_mdl))
    
  }
  
  opt_df <- df_list[which.min(val)]
  
  return(list(opt_df = opt_df))
  
}


compute_fdphat <- function(W, revealed_id, offset = 1) 
{
  p <- length(W)
  if (length(revealed_id) > 0) {
    fdphat <- (offset + sum(W[-revealed_id] < 0))/max(sum(W[-revealed_id] > 
                                                            0), 1)
  } else {
    fdphat <- (offset + sum(W < 0))/max(sum(W > 0), 1)
  }
  return(fdphat)
}