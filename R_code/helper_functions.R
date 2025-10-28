split_folds = function(labels, folds) {
  fold_indices = caret::createFolds(factor(labels), k = folds)
  check_zero = sapply(1:folds, function(x)
    any(labels[fold_indices[[x]]] == 1) &
      any(labels[fold_indices[[x]]] == -1))
  while (any(check_zero == 0) & sum(labels == -1) >= folds) {
    fold_indices = caret::createFolds(factor(labels), k = folds)
    check_zero = sapply(1:folds, function(x)
      any(labels[fold_indices[[x]]] == 1) &
        any(labels[fold_indices[[x]]] == -1))
  }
  return(fold_indices)
}
train_test = function(df,
                      labels,
                      indices,
                      positive,
                      negative,
                      which_var,
                      f_indx,
                      dependent = TRUE,
                      method = 'rf',
                      decoy_int = c(1 / 2, 1),
                      target_upper = 1 / 2,
                      scale_params = NULL) {
  train_indices = intersect(unlist(indices[-f_indx]), which(positive |
                                                              negative))
  validation_indices = indices[[f_indx]]
  train_data = df[train_indices,]
  if (dependent) {
    neg_scores = runif(sum(labels[train_indices] == -1), min = decoy_int[1], max = decoy_int[2])
    neg_scores = (decoy_int[2] - neg_scores) * target_upper / (decoy_int[2] - decoy_int[1])
    neg_scores = abs(qnorm(neg_scores))
    neg_scores = (neg_scores - scale_params$score_center) / scale_params$score_scale
    train_data$score[labels[train_indices] == -1] = neg_scores
  }
  
  validation_data = df[validation_indices,]
  validation_labels = factor(labels[validation_indices])
  
  train_data = train_data[, c(which_var, ncol(train_data))]
  validation_data = validation_data[, c(which_var, ncol(validation_data))]
  
  return(
    list(
      train_data = train_data,
      validation_data = validation_data,
      validation_labels = validation_labels,
      validation_indices = validation_indices
    )
  )
}
TDC_flex_c = function(KF,
                      obs,
                      score = NULL,
                      BC1 = 1,
                      c = 1 / 2,
                      lambda = 1 / 2) {
  if (!is.null(score)) {
    indxs = order(score, decreasing = TRUE)
    KF = KF[indxs]
    obs = obs[indxs]
  }
  nTD = cumsum(obs)
  nDD = cumsum(KF)
  fdps = (pmin(1, ((BC1 + nDD) / nTD) * (c / (1 - lambda)))) #nDD, nTD have same length
  qvals = rev(cummin(rev(fdps)))
  if (!is.null(score)) {
    qvals[indxs] = qvals
  }
  return(qvals)
}
custom_perf = function(obs,
                       decision_values,
                       train_alpha,
                       BC1,
                       fdr = TRUE,
                       conf = 0.1,
                       c_val = 1 / 2) {
  decreasing = TRUE
  indxs = order(decision_values, decreasing = decreasing)
  obs = obs[indxs]
  decoy_wins = (obs == -1)
  target_wins = (obs == 1)
  qvals = TDC_flex_c(
    decoy_wins,
    target_wins,
    BC1 = BC1,
    c = c_val,
    lambda = c_val
  )
  disc = sum((qvals <= train_alpha) & target_wins)
  
  disc_fdp = NA
  if (!fdr) {
    disc_fdp = stepdownfdp::fdp_sd(cbind(length(obs):1, as.numeric(as.character(obs))), 
                               alpha = train_alpha, conf = conf, c = c_val, lambda = c_val)
    disc_fdp = length(disc_fdp$discoveries_ind)
  }
  return(list(disc = disc, disc_fdp = disc_fdp))
}
get_zeros = function(W, x, Labels,
                     get_nn,
                     thresh = 0.01) {
  if (!is.numeric(get_nn) ||
      (sum(Labels == 0) / length(W) <= thresh))
    return(x)
  
  X_scaled = scale(x)
  for (k_idx in seq_along(get_nn)) {
    k_val = get_nn[k_idx]
    nn = get.knn(X_scaled, k_val)$nn.index
    zcnt = rowSums(matrix(Labels[nn] == 0, ncol = k_val))
    x = cbind(x, zcnt)
    colnames(x)[ncol(x)] = paste0("avg_scores", k_idx)
  }
  return(x)
}
filter_predictors = function(W,
                             Labels_train,
                             x,
                             pval_thresh = 0.01,
                             verbose = TRUE) {
  select_df = data.frame(score = abs(W) * Labels_train, x)
  keep_features = character()
  
  # Loop over the predictor columns (skip the first, which is 'score')
  for (col in names(select_df)[-1]) {
    p = tryCatch({
      # try a smooth first
      f = reformulate(paste0("s(", col, ")"), response = "score")
      summary(gam(f, data = select_df))$s.pv
    },
    error = function(e) {
      # fall back to linear
      f  = reformulate(col, response = "score")
      summary(gam(f, data = select_df))$p.pv[2]
    })
    
    if (!is.na(p) && p <= pval_thresh)
      keep_features = c(keep_features, col)
  }
  if (verbose) {
    print(paste('Using columns:', keep_features))
  }
  
  if (length(keep_features) <= 3) {
    # Arbitrary cutoff for computational convenience
    transform_gam = paste('s(', paste(keep_features[!grepl('avg_scores', keep_features)], collapse = ','),  ')', sep = '')
  } else {
    transform_gam = paste(paste('splines::ns(', keep_features[!grepl('avg_scores', keep_features)], ', df = 5)', sep = ''),
                          collapse = '+')
  }
  if (any(grepl('avg_scores', keep_features))) {
    # The additional variables are not compatible with s()
    additional_vars = keep_features[grepl('avg_scores', keep_features)]
    transform_gam = paste(c(transform_gam, additional_vars), collapse = ' + ')
  }
  
  #The chosen side information variables
  which_var = c(1, which(names(select_df) %in% keep_features))
  return(list(which_var = which_var,
              transform_gam = transform_gam))
}
initialise_pos_neg = function(ensemble_df,
                              Labels_train,
                              which_var,
                              mult,
                              s,
                              test_alpha,
                              train_alpha,
                              min_positive,
                              verbose) {
  negative_set = (Labels_train == -1)
  init_pos = -Inf
  for (col in which_var) {
    score_init = ensemble_df[, col]
    name_init = colnames(ensemble_df)[col]
    if (var(Labels_train == 1, score_init) < 0) {
      score_init = -score_init
    }
    q_vals_temp = TDC_flex_c(
      Labels_train == -1 ,
      Labels_train == 1,
      score = score_init,
      BC1 = 0,
      c = (mult * (1 - s) + 1) / (mult + 1),
      lambda = (mult * (1 - s) + 1) / (mult + 1)
    )
    positive_set_temp = (q_vals_temp <= test_alpha) &
      (Labels_train == 1)
    if (sum(positive_set_temp) > init_pos) {
      positive_set = (q_vals_temp <= train_alpha) & (Labels_train == 1)
      q_vals = q_vals_temp
      init_pos = sum(positive_set_temp)
      best_scores = score_init
      best_name = name_init
    }
  }
  
  # If not enough positive labels...
  train_alpha_new = train_alpha
  while (sum(positive_set) < min(min_positive, sum(Labels_train == 1))) {
    train_alpha_new = train_alpha_new + 0.01
    positive_set = ((q_vals <= train_alpha_new) &
                      (Labels_train == 1))
  }
  
  if (verbose) {
    print(paste('best feature:', best_name))
  }
  return(list(positive_set = positive_set, negative_set = negative_set))
}
test_dependency = function(df,
                           labels,
                           dependent,
                           seed,
                           which_var,
                           model_combs,
                           scale_params,
                           cutoff = 0.1,
                           reps_RF = 1000,
                           nns = 0.01) {
  local_dependency = FALSE
  ss_obs = NA
  ss = NA
  if (dependent) {
    rng = array(RNGseq(reps_RF, seed))
    cutoff = (abs(qnorm(cutoff)) - scale_params$score_center)/scale_params$score_scale
    df_neg = df[labels == -1 & df$score >= cutoff, ] #for fast knn take only the top scoring decoys.
    nns = min(ceiling(nns*nrow(df)), 100)
    knn_results = nabor::knn(data = df_neg[, which_var[-1]],
                    query = df_neg[, which_var[-1]],
                    k = nns + 1)
    
    ss = foreach(rep = 1:reps_RF) %dopar% {
      score_temp = sample(df_neg$score)
      mean(sapply(2:(nns + 1), function(i)
        sum((
          score_temp - score_temp[knn_results$nn.idx[, i]]
        ) ^ 2)))
    }
    ss = unlist(ss)
    ss_obs = mean(sapply(2:(nns + 1), function(i)
      sum((
        df_neg$score - df_neg$score[knn_results$nn.idx[, i]]
      ) ^ 2)))
    if (sum(ss < ss_obs) <= 10) { # bootstrap p-value less than 1%
      if ('rf' %in% model_combs$type) {
        print("Possible local dependency detected, removing Random Forest")
        model_combs = model_combs[model_combs$type != 'rf', ]
      }
      local_dependency = TRUE
    }
  }
  return(list(local_dependency = local_dependency, model_combs = model_combs,
              ss_obs = ss_obs, ss = ss))
}