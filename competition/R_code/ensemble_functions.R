#library
library(mgcv)
library(glmnet)
library(randomForest)
library(nnet)
library(foreach)
library(doParallel)
library(doRNG)
library(caret)
library(FNN)

TDC_flex_c = function(KF, obs, score = NULL, BC1 = 1, c = 1/2, lambda = 1/2) {
  if (!is.null(score)) {
    indxs = order(score, decreasing = TRUE)
    KF = KF[indxs]
    obs = obs[indxs]
  }
  nTD = cumsum(obs)
  nDD = cumsum(KF)
  fdps = (pmin(1, ((BC1 + nDD)/ nTD) * (c / (1-lambda)))) #nDD, nTD have same length
  qvals = rev(cummin(rev(fdps)))
  if (!is.null(score)) {
    qvals[indxs] = qvals
  }
  return(qvals)
}
custom_perf = function(obs, decision_values, train_alpha, BC1, c_val = 1 / 2) {
  decreasing = TRUE
  indxs = order(decision_values, decreasing = decreasing)
  obs = obs[indxs]
  decoy_wins = (obs == -1)
  target_wins = (obs == 1)
  qvals = TDC_flex_c(decoy_wins,
                     target_wins,
                     BC1=BC1,
                     c = c_val,
                     lambda = c_val)
  
  disc = sum((qvals <= train_alpha) & target_wins)
  return(disc)
}
do_custom_cv = function(f_indx, fold_indices, df, labels, positive, negative, c, train_alpha, transform, method, n_node, decay, which_var) {
  
  # Split the data into training and validation sets for the current fold
  folds = length(fold_indices)
  if (folds == 1) {
    train_indices = intersect(fold_indices[[f_indx]], which(positive | negative))
    validation_indices = fold_indices[[f_indx]]
    train_data = df[train_indices, ]
    validation_data = df[validation_indices, ]
    validation_labels = factor(labels[validation_indices])
  } else {
    train_indices = intersect(unlist(fold_indices[-f_indx]), which(positive | negative))
    validation_indices = fold_indices[[f_indx]]
    train_data = df[train_indices, ]
    validation_data = df[validation_indices, ]
    validation_labels = factor(labels[validation_indices])
  }
  
  train_data = train_data[,c(which_var, ncol(train_data))]
  validation_data = validation_data[,c(which_var, ncol(validation_data))]
  # Train the model with desired parameters and obtain decision values
  if (method == 'gam') {
    decision_values_func <- function(train_data, validation_data, transform) {
      tryCatch({
        # Fit the GAM model
        formula <- as.formula(paste('labels ~ score + ', transform, sep = ''))
        model <- mgcv::gam(formula, data = train_data, family = binomial(), drop.intercept = TRUE)
        
        # Predict using the GAM model
        predict(model, validation_data, type = 'response')
        
      }, error = function(e) {
        # In case of error, fit the RandomForest model
        model <- randomForest::randomForest(labels ~ ., data = train_data, ntree = 1000, keep.forest = TRUE)
        
        # Predict using the RandomForest model
        dval <- predict(model, validation_data, type = "response", predict.all = TRUE)$individual
        
        # Aggregate the predictions
        apply(dval, 1, function(x) sum(x == 1))
      })
    }
    decision_values <- decision_values_func(train_data, validation_data, transform)
  } else if (method == 'rf') {
    model = randomForest::randomForest(labels~., data = train_data, ntree=1000, keep.forest = TRUE)
    decision_values = predict(model, validation_data, type = "response", predict.all = TRUE)$individual
    decision_values = apply(decision_values, 1, function(x) sum(x == 1))
  } else if (method == 'nn') {
    model = nnet::nnet(labels ~ ., data=train_data, size = n_node, decay = decay, trace = FALSE, MaxNWts = 1e5)
    decision_values = predict(model, validation_data)
    decision_values = decision_values[, ncol(decision_values)]
  }
  
  # Sometimes the positive and negative classes are "swapped"
  if (var(decision_values, validation_data$score) <= 0) {
    decision_values = -decision_values
  }
  
  performance_fold = custom_perf(validation_labels, decision_values, train_alpha, BC1 = 0,
                                 c_val = c)
  
  return(list(decision_values=decision_values, performance_fold=performance_fold, validation_indices = validation_indices))
}
train_ensemble = function(labels, df, positive, negative, transform_gam='.', methods = c('gam','rf','nn'), n_nodes = c(2,5,10), 
                          decays = c(0,0.1,1), folds = 3, train_alpha = 0.1, c, verbose = FALSE,seed=1, reps = 5, which_var = which_var, num_cores = 5) {
  #set seed
  set.seed(seed)
  
  #get all methods to decide
  model_combs_gam = expand.grid(transform_gam, n_nodes = NA, decays = NA, 'gam')
  model_combs_rf = expand.grid('.', n_nodes = NA, decays = NA, 'rf')
  model_combs_nn = expand.grid('.', n_nodes = n_nodes, decays = decays, 'nn')
  model_combs = rbind(model_combs_gam, model_combs_rf, model_combs_nn)
  model_combs = model_combs[model_combs[, 4] %in% methods,]
  
  if (verbose) {
    print(model_combs)
  }
  
  mean_performances = matrix(0, nrow = nrow(model_combs), ncol = reps)
  decision_values_all = array(0, dim = c(length(labels), nrow(model_combs), reps))
  df$labels = factor(labels)
  
  # Stop any previously running clusters
  if(!is.null(getDoParWorkers()) && getDoParWorkers() > 1) {
    registerDoSEQ()
  }
  
  # Do split (ensure positive labels exist in each fold)
  fold_indices_all = vector('list', length = reps)
  for (rep in 1:reps) {
    fold_indices = caret::createFolds(factor(labels), k = folds)
    check_zero = sapply(1:folds, function(x) sum(labels[intersect(fold_indices[[x]], which(positive | negative))] == 1))
    while (any(check_zero == 0)) {
      fold_indices = caret::createFolds(factor(labels), k = folds)
      check_zero = sapply(1:folds, function(x) sum(labels[intersect(fold_indices[[x]], which(positive | negative))] == 1))
    }
    fold_indices_all[[rep]] = fold_indices
  }
  
  
  # Register the parallel backend
  cl = makeCluster(num_cores)
  registerDoParallel(cl)
  
  #registerDoRNG(seed)
  #set.seed(seed, kind = "L'Ecuyer-CMRG")
  #all_seeds = array(seed:(seed + reps*nrow(model_combs)*folds - 1), dim = c(reps, nrow(model_combs), folds))
  rng <- array(RNGseq( reps*nrow(model_combs)*folds, seed), dim = c(reps, nrow(model_combs), folds))
  
  # Export the custom function explicitly
  clusterExport(cl, list("do_custom_cv", "custom_perf", "TDC_flex_c"))
  
  res = foreach(rep = 1:reps) %:%
    foreach(i = 1:nrow(model_combs)) %:%
    foreach(j = 1:folds) %dopar% {
      rngtools::setRNG(rng[rep, i, j][[1]])
      transform = model_combs[i, 1]
      n_node = model_combs[i, 2]
      decay = model_combs[i, 3]
      method = model_combs[i, 4]
      do_custom_cv(j, fold_indices_all[[rep]], df, labels, positive, negative, c, train_alpha, transform, method, n_node, decay, which_var)
    }
  
  # Get average discoveries
  for (i in 1:nrow(model_combs)) {
    for (rep in 1:reps) {
      mean_performances[i, rep] = mean(sapply(1:folds, function(x) res[[rep]][[i]][[x]]$performance_fold))
    }
  }
  
  # Get all decision values
  for (i in 1:folds) {
    for (j in 1:nrow(model_combs)) {
      for (rep in 1:reps) {
        decision_values_all[res[[rep]][[j]][[i]]$validation_indices, j, rep] = res[[rep]][[j]][[i]]$decision_values
      }
    }
  }
  stopCluster(cl)
  registerDoSEQ()
  
  if (verbose) {
    print(mean_performances)
  }
  
  # Average discoveries on each of the rep*folds folds
  mean_performances = apply(mean_performances, 1, mean)
  
  if (verbose) {
    print(mean_performances)
  }
  
  select_indx = which.max(mean_performances)
  
  model_combs$perf = mean_performances
  
  decision_values_all = apply(as.matrix(decision_values_all[, select_indx, ]), 1, mean)
  
  return(list(perf = decision_values_all, method = model_combs[select_indx, 4], max_perf = max(mean_performances)))
}
filter_ensemble_RESET = function(W,x,Labels=NULL,methods=c('gam','rf','nn'),min_positive=50,train_alpha=0.5,
                                 seed=1,s=1/2,test_alpha = 0.1,verbose=FALSE,mult=1, folds=3, n_nodes = c(2,5,10), decays = c(0,0.1,1), reps = 10, num_cores = 5, get_nn = 10) {
  # Set seed
  set.seed(seed)
  
  # Total iterations
  total_iter = 2
  
  # Get labels
  if (is.null(Labels)) {
    Labels = rep(0, length(W))
    Labels[W > 0] = 1
    Labels[W < 0] = -1
  }
  
  # Sample half the decoy labels
  Labels_train = rep(1, length(W))
  rand_indxs = sample(c(TRUE, FALSE), size = sum(Labels == -1), replace = TRUE, prob = c(s, 1-s))
  Labels_train[which(Labels == -1)[rand_indxs]] = -1
  Labels_train[Labels == 0] = 0
  
  x = as.matrix(x)
  
  colnames(x) = paste('x', 1:ncol(x), sep = '')
  
  # Get additional info from zero-scoring hypotheses
  if (is.numeric(get_nn)) {
    # Find the k nearest neighbors
    x_init = x
    for (i in 1:length(get_nn)) {
      knn_results = get.knn(scale(x_init), get_nn[i])
      nn_indices = knn_results$nn.index
      avg_scores = apply(nn_indices, 1, function(x) sum(W[x] == 0))
      x_init = cbind(x_init, avg_scores)
    }
    colnames(x_init)[colnames(x_init) == 'avg_scores'] = paste('avg_scores', 1:length(get_nn), sep = '')
    x = x_init
  } 
  
  # Choose relevant side-info
  pval_thresh = 0.01 # Cutoff
  select_df = data.frame(score = abs(W)*Labels_train, x)
  gam_pvals = lapply(colnames(select_df)[2:ncol(select_df)], function(x) 
    tryCatch({ gam_fit = gam(formula(paste('score ~ ', 's(', x, ')')), data = select_df) #not all variables work with s()
    return(summary(gam_fit)$s.pv)
    },
    error = function(e) {
      cat("An error occurred trying to using a smoothing spline: ", conditionMessage(e), "\n")
      gam_fit = gam(formula(paste('score ~ ',  x)), data = select_df)
      return(summary(gam_fit)$p.pv[2])
    }
    )
  )
  gam_pvals = (which(gam_pvals <= pval_thresh)) + 1
  cols_selected = colnames(select_df)[gam_pvals]
  if (verbose) {
    print(paste('Using columns:', cols_selected))
  }
  if (length(cols_selected) <= 3) { # Arbitrary cutoff for computational convenience
    transform_gam = paste('s(', paste( cols_selected[!grepl('avg_scores', cols_selected)], collapse = ','),  ')', sep = '')
  } else {
    transform_gam = paste(paste('splines::ns(', cols_selected[!grepl('avg_scores', cols_selected)], ', df = 5)', sep = ''), collapse = '+')
  }
  if (any(grepl('avg_scores', cols_selected))) { # The additional variables are not compatible with s()
    additional_vars = cols_selected[grepl('avg_scores', cols_selected)]
    transform_gam = paste(c(transform_gam, additional_vars), collapse = ' + ')
  }
  
  #The chosen side information variables
  which_var = c(1, gam_pvals) 
  
  # Throw out zero-scoring
  Labels_train = Labels_train[Labels!=0]
  x = x[Labels!=0,]
  W = W[Labels!=0]
  Labels = Labels[Labels!=0]
  pseudo_Labels = Labels_train
  ensemble_df = scale(data.frame(score = abs(W), x))
  
  # Get initial positive negative training set
  negative_set = (Labels_train < 1)
  init_pos = -Inf
  for (col in which_var) {
    score_init = ensemble_df[, col]
    name_init = colnames(ensemble_df)[col]
    if (var(Labels_train == 1, score_init) < 0) {
      score_init = -score_init
    }
    q_vals_temp = TDC_flex_c(Labels_train == -1 , Labels_train == 1, score = score_init, BC1 = 0, c = (mult*(1 - s) + 1)/(mult + 1), lambda = (mult*(1 - s) + 1)/(mult + 1))
    positive_set_temp = (q_vals_temp <= test_alpha) & (Labels_train == 1)
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
  while(sum(positive_set) < min(min_positive, sum(Labels_train == 1))) {
    train_alpha_new = train_alpha_new + 0.01
    positive_set = ((q_vals <= train_alpha_new) & (Labels_train == 1)) 
  }
  
  if (verbose) {
    print(paste('best feature:', best_name))
  }
  best_pseudo = -Inf
  for (iter in 1:total_iter) {
    if (verbose) {
      print(paste('iter:', iter))
      print(paste('positive_set:', sum(positive_set)))
    }
    
    res = try(train_ensemble(Labels_train, data.frame(ensemble_df), positive_set, negative_set, transform_gam = transform_gam,
                             methods = methods, n_nodes = n_nodes, decays = decays, folds = folds, train_alpha = test_alpha, 
                             c = (mult*(1 - s) + 1)/(mult + 1), seed = seed + iter, reps = reps, verbose=verbose,
                             which_var = which_var, num_cores = num_cores))
    
    if (class(res) == 'try-error') {
      new_scores = best_scores
      pseudo_disc = -Inf
    } else {
      new_scores = res$perf
      pseudo_disc = res$max_perf
      if (verbose) {
        print(paste('selected model:', res$method))
      }
    }
    
    if (var(new_scores,ensemble_df[, 1]) <= 0) {
      new_scores = -new_scores
    }
    comp_order = order(new_scores, ensemble_df[,1], decreasing = TRUE)
    new_scores[comp_order] = length(new_scores):1
    q_vals = TDC_flex_c(Labels_train == -1 , Labels_train == 1, score = new_scores, BC1 = 0, c = (mult*(1 - s) + 1)/(mult + 1), lambda = (mult*(1 - s) + 1)/(mult + 1))
    # Obtain a new positive and negative set
    positive_set = ((q_vals <= test_alpha) & (Labels_train == 1)) 
    if (verbose) {
      print(paste('pseudo-discoveries:',pseudo_disc))
      q_vals = TDC_flex_c((Labels == -1) & (pseudo_Labels == 1), (Labels == 1) & (pseudo_Labels == 1), score = new_scores, BC1 = 1, c = 1/(mult*(1 - s) + 1), lambda = 1/(mult*(1 - s) + 1))
      print(paste('real discoveries:', sum(q_vals <= test_alpha & Labels == 1 & pseudo_Labels == 1)))
    }
    train_alpha_new = test_alpha
    while(sum(positive_set) < min(min_positive, sum(Labels_train == 1))) {
      train_alpha_new = train_alpha_new + 0.01
      positive_set = ((q_vals <= train_alpha_new) & (Labels_train == 1)) 
    }
    train_alpha_new = test_alpha
    best_pseudo = pseudo_disc
    best_scores = new_scores
  }
  q_vals = TDC_flex_c((Labels == -1) & (pseudo_Labels == 1), (Labels == 1) & (pseudo_Labels == 1), score = best_scores, BC1 = 1, c = 1/(mult*(1 - s) + 1), lambda = 1/(mult*(1 - s) + 1))
  return(list(q_vals = q_vals, Labels = Labels, pseudo_Labels = pseudo_Labels, score = best_scores))
}