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
library(nabor)
library(kknn)

do_custom_cv = function(f_indx,
                        fold_indices,
                        df,
                        labels,
                        positive,
                        negative,
                        c,
                        train_alpha,
                        transform,
                        method,
                        n_node,
                        decay,
                        which_var,
                        skip_gam,
                        dependent,
                        decoy_int,
                        target_upper,
                        fdr,
                        conf,
                        scale_params) {
  # Split the data into training and validation sets for the current fold
  folds = length(fold_indices)
  indices = fold_indices
  get_train_test = train_test(
    df = df,
    labels = labels,
    indices = indices,
    positive = positive,
    negative = negative,
    which_var = which_var,
    f_indx = f_indx,
    dependent = dependent,
    method = method,
    decoy_int = decoy_int,
    target_upper = target_upper,
    scale_params = scale_params
  )
  
  train_data = get_train_test$train_data
  validation_data = get_train_test$validation_data
  validation_labels = get_train_test$validation_labels
  validation_indices = get_train_test$validation_indices
  
  # Train the model with desired parameters and obtain decision values
  if (method == 'gam') {
    formula = as.formula(paste('labels ~ score + ', transform, sep = ''))
    model = try(mgcv::gam(
      formula,
      data = train_data,
      family = binomial(),
      drop.intercept = TRUE
    ))
    if (any(class(model) == "try-error")) {
      if (skip_gam) {
        decision_values = rep(NA, nrow(validation_data))
      } else {
        vars = paste(names(train_data)[1:(ncol(train_data) - 1)], collapse = "+")
        formula = as.formula(paste('labels ~ ', vars))
        model = mgcv::gam(
          formula,
          data = train_data,
          family = binomial(),
          drop.intercept = TRUE
        )
        decision_values = predict(model, validation_data, type = 'response')
      }
    } else {
      decision_values = predict(model, validation_data, type = 'response')
    }
  } else if (method == 'rf') {
    model = randomForest::randomForest(labels ~ .,
                                       data = train_data,
                                       ntree = 1000,
                                       keep.forest = TRUE)
    decision_values = predict(model,
                              validation_data,
                              type = "response",
                              predict.all = TRUE)$individual
    decision_values = apply(decision_values, 1, function(x)
      sum(x == 1))
  } else if (method == 'nn') {
    model = nnet::nnet(
      labels ~ .,
      data = train_data,
      size = n_node,
      decay = decay,
      trace = FALSE,
      MaxNWts = 1e5
    )
    decision_values = predict(model, validation_data)
    decision_values = decision_values[, ncol(decision_values)]
  }
  
  # Sometimes the positive and negative classes are "swapped"
  if (!any(is.na(decision_values)) &&
      var(decision_values, validation_data$score) <= 0) {
    decision_values = -decision_values
  }
  
  if (any(is.na(decision_values))) {
    performance_fold = NA
    performance_fold_fdp = NA
  } else {
    res = custom_perf(
      validation_labels,
      decision_values,
      train_alpha,
      BC1 = 0,
      fdr = fdr,
      conf = conf,
      c_val = c
    )
    performance_fold = res$disc
    performance_fold_fdp = res$disc_fdp
  }
  
  return(
    list(
      decision_values = decision_values,
      performance_fold = performance_fold,
      performance_fold_fdp = performance_fold_fdp,
      validation_indices = validation_indices
    )
  )
}
train_ensemble = function(labels,
                          df,
                          positive,
                          negative,
                          transform_gam = '.',
                          methods = c('gam', 'rf', 'nn'),
                          n_nodes = c(2, 5, 10),
                          decays = c(0, 0.1, 1),
                          folds = 3,
                          train_alpha = 0.1,
                          c = 3 / 4,
                          verbose = FALSE,
                          seed = 1,
                          reps = 5,
                          which_var,
                          num_cores = 5,
                          dependent = TRUE,
                          decoy_int = c(1 / 2, 1),
                          target_upper = 1 / 2,
                          fdr = TRUE,
                          conf = NULL,
                          scale_params) {
  #set seed
  set.seed(seed)
  
  #get all methods to decide
  model_combs_gam = expand.grid(
    transformation = transform_gam,
    n_nodes = NA,
    decays = NA,
    type = 'gam'
  )
  model_combs_rf = expand.grid(
    transformation = '.',
    n_nodes = NA,
    decays = NA,
    type = 'rf'
  )
  model_combs_nn = expand.grid(
    transformation = '.',
    n_nodes = n_nodes,
    decays = decays,
    type = 'nn'
  )
  model_combs = rbind(model_combs_gam, model_combs_rf, model_combs_nn)
  model_combs = model_combs[model_combs[, 4] %in% methods,]
  
  #stop any previously running clusters
  if (!is.null(getDoParWorkers()) && getDoParWorkers() > 1) {
    registerDoSEQ()
  }
  
  #initalise random seeds
  rng = RNGseq(reps, seed)
  
  #register the parallel backend
  cl = makeCluster(num_cores)
  registerDoParallel(cl)
  clusterExport(
    cl,
    list(
      "do_custom_cv",
      "custom_perf",
      "TDC_flex_c",
      "train_test",
      "split_folds"
    )
  )
  
  #get folds
  fold_indices_all = foreach(rep = 1:reps) %dopar% {
    rngtools::setRNG(rng[rep][[1]])
    split_folds(labels, folds)
  }
  
  if (verbose) {
    print("Folds obtained.")
  }
  
  #test if there is any local dependency
  res = test_dependency(df, labels,
                  dependent,
                  seed,
                  which_var,
                  model_combs,
                  scale_params)
  local_dependency = res$local_dependency
  model_combs = res$model_combs
  skip_gam = FALSE
  ss_obs = res$ss_obs
  ss = res$ss
  
  #default to nn if model_combs is empty
  if (nrow(model_combs) == 0){
    model_combs = expand.grid(
      transformation = '.',
      n_nodes = c(2, 5, 10),
      decays = c(0, 0.1, 1),
      type = 'nn'
    )
  }
  
  #print model and parameter combinations
  if (verbose) {
    print('Using the following models:')
    print(model_combs)
  }
  
  #initialisation
  mean_performances = array(0, dim = c(nrow(model_combs),  reps))
  decision_values_all = array(0, dim = c(length(labels), nrow(model_combs), reps))
  df$labels = factor(labels)
  
  #initalise random seeds
  rng = array(RNGseq(reps * nrow(model_combs) * folds, seed),
              dim = c(reps, nrow(model_combs), folds))
  
  if (verbose) {
    print('Training starts now.')
  }
  
  res = foreach(rep = 1:reps) %:%
    foreach(i = 1:nrow(model_combs)) %:%
    foreach(j = 1:folds) %dopar% {
      rngtools::setRNG(rng[rep, i, j][[1]])
      transform = model_combs[i, 1]
      n_node = model_combs[i, 2]
      decay = model_combs[i, 3]
      method = model_combs[i, 4]
      do_custom_cv(
        j,
        fold_indices_all[[rep]],
        df,
        labels,
        positive,
        negative,
        c,
        train_alpha,
        transform,
        method,
        n_node,
        decay,
        which_var,
        skip_gam,
        dependent,
        decoy_int,
        target_upper,
        fdr,
        conf,
        scale_params
      )
    }
  
  stopCluster(cl)
  registerDoSEQ()
  
  # Get average discoveries
  if (!fdr) {
    for (i in 1:nrow(model_combs)) {
      for (rep in 1:reps) {
        mean_performances[i, rep] = mean(sapply(1:folds, function(x)
          res[[rep]][[i]][[x]]$performance_fold_fdp))
      }
    }
  }
  if (all(mean_performances == 0) | fdr) {
    for (i in 1:nrow(model_combs)) {
      for (rep in 1:reps) {
        mean_performances[i, rep] = mean(sapply(1:folds, function(x)
          res[[rep]][[i]][[x]]$performance_fold))
      }
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
  
  if (verbose) {
    print("All performances:")
    print(mean_performances)
  }
  
  mean_performances = apply(mean_performances, 1, mean)
  select_indx = which.max(mean_performances)
  
  if (verbose) {
    print("Average performances:")
    print(mean_performances)
  }
  
  decision_values_all = apply(as.matrix(decision_values_all[, select_indx, ]), 1, mean)
  
  return(list(
    perf = decision_values_all,
    method = model_combs[select_indx, 4],
    max_perf = max(mean_performances, na.rm = TRUE),
    local_dependency = local_dependency,
    ss = ss,
    ss_obs = ss_obs
  ))
}
filter_ensemble_RESET = function(W,
                                 x,
                                 decoy_int = c(1 / 2, 1),
                                 target_upper = 1 / 2,
                                 Labels = NULL,
                                 methods = c('gam', 'rf', 'nn'),
                                 min_positive = 50,
                                 train_alpha = 0.5,
                                 seed = 1,
                                 s = 1 / 2,
                                 test_alpha = 0.1,
                                 verbose = FALSE,
                                 mult = 1,
                                 folds = 3,
                                 n_nodes = c(2, 5, 10),
                                 decays = c(0, 0.1, 1),
                                 reps = 10,
                                 num_cores = 5,
                                 get_nn = 20,
                                 dependent = TRUE,
                                 fdr = TRUE,
                                 conf = NULL) {
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
  rand_indxs = sample(
    c(TRUE, FALSE),
    size = sum(Labels == -1),
    replace = TRUE,
    prob = c(s, 1 - s)
  )
  Labels_train[which(Labels == -1)[rand_indxs]] = -1
  Labels_train[Labels == 0] = 0
  
  # Convert x to dataframe
  x = data.frame(x)
  colnames(x) = paste('x', 1:ncol(x), sep = '')
  
  # Get additional info from zero-scoring hypotheses
  x = get_zeros(
    W = W,
    x = x,
    Labels = Labels,
    get_nn = get_nn
  )
  
  # Choose relevant side-info
  res = filter_predictors(W, Labels_train, x, verbose = verbose)
  which_var = res$which_var
  transform_gam = res$transform_gam
  
  # Throw out zero-scoring
  Labels_train = Labels_train[Labels != 0]
  x = x[Labels != 0, , drop = FALSE]
  W = W[Labels != 0]
  Labels = Labels[Labels != 0]
  pseudo_Labels = Labels_train
  
  # Standardize
  ensemble_df = data.frame(score = abs(W), x)
  ensemble_df = scale(ensemble_df)
  scale_params = list(
    score_center = attr(ensemble_df, "scaled:center")[[1]],
    score_scale = attr(ensemble_df, "scaled:scale")[[1]]
  )
  
  # Get initial positive negative training set
  res = initialise_pos_neg(ensemble_df,
                           Labels_train,
                           which_var,
                           mult,
                           s,
                           test_alpha,
                           train_alpha,
                           min_positive,
                           verbose)
  positive_set = res$positive_set
  negative_set = res$negative_set
  
  best_pseudo = -Inf
  # Start the two iterations
  for (iter in 1:total_iter) {
    if (verbose) {
      print(paste('iter:', iter))
      print(paste('positive_set:', sum(positive_set)))
    }
    
    res = try(train_ensemble(
      Labels_train,
      data.frame(ensemble_df),
      positive_set,
      negative_set,
      transform_gam = transform_gam,
      methods = methods,
      n_nodes = n_nodes,
      decays = decays,
      folds = folds,
      train_alpha = test_alpha,
      c = (mult * (1 - s) + 1) / (mult + 1),
      seed = seed + iter,
      reps = reps,
      verbose = verbose,
      which_var = which_var,
      num_cores = num_cores,
      dependent = dependent,
      decoy_int = decoy_int,
      target_upper = target_upper,
      fdr = fdr,
      conf = conf,
      scale_params = scale_params
    ))
    
    if (class(res) == 'try-error') {
      new_scores = best_scores
      pseudo_disc = -Inf
      local_dependency = NULL
    } else {
      new_scores = res$perf
      pseudo_disc = res$max_perf
      local_dependency = res$local_dependency
      if (verbose) {
        print(paste('selected model:', res$method))
      }
    }
    
    if (iter == 1) {
      flag = local_dependency
    }
    
    if (var(new_scores, ensemble_df[, 1]) <= 0) {
      new_scores = -new_scores
    }
    
    comp_order = order(new_scores, ensemble_df[, 1], decreasing = TRUE)
    new_scores[comp_order] = length(new_scores):1
    
    q_vals = TDC_flex_c(
      Labels_train == -1 ,
      Labels_train == 1,
      score = new_scores,
      BC1 = 0,
      c = (mult * (1 - s) + 1) / (mult + 1),
      lambda = (mult * (1 - s) + 1) / (mult + 1)
    )
    
    # Obtain a new positive and negative set
    positive_set = ((q_vals <= test_alpha) & (Labels_train == 1))
    if (verbose) {
      print(paste('pseudo-discoveries:', pseudo_disc))
      q_vals = TDC_flex_c((Labels == -1) &
                            (pseudo_Labels == 1),
                          (Labels == 1) &
                            (pseudo_Labels == 1),
                          score = new_scores,
                          BC1 = 1,
                          c = 1 / (mult * (1 - s) + 1),
                          lambda = 1 / (mult * (1 - s) + 1)
      )
      print(paste(
        'real discoveries:',
        sum(q_vals <= test_alpha & Labels == 1 & pseudo_Labels == 1)
      ))
    }
    train_alpha_new = test_alpha
    while (sum(positive_set) < min(min_positive, sum(Labels_train == 1))) {
      train_alpha_new = train_alpha_new + 0.01
      positive_set = ((q_vals <= train_alpha_new) &
                        (Labels_train == 1))
    }
    train_alpha_new = test_alpha
    best_pseudo = pseudo_disc
    best_scores = new_scores
  }
  
  q_vals = TDC_flex_c((Labels == -1) &
                        (pseudo_Labels == 1),
                      (Labels == 1) &
                        (pseudo_Labels == 1),
                      score = best_scores,
                      BC1 = 1,
                      c = 1 / (mult * (1 - s) + 1),
                      lambda = 1 / (mult * (1 - s) + 1)
  )
  return(
    list(
      q_vals = q_vals,
      Labels = Labels,
      pseudo_Labels = pseudo_Labels,
      score = best_scores,
      flag = flag
    )
  )
}