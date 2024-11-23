#libraries
library(tidyverse)
library(stepdownfdp)
source("../R_code/ensemble_functions.R")

#for parallelization
RNGkind(kind = "L'Ecuyer-CMRG")

#seed
set.seed(5022024)

df_all = data.frame(discoveries = c(),
                    alpha = c(),
                    time = c(),
                    type = c(),
                    PXID = c(),
                    seed = c(),
                    conf = c()
)

seeds = c(0:4)

for (seed in seeds) {
  df = read_delim(paste('HEK293/crux-output/narrow_', seed, '.make-pin.pin', sep = ''))
  peptide_list = read_delim(paste('HEK293/index-', seed, '/tide-index.peptides.txt', sep = ''))
  
  df = df[sample(nrow(df)), ]
  df = df[order(df$XCorr, decreasing = TRUE), ]
  df$original_target_sequence = substring(df$Peptide, 3, nchar(df$Peptide) - 2)
  df$original_target_sequence = gsub('[0-9]+|\\.|\\[|\\]', '', df$original_target_sequence)
  df$original_target_sequence[df$Label == -1] = peptide_list$target[match(df$original_target_sequence[df$Label == -1], peptide_list$`decoy(s)`)]
  
  df = df[!duplicated(df$original_target_sequence), ]
  
  W = df$XCorr
  Labels = (df$Label)
  z = df %>% dplyr::select('deltCn', 'Charge2', 'Charge3', 'Charge4', 'Charge5', 'PepLen', 'lnNumDSP', 'dM', 'absdM')
  z = as.matrix(z)
  
  start.time = Sys.time()
  res = filter_ensemble_RESET(W, z, Labels = Labels, verbose = TRUE, test_alpha = 0.01, 
                              seed = seed + 5022024, mult = 1, reps = 10, 
                              n_nodes = c(2, 5, 10), decays = c(0, 0.1, 1), num_cores = 10, get_nn = FALSE)
  scores = rank(res$score[res$pseudo_Labels == 1], ties.method = 'random')
  labels = res$Labels[res$pseudo_Labels == 1]
  end.time = Sys.time()
  for (conf in c(0.1, 0.2, 0.5)) {
    res = fdp_sd(cbind(scores, labels), alpha = 0.01, conf = conf, c = 2/3, lambda = 2/3, procedure = 'coinflip')
    total.time = difftime(end.time, start.time, units='mins')
    power = length(res$discoveries)
    print(power)
    df_all = rbind(df_all, list(discoveries = power, alpha = 0.01, time = total.time, type = 'RESET ensemble', PXID = 'HEK293', seed = seed, conf = conf))
  }
  
  for (conf in c(0.1, 0.2, 0.5)) {
    start.time = Sys.time()
    res = fdp_sd(cbind(W, Labels), alpha = 0.01, conf = conf, procedure = 'coinflip')
    end.time = Sys.time()
    total.time = difftime(end.time, start.time, units='mins')
    power = length(res$discoveries_ind)
    
    df_all = rbind(df_all, list(discoveries = power, alpha = 0.01, time = total.time, type = 'KF', PXID = 'HEK293', seed = seed, conf = conf))
  }
}

write.csv(df_all, 'results/RESET_ensemble_HEK293_power_FDP.csv')
