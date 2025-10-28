#libraries
library(tidyverse)
source("../../R_code/ensemble_functions.R")
source("../../R_code/helper_functions.R")

#for parallelization
RNGkind(kind = "L'Ecuyer-CMRG")

#seed
set.seed(5022024)

INDs = c(0:4)

df_all = data.frame(discoveries = c(),
                    alpha = c(),
                    time = c(),
                    type = c()
)

for (ind in INDs) {
  df = read_delim(paste('HEK293/crux-output/narrow_', ind, '.make-pin.pin', sep = ''))
  peptide_list = read_delim(paste('HEK293/index-', ind, '/tide-index.peptides.txt', sep = ''))
  
  df = df[sample(nrow(df)), ]
  df = df[order(df$XCorr, decreasing = TRUE), ]
  df$original_target_sequence = substring(df$Peptide, 3, nchar(df$Peptide) - 2)
  df$original_target_sequence = gsub('[0-9]+|\\.|\\[|\\]', '', df$original_target_sequence)
  df$original_target_sequence[df$Label == -1] = peptide_list$target[match(df$original_target_sequence[df$Label == -1], peptide_list$`decoy(s)`)]
  
  df = df[!duplicated(df$original_target_sequence), ]
  
  W = df$XCorr - min(df$XCorr) #Can be negative so adjusting
  Labels = (df$Label)
  z = df %>% dplyr::select('deltCn', 'Charge2', 'Charge3', 'Charge4', 'Charge5', 'PepLen', 'lnNumDSP', 'dM', 'absdM')
  z = as.matrix(z)
  
  start.time = Sys.time()
  res = filter_ensemble_RESET(W, z, Labels = Labels, verbose = TRUE, test_alpha = 0.01, 
                              seed = ind + 5022024, mult = 1, reps = 10, 
                              num_cores = 20, dependent = FALSE)
  end.time = Sys.time()
  total.time = difftime(end.time, start.time, units='mins')
  power = sum(res$q_vals <= 0.01 & res$pseudo_Labels == 1 & res$Labels == 1)
  
  df_all = rbind(df_all, list(discoveries = power, alpha = 0.01, time = total.time, type = 'RESET ensemble'))
  
  start.time = Sys.time()
  q_vals = TDC_flex_c(df$Label == -1, df$Label == 1)
  end.time = Sys.time()
  total.time = difftime(end.time, start.time, units='mins')
  power = sum(q_vals <= 0.01 & df$Label == 1)
  
  df_all = rbind(df_all, list(discoveries = power, alpha = 0.01, time = total.time, type = 'KF'))
}
write.csv(df_all, 'results/HEK293_RESET_ensemble.csv')
