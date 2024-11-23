#libraries
library(tidyverse)
source("../R_code/ensemble_functions.R")

#for parallelization
RNGkind(kind = "L'Ecuyer-CMRG")

#seed
set.seed(5022024)

#read in example data

INDs = c(0:0)
PXIDs = c(
  'PXD008920', 'PXD008996', 'PXD010504',
  'PXD006856', 'PXD012528', 'PXD012611',
  'PXD013274', 'PXD016724', 'PXD019186',
  'PXD022257', 'PXD023571', 'PXD026895',
  'PXD002470'
)

df_all = data.frame(discoveries = c(),
                    alpha = c(),
                    time = c(),
                    type = c(),
                    PXID = c(),
                    seed = c()
)

seeds = c(0:9)

#these are the 'smaller' datasets
#PXD008920, PXD008996, PXD010504,
#PXD006856, PXD0012528, PXD012611,
#PXD013274, PXD016724, PXD019186,
#PXD022257,PXD023571,PXD026895,
#PXD002470

for (PXID in PXIDs) {
  df = read_delim(paste('data/', PXID, '/narrow_1_0.make-pin.pin', sep = ''))
  peptide_list = read_delim(paste('data/', PXID, '/index-0/tide-index.peptides.txt', sep = ''))
  
  df = df[sample(nrow(df)), ]
  df = df[order(df$TailorScore, decreasing = TRUE), ]
  df$original_target_sequence = substring(df$Peptide, 3, nchar(df$Peptide) - 2)
  df$original_target_sequence = gsub('[0-9]+|\\.|\\[|\\]', '', df$original_target_sequence)
  df$original_target_sequence[df$Label == -1] = peptide_list$target[match(df$original_target_sequence[df$Label == -1], peptide_list$`decoy(s)`)]
  
  df = df[!duplicated(df$original_target_sequence), ]
  
  W = df$TailorScore
  Labels = df$Label
  
  z = df %>% dplyr::select(any_of(c('deltLCn', 'deltCn', 'Charge2', 'Charge3', 'Charge4', 'Charge5', 'PepLen', 'lnNumDSP', 'dM', 'absdM', 'XCorr')))
  z = as.matrix(z)
  
  for (seed in seeds) {
    start.time = Sys.time()
    res = filter_ensemble_RESET(W, z, Labels = Labels, verbose = TRUE, test_alpha = 0.01, 
                                seed = seed + 5022024, mult = 1, reps = 10, 
                                n_nodes = c(2, 5, 10), decays = c(0, 0.1, 1), num_cores = 20, get_nn = FALSE)
    end.time = Sys.time()
    total.time = difftime(end.time, start.time, units='mins')
    power = sum(res$q_vals <= 0.01 & res$pseudo_Labels == 1 & res$Labels == 1)
    print(power)
    df_all = rbind(df_all, list(discoveries = power, alpha = 0.01, time = total.time, type = 'RESET ensemble', PXID = PXID, seed = seed))
  }
  
  start.time = Sys.time()
  q_vals = TDC_flex_c(df$Label == -1, df$Label == 1)
  end.time = Sys.time()
  total.time = difftime(end.time, start.time, units='mins')
  power = sum(q_vals <= 0.01 & df$Label == 1)
  
  df_all = rbind(df_all, list(discoveries = power, alpha = 0.01, time = total.time, type = 'KF', PXID = PXID, seed = seed))
}

write.csv(df_all, 'results/RESET_ensemble_power_narrow_10_avg.csv')
