#libraries
library(adaptiveKnockoff)
library(randomForest)
library(tidyverse)
library(gam)
source("filter_randomForest.R")
source("filter_EM.R")
source("filter_gam.R")
source("filter_glm.R")

#seed
set.seed(5022024)


INDs = c(0:4)

df_all = data.frame(time = numeric(),
                    type = character(),
                    prop = numeric(),
                    ind = numeric()
)

reset_results = read.csv('results/HEK293_RESET_ensemble.csv')
target_discoveries = reset_results$discoveries[reset_results$type == 'RESET ensemble']

for (ind in INDs) {
  df = read_delim(paste('HEK293/crux-output/narrow_', ind, '.make-pin.pin', sep = ''))
  peptide_list = read_delim(paste('HEK293/index-', ind, '/tide-index.peptides.txt', sep = ''))
  
  df = df[sample(nrow(df)), ]
  df = df[order(df$XCorr, decreasing = TRUE), ]
  df$original_target_sequence = substring(df$Peptide, 3, nchar(df$Peptide) - 2)
  df$original_target_sequence = gsub('[0-9]+|\\.|\\[|\\]', '', df$original_target_sequence)
  df$original_target_sequence[df$Label == -1] = peptide_list$target[match(df$original_target_sequence[df$Label == -1], peptide_list$`decoy(s)`)]
  
  df = df[!duplicated(df$original_target_sequence), ]
  
  W = df$XCorr - min(df$XCorr)
  W = abs(W)*(df$Label)
  z = df %>% dplyr::select('deltCn', 'Charge2', 'Charge3', 'Charge4', 'Charge5', 'PepLen', 'lnNumDSP', 'dM', 'absdM')
  z = as.matrix(z)
  
  ############
  #RESET
  ############
  targets = target_discoveries[ind + 1]
  decoys = floor(0.01*targets -1)
  reset_disc = targets + decoys
  
  ###########
  #obtain k
  ###########
  k = floor( (nrow(df) - reset_disc)*10/nrow(df) )
  
  ############
  #run filters
  ############
  props = (1:k)/10
  for (prop in props) {
    print(paste('running:', prop))
    res_rf = filter_randomForest_time(W, z, reveal_prop = prop)
    res_glm = filter_glm_time(W, z, reveal_prop = prop)
    res_gam = filter_gam_time(W, z, reveal_prop = prop)
    
    W1 = df$XCorr - min(df$XCorr)
    W1 = abs(W1)*(df$Label)
    W1[abs(W1) <= 1e-3] = 0
    s0 = quantile(abs(W1)[W1!=0], 0.1)
    res_em = filter_EM_time(W1, z, s0 = s0)
    
    df_all = rbind(df_all, list(time = res_rf*(sum(abs(W) > quantile(abs(W)[W!=0], 0.1)) - reset_disc), type = 'RF', prop = prop, ind = ind))
    df_all = rbind(df_all, list(time = res_glm*(sum(abs(W) > quantile(abs(W)[W!=0], 0.1)) - reset_disc), type = 'GLM', prop = prop, ind = ind))
    df_all = rbind(df_all, list(time = res_gam*(sum(abs(W) > quantile(abs(W)[W!=0], 0.1)) - reset_disc), type = 'GAM', prop = prop, ind = ind))
    df_all = rbind(df_all, list(time = res_em*(sum(abs(W1) > s0) - reset_disc), type = 'EM', prop = prop, ind = ind))
  }
}

write.csv(df_all, 'results/HEK293_AdaKF.csv')
