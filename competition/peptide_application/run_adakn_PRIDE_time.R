#example on running Adaptive KFs

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


df_all = data.frame(time = numeric(),
                    type = character(),
                    prop = numeric(),
                    PXID = character()
)


PXIDs = c('PXD008920', 'PXD008996', 'PXD010504',
          'PXD006856', 'PXD012528', 'PXD012611',
          'PXD013274', 'PXD016724', 'PXD019186',
          'PXD022257', 'PXD023571', 'PXD026895',
          'PXD002470')



reset_results = read.csv('results/RESET_ensemble_power_narrow_10_avg.csv')

for (PXID in PXIDs) {
  df = read_delim(paste('data/', PXID, '/narrow_1_0', '.make-pin.pin', sep = ''))
  peptide_list = read_delim(paste('data/', PXID, '/index-0', '/tide-index.peptides.txt', sep = ''))
  
  df = df[sample(nrow(df)), ]
  df = df[order(df$TailorScore, decreasing = TRUE), ]
  df$original_target_sequence = substring(df$Peptide, 3, nchar(df$Peptide) - 2)
  df$original_target_sequence = gsub('[0-9]+|\\.|\\[|\\]', '', df$original_target_sequence)
  df$original_target_sequence[df$Label == -1] = peptide_list$target[match(df$original_target_sequence[df$Label == -1], peptide_list$`decoy(s)`)]
  
  df = df[!duplicated(df$original_target_sequence), ]
  
  W = abs(df$TailorScore)*(df$Label)
  
  z = df %>% dplyr::select(any_of(c('deltLCn', 'deltCn', 'Charge2', 'Charge3', 'Charge4', 'Charge5', 'PepLen', 'lnNumDSP', 'dM', 'absdM', 'XCorr')))
  z = as.matrix(z)
  
  ############
  #RESET
  ############
  targets = reset_results$discoveries[reset_results$PXID == PXID & reset_results$type == 'RESET ensemble']
  decoys = floor(0.01*targets -1)
  reset_disc = mean(targets + decoys)
  
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
    
    W1 = df$TailorScore - min(df$TailorScore)
    W1 = abs(W1)*(df$Label)
    W1[abs(W1) <= 1e-3] = 0
    s0 = quantile(abs(W1)[W1!=0],reveal_prop=prop, 0.1)
    res_em = filter_EM_time(W1, z, s0 = s0)
    
    df_all = rbind(df_all, list(time = res_rf*(sum(abs(W) > quantile(abs(W)[W!=0], 0.1)) - reset_disc), type = 'RF', prop = prop, PXID = PXID))
    df_all = rbind(df_all, list(time = res_glm*(sum(abs(W) > quantile(abs(W)[W!=0], 0.1)) - reset_disc), type = 'GLM', prop = prop, PXID = PXID))
    df_all = rbind(df_all, list(time = res_gam*(sum(abs(W) > quantile(abs(W)[W!=0], 0.1)) - reset_disc), type = 'GAM', prop = prop, PXID = PXID))
    df_all = rbind(df_all, list(time = res_em*(sum(abs(W1) > s0) - reset_disc), type = 'EM', prop = prop, PXID = PXID))
  }
}

write.csv(df_all, 'results/PRIDE_AdaKF_est_times.csv')
