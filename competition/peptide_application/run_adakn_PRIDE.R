args <- commandArgs(trailingOnly = TRUE)
PXID <- args[1]
#example on running Adaptive KFs

#libraries
library(randomForest)
library(tidyverse)
library(gam)
source("../R_code/filter_randomForest.R")
source("../R_code/filter_EM.R")
source("../R_code/filter_gam.R")
source("../R_code/filter_glm.R")
source("../R_code/TDC_flex_c.R")

#seed
set.seed(5022024)

#read in example data

df_all = data.frame(discoveries = c(),
                    alpha = c(),
                    time = c(),
                    type = c(),
                    PXID = c()
)

df = read_delim(paste('data/', PXID, '/narrow_1_0.make-pin.pin', sep = ''))
peptide_list = read_delim(paste('data/', PXID, '/index-0/tide-index.peptides.txt', sep = ''))

df = df[sample(nrow(df)), ]
df = df[order(df$TailorScore, decreasing = TRUE), ]
df$original_target_sequence = substring(df$Peptide, 3, nchar(df$Peptide) - 2)
df$original_target_sequence = gsub('[0-9]+|\\.|\\[|\\]', '', df$original_target_sequence)
df$original_target_sequence[df$Label == -1] = peptide_list$target[match(df$original_target_sequence[df$Label == -1], peptide_list$`decoy(s)`)]

df = df[!duplicated(df$original_target_sequence), ]

W = abs(df$TailorScore)*(df$Label)

z = df %>% dplyr::select(any_of(c('deltLCn', 'deltCn', 'Charge2', 'Charge3', 'Charge4', 'Charge5', 'PepLen', 'lnNumDSP', 'dM', 'absdM', 'XCorr')))
z = as.matrix(z)

#########
#filters
#########
print('running AdaKF filters now')

prop = 0.1
print('running AdaKF EM')
W1 = df$TailorScore - min(df$TailorScore)
W1 = abs(W1)*(df$Label)
W1[abs(W1) <= 1e-3] = 0
s0 = quantile(abs(W1)[W1!=0],reveal_prop=prop, 0.1)
start.time = Sys.time()
res_em = filter_EM(W1, z, alpha = 0.01, s0 = s0)
end.time = Sys.time()
total.time_em = difftime(end.time, start.time, units='mins')
print(total.time_em)
print('running AdaKF GLM')
start.time = Sys.time()
res_glm = filter_glm(W, z, alpha = 0.01, reveal_prop = prop)
end.time = Sys.time()
total.time_glm = difftime(end.time, start.time, units='mins')
print(total.time_glm)
print('running AdaKF GAM')
start.time = Sys.time()
res_gam = filter_gam(W, z, alpha = 0.01, reveal_prop = prop)
end.time = Sys.time()
total.time_gam = difftime(end.time, start.time, units='mins')
print(total.time_gam)
print('running AdaKF RF')
start.time = Sys.time()
res_rf = filter_randomForest(W, z, alpha = 0.01, reveal_prop = prop)
end.time = Sys.time()
total.time_rf = difftime(end.time, start.time, units='mins')
print(total.time_rf)

df_all = rbind(df_all, list(discoveries = res_rf$nrejs[1], alpha = 0.01, time = total.time_rf, type = 'RF', PXID = PXID))
df_all = rbind(df_all, list(discoveries = res_glm$nrejs[1], alpha = 0.01, time = total.time_glm, type = 'GLM', PXID = PXID))
df_all = rbind(df_all, list(discoveries = res_gam$nrejs[1], alpha = 0.01, time = total.time_gam, type = 'GAM', PXID = PXID))
df_all = rbind(df_all, list(discoveries = res_em$nrejs[1], alpha = 0.01, time = total.time_em, type = 'EM', PXID = PXID))

write.csv(df_all, paste('results/adakn_power_narrow_', PXID, '.csv', sep = ''))
