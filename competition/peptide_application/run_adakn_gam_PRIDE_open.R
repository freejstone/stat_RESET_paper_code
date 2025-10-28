args <- commandArgs(trailingOnly = TRUE)
PXID <- args[1]
#example on running Adaptive KFs

#libraries
library(randomForest)
library(tidyverse)
library(gam)
source("../../R_code/filter_randomForest.R")
source("../../R_code/filter_EM.R")
source("../../R_code/filter_gam.R")
source("../../R_code/filter_glm.R")
source("../../R_code/TDC_flex_c.R")

#seed
set.seed(5022024)

#read in example data

df_all = data.frame(discoveries = c(),
                    alpha = c(),
                    time = c(),
                    type = c(),
                    PXID = c()
)

df = read_delim(paste('data/', PXID, '/open_1_0.make-pin.pin', sep = ''))
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
print('running AdaKF GAM')
start.time = Sys.time()
res_gam = filter_gam(W, z, alpha = 0.01, reveal_prop = prop)
end.time = Sys.time()
total.time_gam = difftime(end.time, start.time, units='mins')
print(total.time_gam)

df_all = rbind(df_all, list(discoveries = res_gam$nrejs[1], alpha = 0.01, time = total.time_gam, type = 'GAM', PXID = PXID))

write.csv(df_all, paste('results/adakn_gam_power_open', PXID, '.csv', sep = ''))
