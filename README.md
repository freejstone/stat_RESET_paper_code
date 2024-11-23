### RESET repository

This repository: contains the code used for our paper (link to be added). A package for our method is on the way! A description of the files/folders are as follows.

* `competition/R_code`: contains our RESET code plus the AdaKO code from [here](https://github.com/zhimeir/adaptiveKnockoffs) however we fixed a small issue that made the procedure slightly more conservative as well as made the `reveal_prop` argument consistent for each of the filters.

The `competition` folder: contains code related to the competition-based testing part of our paper. For our competition simulations:
* `simulation1_vary_ensemble/adakn.R`: contains the code for AdaKO plus KO in Simulation 1 in our paper
* `simulation1_vary_ensemble/RESET.R`: contains the code for RESET in Simulation 1 in our paper
* `simulation1_vary_ensemble/RESET_adakn.sh` runs the above two scripts
* `simulation1_vary_ensemble_beta/adakn.R`: contains the code for AdaKO plus KO in Simulation 2 in our paper
* `simulation1_vary_ensemble_beta/adakn.sh` runs the above script
* `simulation1_vary_ensemble_beta/RESET.R`: contains the code for RESET in Simulation 2 in our paper
* `simulation1_vary_ensemble_beta/RESET.sh` runs the above script
* `simulation2_vary_ensemble/adakn.R`: contains the code for AdaKO plus KO in Simulation 3 in our paper
* `simulation2_vary_ensemble/RESET.R`: contains the code for RESET in Simulation 3 in our paper
* `simulation2_vary_ensemble/RESET_adakn.sh` runs the above two scripts
* `simulation3_ensemble_benchtime_ind/adakn.R`: contains the code for AdaKO in Simulation 4 in our paper
* `simulation3_ensemble_benchtime_ind/adakn.sh` runs the above AdaKO script
* `simulation3_ensemble_benchtime_ind/RESET.R`: contains the code for RESET in Simulation 4 in our paper
* `simulation3_ensemble_benchtime_ind/RESET.sh` runs the above RESET script
* `simulation3_ensemble_benchtime_ind/KO.R`: contains the code for KO in Simulation 4 in our paper
* `simulation3_ensemble_benchtime_ind/KO.sh` runs the above KO script
(The naming of folders is not the greatest. I named them before the paper was assembled...)

For our `competition` real data experiments:
* `peptide_application/run_adakn_HEK293_time.R`: contains the code for getting the estimated run times for AdaKO on the HEK293 dataset
* `peptide_application/RESET_ensemble_HEK293_FDR.R`: contains the code for getting the run times and discoveries for RESET on the HEK293 dataset
* `peptide_application/run_adakn_PRIDE.R`: contains the code for getting the run times and discoveries for AdaKO on 13 of the PRIDE datasets
* `peptide_application/run_adakn_PRIDE_time.R`: contains the code for getting the estimated run times for AdaKO on 13 of the the PRIDE datasets
* `peptide_application/RESET_ensemble_PRIDE_FDR.R`: contains the code for getting the run times and discoveries for RESET on 13 of the PRIDE datasets
* `peptide_application/RESET_ensemble_HEK293_FDP.R`: contains the code for getting the discoveries for RESET with FDP on the HEK293 dataset
* `peptide_application/RESET_ensemble_PRIDE_FDP_narrow.R`: contains the code for getting the discoveries for RESET with FDP on all the PRIDE datasets
* `peptide_application/RESET_ensemble_PRIDE_FDP_open.R`: contains the code for getting the discoveries for RESET with FDP on all the PRIDE datasets using an open search
* `peptide_application/filter_x.R` where `x` denotes the type of AdaKO filter: contains the code for determining how long each update takes for AdaKO

* `extras/delta_comparison.R`: contains the analysis of the delta upper bound in the discussion.

The `pvalue` folder: contains code related to the p-value-based testing part of our paper. For our p-value simulations:
* `simulations/adapt_paper_simulations/produce_data1.R`: contains the code for producing the simulation for Simulation 5 (the script is based on code from the [AdaPT](https://github.com/lihualei71/adaptPaper) repository)
* `simulations/adapt_paper_simulations/produce_data2.R`: contains the code for producing the simulation for Simulation 6 (the script is based on code from the [AdaPT](https://github.com/lihualei71/adaptPaper) repository)
* `simulations/adapt_paper_simulations/RESET_ensemble.R`: contains the code for running RESET on Simulation 5 and 6
* `simulations/adapt_paper_simulations/RESET_ensemble.sh` runs the above script
* `simulations/adapt_paper_simulations/AdaPT.R`: contains the code for running AdaPT on Simulation 5 and 6
* `simulations/adapt_paper_simulations/AdaPT.sh` runs the above script using AdaPT GAM
* `simulations/adapt_paper_simulations/AdaPT_glmnet.sh` runs the above script using AdaPT GLMnet
* `simulations/adapt_paper_simulations/AdaPTg.R`: contains the code for running AdaPTg on Simulation 5 and 6
* `simulations/adapt_paper_simulations/AdaPTg.sh` runs the above script
* `simulations/adapt_paper_simulations/AdaPT-GMM.R`: contains the code for running AdaPTg-GMM on Simulation 5 and 6
* `simulations/adapt_paper_simulations/AdaPT-GMM.sh` runs the above script
* `simulations/adapt_paper_simulations/ZAP.R`: contains the code for running ZAP on Simulation 5 and 6
* `simulations/adapt_paper_simulations/AdaFDR.py`: contains the code for running AdaFDR on Simulation 5

For our `pvalue` real data experiments:
* `real_data/run_RESET_ensemble.R`: contains the code for running RESET on all real datasets except for the estrogen data
* `real_data/run_estrogen_RESET_ensemble.R`:  contains the code for running RESET on the estrogen data
* `real_data/run_AdaPT.R`: contains the code for running AdaPT on all real datasets except for the estrogen data
* `real_data/run_estrogen_AdaPT.R`: contains the code for running AdaPT on the estrogen data (the script is based on code from the [AdaPT](https://github.com/lihualei71/adaptPaper) repository)
* `real_data/run_AdaPTg.R`: contains the code for running AdaPTg on all real datasets except for the estrogen data
* `real_data/run_estrogen_AdaPTg.R`: contains the code for running AdaPTg on the estrogen data 
* `real_data/run_all_exp_AdaPTGMM.R`: contains the code for running AdaPTg-GMM on all real datasets except for the estrogen data (note that this script is taken from [here](https://github.com/patrickrchao/AdaPTGMM_Experiments/tree/main/AdaFDR_experiments)
* `real_data/run_estrogen_AdaPTGMM.R`: contains the code for running AdaPTg-GMM on the estrogen data
* `real_data/run_ZAP.R`: contains the code for running ZAP on all real datasets except for the estrogen data
* `real_data/run_estrogen_zap.R`: contains the code for running ZAP on the estrogen data
* `real_data/run_adafdr.py`: contains the code for running AdaFDR on all real datasets except for the estrogen data (note that this script is based on the scripts [here](https://github.com/martinjzhang/AdaFDRpaper/tree/master/vignettes))
* `real_data/run_adafdr_estrogen.py`: contains the code for running AdaFDR on the estrogen data
* `real_data/run_RESET_ensemble_FDP.R`: contains the code for running RESET on all real datasets with FDP control except for the estrogen data
* `real_data/run_estrogen_RESET_ensemble_fdp.R`:  contains the code for running RESET on the estrogen data with FDP control

* `functions/guo.R`: GR-SD method (since there is no package that we are aware of). See the paper [here](https://web.njit.edu/~wguo/Guo%20&%20Romano%202007.pdf).
* `functions/summarize_methods.R`: Helper functions from the [AdaPT](https://github.com/lihualei71/adaptPaper) repository)
