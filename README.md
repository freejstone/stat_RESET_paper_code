### RESET repository

This repository: contains the code used for our [paper](https://arxiv.org/abs/2411.15771). A package for our method is on the way! A description of the files/folders are as follows. The data for the peptide application is very large, so we have omitted it here. Please contact me if you would like it.

* `R_code`: contains our RESET code plus the AdaKO code from [here](https://github.com/zhimeir/adaptiveKnockoffs) where we fixed a small bug that made the procedure slightly more conservative as well as made the `reveal_prop` argument consistent for each of the filters for our convenience.

The `competition` folder: contains code related to the competition-based testing part of our paper. For our competition simulations:
* `simulation1/adakn.R`: contains the code for AdaKO plus KO in Simulation 1 in our paper
* `simulation1/adakn.sh` runs the above script
* `simulation1/RESET.R`: contains the code for RESET in Simulation 1 in our paper
* `simulation1/RESET.sh` runs the above script
* `simulation2/adakn.R`: contains the code for AdaKO plus KO in Simulation 2 in our paper
* `simulation2/adakn.sh` runs the above script
* `simulation2/RESET.R`: contains the code for RESET in Simulation 2 in our paper
* `simulation2/RESET.sh` runs the above script
* `simulation3/adakn.R`: contains the code for AdaKO plus KO in Simulation 3 in our paper
* `simulation3/adakn.sh` runs the above script
* `simulation3/RESET.R`: contains the code for RESET in Simulation 3 in our paper
* `simulation3/RESET.sh` runs the above script
* `simulation4/adakn.R`: contains the code for AdaKO in Simulation 4 in our paper
* `simulation4/adakn.sh` runs the above AdaKO script
* `simulation4/RESET.R`: contains the code for RESET in Simulation 4 in our paper
* `simulation4/RESET.sh` runs the above RESET script
* `simulation4/KO.R`: contains the code for KO in Simulation 4 in our paper
* `simulation4/KO.sh` runs the above KO script

For our `competition` real data experiments using mass spectrometry based proteomics data:
* `peptide_application/run_adakn_HEK293_time.R`: contains the code for getting the estimated run times for AdaKO on the HEK293 dataset using a narrow-search
* `peptide_application/RESET_ensemble_HEK293_FDR.R`: contains the code for getting the run times and discoveries for RESET on the HEK293 dataset using a narrow-search
* `peptide_application/run_adakn_PRIDE.R`: contains the code for getting the run times and discoveries for AdaKO on 13 of the PRIDE datasets using a narrow-search
* `peptide_application/run_adakn_gam_PRIDE_open.R`: contains the code for the extra results on the run times and discoveries for AdaKO GAM on 13 of the PRIDE datasets using an open-search (only did AdaKO GAM because it is too slow to repeat the analysis using an open-search for all other AdaKO methods).
* `peptide_application/run_adakn_PRIDE_time.R`: contains the code for getting the estimated run times for AdaKO on 13 of the the PRIDE datasets using a narrow-search
* `peptide_application/RESET_ensemble_PRIDE_FDR_narrow.R`: contains the code for getting the run times and discoveries for RESET on 13 of the PRIDE datasets using a narrow-search
* `peptide_application/RESET_ensemble_PRIDE_FDR_open.R`: contains the code for getting the run times and discoveries for RESET on 13 of the PRIDE datasets using an open-search
* `peptide_application/RESET_ensemble_HEK293_FDP.R`: contains the code for getting the discoveries for RESET with FDX on the HEK293 dataset using a narrow-search
* `peptide_application/RESET_ensemble_PRIDE_FDP_narrow.R`: contains the code for getting the discoveries for RESET with FDX on all the PRIDE datasets using a narrow-search
* `peptide_application/RESET_ensemble_PRIDE_FDP_open.R`: contains the code for getting the discoveries for RESET with FDX on all the PRIDE datasets using an open-search
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

* `simulations/adadetect_no_control/produce_data.R`: contains the code for producing the simulation that demonstrates AdaDetect does not control the FDR in multiple testing.
* `simulations/adadetect_no_control/AdaDetect_RESET.R`: contains the code for running RESET with AdaDetect's empirical p-values.
* `simulations/adadetect_no_control/RESET_AdaDetect.sh`: runs the above script.
* `simulations/adadetect_no_control/RESET.R`: contains the code for running RESET.
* `simulations/adadetect_no_control/RESET.R`: runs the above script.

* `simulations/depedency/produce_data1.R`: contains the code for producing the simulation for Simulation 7.
* `simulations/depedency/produce_data2.R`: contains the code for producing the simulation for Simulation 8.
* `simulations/depedency/RESET_ensemble.R`: contains the code for running RESET on Simulation 7 and 8.
* `simulations/depedency/RESET_ensemble.sh` runs the above script.
* `simulations/depedency/AdaPT.R`: contains the code for running AdaPT on Simulation 7 and 8.
* `simulations/depedency/AdaPT.sh` runs the above script using AdaPT GAM
* `simulations/depedency/AdaPTg.R`: contains the code for running AdaPTg on Simulation 7 and 8.
* `simulations/depedency/AdaPTg.sh` runs the above script
* `simulations/depedency/AdaPT-GMM.R`: contains the code for running AdaPTg-GMM on Simulation 7 and 8.
* `simulations/depedency/AdaPT-GMM.sh` runs the above script
* `simulations/depedency/ZAP.R`: contains the code for running ZAP on Simulation 7 and 8.
* `simulations/depedency/AdaFDR.py`: contains the code for running AdaFDR on Simulation 7.
* `simulations/depedency/AdaFDR_2.py`: contains the code for running AdaFDR on Simulation 8.

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
* `real_data/run_RESET_ensemble_FDP.R`: contains the code for running RESET on all real datasets with FDX control except for the estrogen data
* `real_data/run_estrogen_RESET_ensemble_fdp.R`:  contains the code for running RESET on the estrogen data with FDX control
* `real_data/airway_repeat.R`:  contains the code for running RESET on the Airway data 100 times.

* `functions/guo.R`: GR-SD method (since there is no package that we are aware of). See the paper [here](https://web.njit.edu/~wguo/Guo%20&%20Romano%202007.pdf).
* `functions/summarize_methods.R`: Helper functions from the [AdaPT](https://github.com/lihualei71/adaptPaper) repository)
