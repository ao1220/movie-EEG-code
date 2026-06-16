# README

This file describes the input and output of each file.

## `readANOVAorgData.m`

### Input

This function uses the following input arguments:

* `N1`: number of subjects in group 1.
* `N2`: number of subjects in group 2.
* `numofchan`: number of EEG channels.
* `numofbands`: number of frequency bands.

The following folder paths need to be set inside the function before running:

* `inputpath1_1`: group 1, condition 1 data folder.
* `inputpath1_2`: group 1, condition 2 data folder.
* `inputpath2_1`: group 2, condition 1 data folder.
* `inputpath2_2`: group 2, condition 2 data folder.

Each input folder should contain `.mat` files. Each `.mat` file should contain:

* `EEG_results.M_zscore_mean`
* `EEG_results.parameter.chanlocs`

### Output

The function returns:

* `network1_1`: group 1, condition 1 network data.
* `network1_2`: group 1, condition 2 network data.
* `network2_1`: group 2, condition 1 network data.
* `network2_2`: group 2, condition 2 network data.
* `y1_1`: reshaped group 1, condition 1 data.
* `y1_2`: reshaped group 1, condition 2 data.
* `y2_1`: reshaped group 2, condition 1 data.
* `y2_2`: reshaped group 2, condition 2 data.
* `chan`: EEG channel location information.

## `anova_code.m`

### Input

This script uses either raw data or a saved MATLAB data file.

If `readOrgDataFlag = true`, the script calls:

* `readANOVAorgData(N1, N2, numofchan, numofbands)`

In this case, the input data folders are set in `readANOVAorgData.m`.

If `readOrgDataFlag = false`, the script reads:

* `orgNetData.mat`

The saved input `.mat` file should contain:

* `network1_1`: group 1, condition 1 network data.
* `network1_2`: group 1, condition 2 network data.
* `network2_1`: group 2, condition 1 network data.
* `network2_2`: group 2, condition 2 network data.
* `y1_1`: reshaped group 1, condition 1 data.
* `y1_2`: reshaped group 1, condition 2 data.
* `y2_1`: reshaped group 2, condition 1 data.
* `y2_2`: reshaped group 2, condition 2 data.

### Output

This script creates one result folder:

* `results/<selected_effect>/`

The selected effect folder is based on `flagOfConcern`:

* `patientEffect`
* `movieEffect`
* `interEffect`

The script saves one line-data `.mat` file:

* `line_positive.mat` or `line_negative.mat`

This file contains:

* `line_mci_m1_PorN`: selected MCI movie 1 edge data.
* `line_mci_m2_PorN`: selected MCI movie 2 edge data.
* `line_scd_m1_PorN`: selected SCD movie 1 edge data.
* `line_scd_m2_PorN`: selected SCD movie 2 edge data.

The script also saves one ANOVA result `.mat` file:

* `anova_Results_positive.mat` or `anova_Results_negative.mat`

This file contains the structure `Result`, including:

* `Result.f`: thresholded F values.
* `Result.p_marked`: significant ANOVA p-value mask.
* `Result.index_PorN`: indices of selected positive or negative edges.
* `Result.num_sidePorN`: number of selected edges for each band.
* `Result.markedSide1`: first node of each selected edge.
* `Result.markedSide2`: second node of each selected edge.
* `Result.p_ttest`: post-hoc p values.
* `Result.t_ttest`: post-hoc t values.
* `Result.t_marked`: post-hoc t values for selected edges.

## `LOOCV_nested_feature_save.m`

### Input

The following folder paths need to be set before running:

* `inputpath_MCI_m1`: MCI group, movie 1 data folder.
* `inputpath_MCI_m2`: MCI group, movie 2 data folder.
* `inputpath_HC_m1`: HC group, movie 1 data folder.
* `inputpath_HC_m2`: HC group, movie 2 data folder.

Each input folder should contain `.mat` files. Each `.mat` file should contain:

* `EEG_results.M_zscore_mean`

The following parameters are also used:

* `outputpath`: folder for saving the output `.mat` file.
* `alpha`: statistical threshold.
* `fdrMode`: multiple-comparison correction mode.
* `useBands`: frequency bands used for feature extraction.

### Output

This script saves one `.mat` file:

* `LOOCV_nested_features_<fdrMode>.mat`

The saved file contains:

* `X_train_all`: training features for each LOOCV fold.
* `y_train_all`: training labels for each LOOCV fold.
* `X_test_all`: test features for each LOOCV fold.
* `y_test_all`: test labels for each LOOCV fold.
* `train_index_all`: training subject indices for each fold.
* `test_index_all`: test subject index for each fold.
* `edge_count_candidate`: number of candidate edges in each fold.
* `edge_count_negative`: number of negative edges in each fold.
* `edge_count_positive`: number of positive edges in each fold.
* `mask_candidate_all`: candidate edge masks.
* `mask_negative_all`: negative edge masks.
* `mask_positive_all`: positive edge masks.
* `feature_names`: names of the extracted features.
* `band_name`: names of the frequency bands.
* `N_HC`, `N_MCI`, `N`: subject numbers.
* `alpha`, `fdrMode`, `useBands`: analysis settings.
* `files_HC_m1`, `files_HC_m2`, `files_MCI_m1`, `files_MCI_m2`: input file lists.

## `svm_LOOCV_nested.py`

### Input

This script reads one MATLAB feature file:

* `LOOCV_nested_features_candidate_edges_fdr.mat`

The input `.mat` file should contain:

* `X_train_all`: training features for each LOOCV fold.
* `y_train_all`: training labels for each LOOCV fold.
* `X_test_all`: test features for each LOOCV fold.
* `y_test_all`: test labels for each LOOCV fold.
* `edge_count_negative`: number of negative edges in each fold.
* `edge_count_positive`: number of positive edges in each fold.
* `edge_count_candidate`: number of candidate edges in each fold.

### Output

This script creates the output folder:

* `SVM_results_candidate_edges_fdr`

The folder contains:

* `nested_LOOCV_predictions.xlsx`: true labels, predicted labels, and decision scores.
* `nested_LOOCV_performance.xlsx`: accuracy, balanced accuracy, sensitivity, specificity, AUC, and confusion matrix values.
* `nested_LOOCV_feature_weights.xlsx`: mean absolute SVM feature weights and standard errors.
* `nested_LOOCV_edge_count_summary.xlsx`: summary of candidate, negative, and positive edge counts.

The script also prints the SVM performance and AUC confidence interval in the console.

