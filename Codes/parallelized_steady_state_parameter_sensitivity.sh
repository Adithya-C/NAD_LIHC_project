# directory containing patient models
data_dir="/home/NAD_project/TCGA_LIHC_ncount_models/TCGA_LIHC_patient_models/"

# command to find patine tmodels and run steady state sensitivity for each model
find "$data_dir" -name "*.xml" | xargs -n1 -P6 Rscript /home/NAD_project/Rscripts/Steady_state_parameter_sensitivity_analysis.R
