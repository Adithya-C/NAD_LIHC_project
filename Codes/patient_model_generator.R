# this function is written to generate tumour specific models 
# input data for this function is patient F_Factor data.frame, base model and output directory
# description of f_factor data frame : f_factors as rows and patients as columns, numeric data entries
patient_model_generator<-function(f_factor_data,base_model_file,model_dir)
{
    # checking necessary libraries
    required_libraries<-c("stringi","readr")
    check_libraries<-(required_libraries %in% installed.packages())
    libraries_to_install<-required_libraries[!check_libraries]
    if (length(libraries_to_install)>1)
    {
        install.packages(libraries_to_install)
    }
    # Loading necessary libraries
    for( p in required_libraries) require(p,character.only = T) 
    
    #model_dir<-"/home/NAD_project/TCGA_LIHC_ncount_models/TCGA_LIHC_patient_models/"
    
    # reading the base model 
    #base_model_file<-"/home/NAD_project/Normal_base_model.xml"
    base_model<-read_lines(base_model_file)
    
    f_factors<-row.names(f_factor_data)
    patient_ids<-colnames(f_factor_data)
    
    f_factor_rex<-paste("\"",f_factors,"\"",sep = "")
    f_factor_rows<-c()
    for (i in 1:length(f_factor_rex))
    {
        f_factor_rows<-c(f_factor_rows,grep(f_factor_rex[i],base_model))
    }
    
    for (p in 1:ncol(f_factor_data))
    {
        patient_model<-base_model
        for (f in 1:nrow(f_factor_data))
        {
            patient_model[f_factor_rows[f]]<-stri_replace_last_regex(patient_model[f_factor_rows[f]],"1",f_factor_data[f,p])
        }
        patient_file<-paste(model_dir,patient_ids[p],"_tumour_model.xml",sep = "")
        write_lines(patient_model,patient_file)
    }
}