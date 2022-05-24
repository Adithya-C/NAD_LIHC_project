# this code is written to simulate the time course for patient models
# checking for necessary libraries
if(!("CoRC" %in% installed.packages()))
{
    # Install the CoRC package directly from GitHub:
    install.packages("remotes")
    remotes::install_github("jpahle/CoRC")
}
# loading necessary libraries
require(CoRC)

# input patient model directory
patient_models_dir<-"/home/NAD_project/TCGA_LIHC_ncount_models/TCGA_LIHC_patient_models/"

# output time series directory
time_series_dir<-"/home/NAD_project/TCGA_LIHC_ncount_models/time_course/"

# outpur metabolite concentrations directory
met_dir<-"/home/NAD_project/TCGA_LIHC_ncount_models/time_course/metabolites_tc_data/"

####  COPASI Time course analysis ####
setwd(patient_models_dir)
patient_ids<-c()
patient_models<- list.files(path = patient_models_dir, pattern = '.xml', full.names = T)
# simulating each patient model for 1E16 time
for (p in 1:length(patient_models))
{
    model_name <- unlist(strsplit(basename(patient_models[p]),split = "_"))[1]
    patient_ids<-c(patient_ids,model_name)
    print(model_name)
    #  load the model  
    loadSBML(patient_models[p])
    # change the time course setting and simulate the model
    try(setTC(duration = 1e04,dt=1e02)) 
    resTC<-runTC() 
    temp1<-resTC[["result"]][-1,]
    
    try(setTC(duration = 1e08,dt=1e04)) 
    resTC<-runTC() 
    temp2<-resTC[["result"]][-1,]
    
    try(setTC(duration = 1e012,dt=1e08)) 
    resTC<-runTC()
    temp3<-resTC[["result"]][-1,]
    
    try(setTC(duration = 1e016,dt=1e12)) 
    resTC<-runTC()
    temp4<-resTC[["result"]][-1,]
    
    # combining all the time course data
    tc<-rbind(temp1,temp2,temp3,temp4)
    
    # writting the time corse simulation data into file
    ts_name<-paste(time_series_dir,model_name,".txt",sep = "")
    write.table(tc,ts_name,sep = "\t",col.names = T,row.names = F)

    # selecting the time points to store the data
    tp<-10^(2:16)
    tp_rows<-match(tp,tc$Time)
    
    # subsetting data and converting the data into log10 scale
    tp_conc<-tc[tp_rows,]
    tp_conc_log<-apply(tp_conc,2,log10)
    final_conc<-tail(tc,1)[,-1] 
    # storing the species concentrations into data frame
    if(length(patient_ids)<2)
    {
        patient_species_conc<-matrix(0,ncol = length(patient_models),nrow = length(final_conc))
        rownames(patient_species_conc)<-names(final_conc)
        metabolite_tc<-paste(names(final_conc),"_tc",sep="")
        metabolite_log_tc<-paste(names(final_conc),"_log_tc",sep="")
        for(i in 1:length(metabolite_tc))
        {
            temp<-tp_conc[,c(1,i+1)]
            colnames(temp)[2]<-model_name
            assign(metabolite_tc[i],temp)
            temp1<-tp_conc_log[,c(1,i+1)]
            colnames(temp1)[2]<-model_name
            assign(metabolite_log_tc[i],temp1)
        }
    }
    else
    {
        for(i in 1:length(metabolite_tc))
        {
            temp<-tp_conc[,c(1,i+1)]
            colnames(temp)[2]<-model_name
            temp_1<-merge(get(metabolite_tc[i]),temp,by=1)
            assign(metabolite_tc[i],temp_1)
            temp1<-tp_conc_log[,c(1,i+1)]
            colnames(temp1)[2]<-model_name
            temp_2<-merge(get(metabolite_log_tc[i]),temp1,by=1)
            assign(metabolite_log_tc[i],temp_2)
        }
    }
    patient_species_conc[,p]<-t(final_conc)
    # removing recurring variables
    rm(temp1,temp2,temp3,temp4,tc,resTC,final_conc,tp_conc,tp_conc_log)
}
colnames(patient_species_conc)<-patient_ids

# writing patient species data frame into files 
patient_species_conc_name<-paste(met_dir,"TCGA_LIHC_metabolite_final_time_series_data.txt",sep = "")
write.table(patient_species_conc,patient_species_conc_name,sep = "\t",col.names = T,row.names = T)

for(i in 1:length(metabolite_tc))
{
    met_filename<-paste(met_dir,"TCGA_LIHC_patients_",metabolite_tc[i],".txt",sep = "")
    met_log_filename<-paste(met_dir,"TCGA_LIHC_patients_",metabolite_tc[i],"_(log10).txt",sep = "")
    write.table(get(metabolite_tc[i]),met_filename,col.names = T,row.names = T,sep = "\t")
    write.table(get(metabolite_log_tc[i]),met_log_filename,col.names = T,row.names = T,sep = "\t")
}
# End