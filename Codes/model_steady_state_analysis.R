# This code is written to calculate the steady states for patient models
# checking for necessary libraries

if(!("CoRC" %in% installed.packages()))
{
    # Install the CoRC package directly from GitHub:
    install.packages("remotes")
    remotes::install_github("jpahle/CoRC")
}
if(!("ggplot2" %in% installed.packages()))
{
    install.packages("ggplot2")
}
# loading necessary libraries
require(CoRC)
require(ggplot2)

# input files and folders
# directory to store the data
data_dir<-"/home/NAD_project/TCGA_LIHC_ncount_models/"
# directory of patient models
patient_models_dir<-"/home/NAD_project/TCGA_LIHC_ncount_models/TCGA_LIHC_patient_models/"
# directory to store the plots
plots_dir<-"/home/NAD_project/TCGA_LIHC_ncount_models/plots/"

####  COPASI steady state analysis for time histogram ####
setwd(patient_models_dir)
patient_ids<-c() # vector to save patient ids
patient_models<- list.files(path = patient_models_dir, pattern = '.xml', full.names = T)
unstable_states<-c() # vector to save unstable models
stable_states<-c() # vector to save stable models
steady_state_model_times<-data.frame(matrix(NA,ncol = 2)) # data frame to store stable state times
dur_t<-10^seq(0,20,by=1) # range of duration to check for stablitiy
# running patient specific models for stability
for (p in 1:length(patient_models))
{
    model_name <- unlist(strsplit(basename(patient_models[p]),split = "_"))[1]
    patient_ids<-c(patient_ids,model_name)
    cat("\n",model_name)
    # load SBML model  
    loadSBML(patient_models[p])
    # get steady state settings
    methods_1<-getSS(model = getCurrentModel())$method
    methods_1$resolution<-1e-03 # resolution
    methods_1$iteration_limit<-100 # iteration limit
    for(j in 1:length(dur_t))
    {
        methods_1$maximum_duration_for_forward_integration <- dur_t[j] # updating the forward integration time
        setSS(method = methods_1,model = getCurrentModel()) # setting the integration time
        if (exists("result")) { rm(result)}
        try(result <- runSteadyState(),silent = T) # calculating the model for steady states
        # checking the model for steady state
        if (exists("result")) 
        {
            stable_states<-c(stable_states,p) # updating the vector of stable models
            sst<-c(model_name,dur_t[j])
            cat("\t",dur_t[j])
            # updating the dataframe of steady state times
            steady_state_model_times<-rbind(steady_state_model_times,sst) 
            break
        } 
        # if the model did not attain steady state by 1E20 and resolution 1E-03, its a unstable model 
        if (j==length(dur_t))
        {
            if(!exists("result"))
            {
                cat("\t",dur_t[j]," unstable")
                unstable_states<-c(unstable_states,p) # updating the list of unstable models
            }
        }
    }
}
# plotting steady state times as histogram
steady_state_model_times<-na.omit(steady_state_model_times)
# converting the time into log scale
steady_state_model_times$logtimes<-log10(as.numeric(steady_state_model_times$X2)) 
h<-hist(steady_state_model_times$logtimes,breaks = seq(0,20),plot = F)
y_max=round(max(h$counts)/100)*100

# generating the histogram
t_hist<-ggplot(steady_state_model_times,aes(x=logtimes))+geom_histogram(binwidth = 1,col="grey",fill="pink")+
    stat_bin(aes(y=..count..,label=ifelse(..count..==0,"",..count..)),geom = "text",vjust=-1,hjust=0,bins = 21)+
    labs(x="Times (log10)",title = paste("Stable state times,n = ",nrow(steady_state_model_times),",e = ",methods_1$resolution,sep = ""))+
    ylim(0,y_max)+xlim(0,20)+
    theme(axis.line = element_line(linetype = 1,color = "black"),
          axis.title = element_text(family = "Arial",face = "bold",size = 14,color="black"),
          axis.text = element_text(family = "Arial",face = "bold",size = 12,color="black"),
          axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
          panel.grid = element_blank(),panel.background = element_blank())

# saving the histogram
hist_name<-paste(plots_dir,"steady_state_time_histogram_(e=",methods_1$resolution,").png",sep="")
png(hist_name,width = 8,height = 6,units = "in",type = "cairo",res = 600)
par(xpd=T,mar=c(5,6,4,4))
plot(t_hist)
dev.off()

# print no. of stable and unstable models
print(paste("No. of stable models : ", length(stable_states),sep = ""))
print(paste("No. of unstable models : ", length(unstable_states),sep = ""))

# writing steady state times data into files
stable_state_times_name<-paste(data_dir,"stable_model_times_(e=",methods_1$resolution,").txt",sep="")
write.table(steady_state_model_times,stable_state_times_name,sep = "\t",col.names = F,row.names = F)

stable_models<-patient_ids[stable_states]
unstable_models<-patient_ids[unstable_states]

# writing stable and unstable models into files
stable_model_files<-paste(data_dir,"stable_models_(e=",methods_1$resolution,").txt",sep="")
unstable_model_files<-paste(data_dir,"unstable_models_(e=",methods_1$resolution,").txt",sep="")
write.table(stable_models,stable_model_files,sep = "\t",col.names = F,row.names = F)
write.table(unstable_models,unstable_model_files,sep = "\t",col.names = F,row.names = F)

#### Model stability analysis ####
# steady state analysis for stable and unstable model classification
setwd(patient_models_dir)
patient_ids<-c()
# list of patient models
patient_models<- list.files(path = patient_models_dir, pattern = '.xml', full.names = T)
# data frame to store metabolite concentrations
patients_species_data<-matrix(NA,ncol = length(patient_models),nrow = 14) 
# data frame to store rate of conversion of species
patients_species_rate<-matrix(NA,ncol = length(patient_models),nrow = 14)
# data frame to store reaction flux data
patients_flux_data<-matrix(NA,ncol = length(patient_models),nrow = 26)
# vectors to store stable and unstable models
unstable_state<-c()
stable_state<-c()
model_status<-data.frame(matrix(NA,ncol = 3))

# checkng for model steady state for patient models
for(p in 1:length(patient_models))
{
    model_name <- unlist(strsplit(basename(patient_models[p]),split = "_"))[1]
    # patient_ids<-c(patient_ids,model_name)
    cat("\n",model_name,"\t")
    #load SBML model  
    loadSBML(patient_models[p])
    #  get steady state settings
    methods_1<-getSS(model = getCurrentModel())$method
    # updating the resolution,iteration limit, integration times
    methods_1$resolution<-1e-03
    methods_1$iteration_limit<-100
    methods_1$maximum_duration_for_forward_integration<-1e20
    methods_1$maximum_duration_for_backward_integration<-1e20
    # setting the methods for stady state function
    setSS(method = methods_1,model = getCurrentModel())

    if (exists("result")) { rm(result)}
    try(result <- runSteadyState(),silent = T) # simulating the model for steady state
    if (exists("result"))
    {
        steadystate_result<-result$result;
        f_res<-grep("found",steadystate_result)
        if(length(f_res))
        {
            ssr<-T
        }
    }
    else {ssr<-F;}#break}
    # }
    if(ssr&length(stable_state)<2)
    {
        species_names<-result[["species"]][,2]
        reaction_names<-result[["reactions"]][,2]
    }
    if(ssr)
    {
        # updating the stable models and storing the steady state details
        stable_state<-c(stable_state,model_name)
        patients_flux_data[,p]<-as.matrix(result[["reactions"]][,3])
        patients_species_data[,p]<-as.matrix(result[["species"]][,4])
        patients_species_rate[,p]<-as.matrix(result[["species"]][,6])
        model_status<-rbind(model_status,c(model_name,T,steadystate_result))
    } else
    {
        # updating the unstable models
        unstable_state<-c(unstable_state,model_name)
        model_status<-rbind(model_status,c(model_name,F,"Notfound"))
    }
    if(nrow(model_status)==2)
    {
        model_status<-na.omit(model_status)
    }
    cat(model_status[p,3])
}
model_status<-na.omit(model_status)
# model_status<-model_status[-1,] # removing the normal base model
colnames(model_status)<-c("patient_id","status","result")
rownames(model_status)<-c(1:length(patient_models))
patients_species_data1<-t(na.omit(t(patients_species_data)))
colnames(patients_species_data1)<-stable_state
row.names(patients_species_data1)<-unlist(species_names)
patients_species_rate1<-t(na.omit(t(patients_species_rate)))
colnames(patients_species_rate1)<-stable_state
row.names(patients_species_rate1)<-unlist(species_names)
patients_flux_data1<-t(na.omit(t(patients_flux_data)))
colnames(patients_flux_data1)<-stable_state
row.names(patients_flux_data1)<-unlist(reaction_names)

print(table(model_status$status))
print(table(model_status$result))

# writing the stable model steady state information into files
patients_species_data_file<-paste(data_dir,"TCGA_LIHC_patient_metabolite_stable_state_concentration_(e=",methods_1$resolution,").txt",sep="")
write.table(patients_species_data1,patients_species_data_file,sep = "\t",col.names = T,row.names = T)
patients_species_rate_file<-paste(data_dir,"TCGA_LIHC_patient_metabolite_stable_state_species_rate_(e=",methods_1$resolution,").txt",sep="")
write.table(patients_species_rate1,patients_species_rate_file,sep = "\t",col.names = T,row.names = T)
patients_flux_data_file<-paste(data_dir,"TCGA_LIHC_patient_metabolite_stable_state_fluxes_(e=",methods_1$resolution,").txt",sep="")
write.table(patients_flux_data1,patients_flux_data_file,sep = "\t",col.names = T,row.names = T)
model_status_file<-paste(data_dir,"patient_model_stability_status_(e=",methods_1$resolution,").txt",sep="")
write.table(model_status[-1,],model_status_file,sep = "\t",col.names = T,row.names = F)

stable_models_file<-paste(data_dir,"patient_stable_models_(e=",methods_1$resolution,").txt",sep="")
unstable_models_file<-paste(data_dir,"patient_unstable_models_(e=",methods_1$resolution,").txt",sep="")
write.table(stable_state[-1],stable_models_file,sep = "\t",col.names = F,row.names = F)
write.table(unstable_state,unstable_models_file,sep = "\t",col.names = F,row.names = F)

# printing the no. of stable and unstable models
print(paste("No. of stable models : ", length(stable_state),sep = ""))
print(paste("No. of unstable models : ", length(unstable_state),sep = ""))
