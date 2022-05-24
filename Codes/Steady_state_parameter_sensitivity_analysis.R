# This code is written to compute the parameter sensitivity of models by changing the paramters by ni% 

# checking for necessary libraries
if(!("CoRC" %in% installed.packages()))
{
    # Install the CoRC package directly from GitHub:
    install.packages("remotes")
    remotes::install_github("jpahle/CoRC")
}
if(!("tibble" %in% installed.packages()))
{
    install.packages("tibble")
}
# loading necessary libraries
require(CoRC)
require(tibble)

#### steady state parameter sensitivity analysis ####
ni_ps_dir<-"/media/adithya03/data/HCC/NAD_project/TCGA_LIHC_ncount_models/parameter_sensitivity/steady_state_sensitivity/parameter_sensitivity_in_steady_state_with_1_%_change/"
#"/home/adithya03/HCC/NAD_project/TCGA_LIHC_patient_ncount_models/parameter_sensitivity/Steady_state_sensitivity/parameter_sensitivity_in_steady_state_with_5_%_change/"
# ni is the percentage value to change each parameter
ni <- 1 
patient_model <- commandArgs(trailingOnly = T)[1]
# patient_model <-"/home/adithya03/HCC/NAD_project/TCGA_LIHC_patient_ncount_models/patient_models/Normal_base_NAD_model.xml"
# print(patient_model)

model_name <- unlist(strsplit(basename(patient_model),split = "_"))[1]
print(model_name)
# load the model  
loadSBML(patient_model)
d_pars<-as.data.frame(getParameters()) # default parameters 
d_globalQ<-as.data.frame(getGlobalQuantities()) # default global quantities
# globalQR<-getGlobalQuantityReferences()
# pars_na<-which(is.na(pars$value))  # global quantaties in parameters
# View(pars[pars_na,])

# unaltered model steady state analysis ####
methods_1<-getSS(model = getCurrentModel())$method
methods_1$resolution<-1e-03 # setting the resolution to 1e-03
methods_1$maximum_duration_for_forward_integration<-1e20 # setting the maximum simulation duration to 1e20 
setSS(method = methods_1,model = getCurrentModel())
if (exists("result")) { rm(result)}
try(result <- runSteadyState(),silent = T) # running steady state analysis
if (exists("result"))
{
    # stable_state<-c(stable_state,p)
    if(result$result == "found")
    {
        species_names<-result[["species"]][,2] # species names
        reaction_names<-result[["reactions"]][,2] # reaction names
        patients_flux_data<-as.matrix(result[["reactions"]][,3]) # patient flux data
        patients_species_data<-as.matrix(result[["species"]][,4]) # patient species concentration
        patients_species_rate<-as.matrix(result[["species"]][,6]) # patient rate 
    }
    # model_status<-rbind(model_status,c(model_name,T))
    
    # final_conc1<-tail(result[["result"]],1)[,-1]  # extracting the final conc. of metabolites
    # creating a matrix to store final conc. values of vector parameter changes
    sensetivities_final_conc<-data.frame(matrix(NA,nrow = 2*nrow(d_pars),ncol = nrow(patients_species_data)+1))
    colnames(sensetivities_final_conc)<-c("Parameters",unlist(species_names))
    # creating a matrix to store species rate values of vector parameter changes
    sensetivities_final_rate<-data.frame(matrix(NA,nrow = 2*nrow(d_pars),ncol = nrow(patients_species_data)+1))
    colnames(sensetivities_final_rate)<-c("Parameters",unlist(species_names))
    # creating a matrix to store flux values of vector parameter changes
    sensetivities_final_flux<-data.frame(matrix(NA,nrow = 2*nrow(d_pars),ncol = nrow(patients_flux_data)+1))
    colnames(sensetivities_final_flux)<-c("Parameters",unlist(reaction_names))
    
    # sensetivities_final_conc[1,]<-c("Original",final_conc1)
    # changing the parameters with + n1 % and performing steady state analysis ####
    for (i in 1:nrow(d_pars))
    {
        pars<-d_pars
        globalQ<-d_globalQ
        if(is.na(pars[i,4]))
        {
            k= pars[i,5]
            n=match(k,globalQ[,1])
            globalQ[n,5]<-globalQ[n,5]*(1+(ni/100))
        } else
        {
            pars[i,4]<-pars[i,4]*(1+(ni/100))
        }
        pars<-as_tibble(pars)
        sens_par_name<-paste(pars[i,2],"+",ni,"%")
        setParameters(data = pars,model = getCurrentModel())
        if (exists("result")) { rm(result)}
        try(result <- runSteadyState(),silent = T)
        if (exists("result"))
        {
            final_conc<-as.matrix(result[["species"]][,4])
            final_rate<-as.matrix(result[["species"]][,6])
            final_flux<-as.matrix(result[["reactions"]][,3])
            sensetivities_final_conc[i,]<-c(sens_par_name,final_conc)
            sensetivities_final_rate[i,]<-c(sens_par_name,final_rate)
            sensetivities_final_flux[i,]<-c(sens_par_name,final_flux)
        } else
        {
            sensetivities_final_conc[i,1]<-c(sens_par_name)
            sensetivities_final_rate[i,1]<-c(sens_par_name)
            sensetivities_final_flux[i,1]<-c(sens_par_name)
        }
        
    }
    
    # changing the parameters with - n1 % and performing steady state analysis ####
    for (i in 1:nrow(d_pars))
    {
        pars<-d_pars
        globalQ<-d_globalQ
        if(is.na(pars[i,4]))
        {
            k= pars[i,5]
            n=match(k,globalQ[,1])
            globalQ[n,5]<-globalQ[n,5]*(1-(ni/100))
        } else
        {
            pars[i,4]<-pars[i,4]*(1-(ni/100))
        }
        pars<-as_tibble(pars)
        sens_par_name<-paste(pars[i,2],"-",ni,"%")
        setParameters(data = pars,model = getCurrentModel())
        if (exists("result")) { rm(result)}
        try(result <- runSteadyState(),silent = T)
        if (exists("result"))
        {
            final_conc<-as.matrix(result[["species"]][,4])
            final_rate<-as.matrix(result[["species"]][,6])
            final_flux<-as.matrix(result[["reactions"]][,3])
            sensetivities_final_conc[i+nrow(pars),]<-c(sens_par_name,final_conc)
            sensetivities_final_rate[i+nrow(pars),]<-c(sens_par_name,final_rate)
            sensetivities_final_flux[i+nrow(pars),]<-c(sens_par_name,final_flux)
        } else
        {
            sensetivities_final_conc[i+nrow(pars),1]<-c(sens_par_name)
            sensetivities_final_rate[i+nrow(pars),1]<-c(sens_par_name)
            sensetivities_final_flux[i+nrow(pars),1]<-c(sens_par_name)
        }
    }
    ### calculating percentage changes from unaltered model ####
    sensetivities_percentage_final_conc<-sensetivities_final_conc
    sensetivities_percentage_final_conc[,2:ncol(sensetivities_percentage_final_conc)]<-apply(sensetivities_percentage_final_conc[,2:ncol(sensetivities_percentage_final_conc)],2,as.numeric)
    sensetivities_percentage_final_rate<-sensetivities_final_rate
    sensetivities_percentage_final_rate[,2:ncol(sensetivities_percentage_final_rate)]<-apply(sensetivities_percentage_final_rate[,2:ncol(sensetivities_percentage_final_rate)],2,as.numeric)
    sensetivities_percentage_final_flux<-sensetivities_final_flux
    sensetivities_percentage_final_flux[,2:ncol(sensetivities_percentage_final_flux)]<-apply(sensetivities_percentage_final_flux[,2:ncol(sensetivities_percentage_final_flux)],2,as.numeric)
    # norm_data<-unlist(patients_species_data)
    # calculating percentage changes in concentrations
    for(j in 1:nrow(sensetivities_final_conc))
    {
        for (i in 2:ncol(sensetivities_final_conc))
        {
            sensetivities_percentage_final_conc[j,i]<-((sensetivities_percentage_final_conc[j,i]-patients_species_data[i-1,1])/patients_species_data[i-1,1])*100
        }
    }
    # calculating percentage changes in rates
    for(j in 1:nrow(sensetivities_final_rate))
    {
        for (i in 2:ncol(sensetivities_final_rate))
        {
            if(patients_species_rate[i-1,1] != 0 & sensetivities_percentage_final_rate[j,i] !=0)
            {
                sensetivities_percentage_final_rate[j,i]<-((sensetivities_percentage_final_rate[j,i]-patients_species_rate[i-1,1])/patients_species_rate[i-1,1])*100
            } else
            {
                sensetivities_percentage_final_rate[j,i]<-ifelse(patients_species_rate[i-1,1] == 0,
                                                                 ifelse(sensetivities_percentage_final_rate[j,i] ==0,1,
                                                                        sensetivities_percentage_final_rate[j,i]),
                                                                 patients_species_rate[i-1,1])
            }
        }
    }

    # calculating percentage changes in fluxes
    for(j in 1:nrow(sensetivities_final_flux))
    {
        for (i in 2:ncol(sensetivities_final_flux))
        {
            sensetivities_percentage_final_flux[j,i]<-((sensetivities_percentage_final_flux[j,i]-patients_flux_data[i-1,1])/patients_flux_data[i-1,1])*100
        }
    }
    ## listing the perturbed parameter with concentration
    pertubed_parameter_conc_data<-matrix(NA,ncol=3)
    cns<-colnames(sensetivities_percentage_final_conc)[-1]
    for(j in 1:nrow(sensetivities_percentage_final_conc))
    {
        per_conc<-which(abs(sensetivities_percentage_final_conc[j,-1])>ni)
        rn<-sensetivities_percentage_final_conc[j,1]
        
        if(length(per_conc)>1)
        {
            for (m in 1:length(per_conc))
            {
                temp3<-c(rn,cns[per_conc[m]],signif(sensetivities_percentage_final_conc[j,per_conc[m]+1],6))
                pertubed_parameter_conc_data<-rbind(pertubed_parameter_conc_data,temp3)
            }
        }
    }
    pertubed_parameter_conc_data<-na.omit(pertubed_parameter_conc_data)
    pertubed_parameter_conc_data1<-cbind(rep(model_name,nrow(pertubed_parameter_conc_data)),pertubed_parameter_conc_data)
    pertubed_parameter_conc_data1<-as.data.frame(pertubed_parameter_conc_data1)
    colnames(pertubed_parameter_conc_data1)<-c("Model","Parameter","Metabolite","Percentage_change")
    
    ## listing the perturbed parameter with flux
    pertubed_parameter_flux_data<-matrix(NA,ncol=3)
    flux_names<-colnames(sensetivities_percentage_final_flux)[-1]
    for(j in 1:nrow(sensetivities_percentage_final_flux))
    {
        per_flux<-which(abs(sensetivities_percentage_final_flux[j,-1])>ni)
        par_name<-sensetivities_percentage_final_flux[j,1]
        
        if(length(per_flux)>1)
        {
            for (m in 1:length(per_flux))
            {
                temp2<-c(par_name,flux_names[per_flux[m]],signif(sensetivities_percentage_final_flux[j,per_flux[m]+1],6))
                pertubed_parameter_flux_data<-rbind(pertubed_parameter_flux_data,temp2)
            }
        }
    }
    pertubed_parameter_flux_data<-na.omit(pertubed_parameter_flux_data)
    pertubed_parameter_flux_data1<-cbind(rep(model_name,nrow(pertubed_parameter_flux_data)),pertubed_parameter_flux_data)
    pertubed_parameter_flux_data1<-as.data.frame(pertubed_parameter_flux_data1)
    colnames(pertubed_parameter_flux_data1)<-c("Model","Parameter","Flux","Percentage_change")
    
    # pppc - patient parameter perturbed changes
    # for concentrations
    raw_c_filename<-paste(ni_ps_dir,model_name,"_model_",ni,"_percent_parameter_senstivity_at_ss_(conc_raw_data).txt",sep = "")
    per_c_filename<-paste(ni_ps_dir,model_name,"_model_",ni,"_percent_parameter_senstivity_at_ss_(conc_percentage_data).txt",sep = "")
    pertubed_conc_filename<-paste(ni_ps_dir,model_name,"_model_",ni,"_percent_parameter_senstivity_at_ss_(conc_perturbed_list).txt",sep = "")
    write.table(sensetivities_final_conc,raw_c_filename,sep = "\t",col.names = T,row.names = F)
    write.table(sensetivities_percentage_final_conc,per_c_filename,sep = "\t",col.names = T,row.names = F)
    write.table(pertubed_parameter_conc_data1,pertubed_conc_filename,sep = "\t",col.names = T,row.names = F,quote = F)
    
    # for Flux
    raw_f_filename<-paste(ni_ps_dir,model_name,"_model_",ni,"_percent_parameter_senstivity_at_ss_(flux_raw_data).txt",sep = "")
    per_f_filename<-paste(ni_ps_dir,model_name,"_model_",ni,"_percent_parameter_senstivity_at_ss_(flux_percentage_data).txt",sep = "")
    pertubed_flux_filename<-paste(ni_ps_dir,model_name,"_model_",ni,"_percent_parameter_senstivity_at_ss_(flux_perturbed_list).txt",sep = "")
    write.table(sensetivities_final_flux,raw_f_filename,sep = "\t",col.names = T,row.names = F)
    write.table(sensetivities_percentage_final_flux,per_f_filename,sep = "\t",col.names = T,row.names = F)
    write.table(pertubed_parameter_flux_data1,pertubed_flux_filename,sep = "\t",col.names = T,row.names = F,quote = F)
    
    # for rate
    raw_r_filename<-paste(ni_ps_dir,model_name,"_model_",ni,"_percent_parameter_senstivity_at_ss_(rate_raw_data).txt",sep = "")
    per_r_filename<-paste(ni_ps_dir,model_name,"_model_",ni,"_percent_parameter_senstivity_at_ss_(rate_percentage_data).txt",sep = "")
    write.table(sensetivities_final_rate,raw_r_filename,sep = "\t",col.names = T,row.names = F)
    write.table(sensetivities_percentage_final_rate,per_r_filename,sep = "\t",col.names = T,row.names = F)
    
}
res<-paste(model_name," sensitivity analysis is completed.")
print(res)
