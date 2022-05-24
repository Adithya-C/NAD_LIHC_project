# This code is written to perform survival analysis on NADnet genes and gene combinations
# checking for necessary package availability 
required_packages<-c("survival","survminer","RColorBrewer","ggplot2","ggfortify","corrplot")
check_packages<-(required_packages %in% installed.packages())
packages_to_install<-required_packages[!check_packages]
# Installing package
if(length(packages_to_install)>1)
{
    install.packages(packages_to_install)
}
# Loading necessary libraries
for( p in required_packages) require(p,character.only = T) 

# setting colors for models/groups
color_schema<-brewer.pal(n = 9, name = 'Spectral')
survival_colors<-color_schema[c(1,9)] # colour schema for survival groups

# necessary input files
survival_data_file<-"/home/NAD_project/TCGA_HCC/TCGA-LIHC.survival.tsv"
survival_data<-read.delim(survival_data_file)
plots_dir<-"/home/NAD_project/TCGA_LIHC_ncount_models/plots/"
# plots_dir<-"/media/adithya03/data/HCC/NAD_project/TCGA_LIHC_ncount_models/Final_plots/";
data_dir<-"/home/NAD_project/TCGA_LIHC_ncount_models/"

model_metabolites_file<-"/home/NAD_project/NAD_model_metabolites.csv"
model_metabolites<-read.delim(model_metabolites_file,header = F,quote = "")
# model status file has two colums, one woth patinet ids and second with steady state status
model_status_file<-"/home/NAD_project/TCGA_LIHC_ncount_models/patient_model_stability_status_(e=0.001).txt"
model_status<-read.delim(model_status_file,header = T)

# phenotypic data selection ####
# phenotype file has phenotype data provided by TCGA-LIHC
phenotype_file<-"/home/NAD_project/TCGA_HCC/TCGA-LIHC.GDC_phenotype.tsv"
phenotype_data<-read.delim(phenotype_file,header = T)
phenotypic_data<-phenotype_data[match(model_status[,1],phenotype_data[,1]),]
phenotypic_data_final<-phenotypic_data[,c(1,3,2,15,25,37:41,75:77,95,99:101,73,83)]
row.names(phenotypic_data_final)<-phenotypic_data_final[,1]
phenotypic_data_final[,3:14]<-apply(phenotypic_data_final[,3:14],2,as.factor)
phenotypic_data_final$OS.time<-apply(phenotypic_data_final[,18:19],1,max,na.rm=T) 
phenotypic_data_final$OS<-ifelse(phenotypic_data_final$vital_status.demographic=="Alive",0,1)

phenotypic_data_final$Stage<-phenotypic_data_final$tumor_stage.diagnoses
phenotypic_data_final$Stage[grep("not reported",phenotypic_data_final$Stage)]<-NA
phenotypic_data_final$Stage<-gsub("stage ","",phenotypic_data_final$Stage)
# clanging the staging information to major stages
phenotypic_data_final$Stage[grep("^i$",phenotypic_data_final$Stage)]<-paste("S",as.roman(1),sep = "_")
phenotypic_data_final$Stage[grep("^ii$",phenotypic_data_final$Stage)]<-paste("S",as.roman(2),sep = "_")
phenotypic_data_final$Stage[grep("^iii",phenotypic_data_final$Stage)]<-paste("S",as.roman(3),sep = "_")
phenotypic_data_final$Stage[grep("^iv",phenotypic_data_final$Stage)]<-paste("S",as.roman(4),sep = "_")

# patient model data file
logfc_data_file<-"/home/NAD_project/TCGA_LIHC_ncount_models/TCGA-LIHC_tumour_NAD_logFC_data(edger_normalized_counts_by_mean).tsv"
metabolite_status_file<-"/home/NAD_project/TCGA_LIHC_ncount_models/TCGA_LIHC_Stable_model_metabolite_groups_(e=0.001).txt"
metabolite_lfc_file<-"/home/NAD_project/TCGA_LIHC_ncount_models/TCGA_LIHC_Stable_model_log2fc_metabolites_(e=0.001).txt"

logfc_data<-read.delim(logfc_data_file)
metabolite_status<-read.delim(metabolite_status_file)
metabolite_lfc<-read.delim(metabolite_lfc_file)

colnames(logfc_data)<-gsub("\\.","-",colnames(logfc_data))
colnames(metabolite_lfc)<-gsub("\\.","-",colnames(metabolite_lfc))
row.names(metabolite_status)<-metabolite_status[,1]
metabolite_status<-metabolite_status[,-1]

# classifying log2fc data
logfc_data1<-apply(logfc_data,2,function(x) ifelse( is.finite(x),ifelse(x>(1),"Up",ifelse(x<(-1),"Down","No_change")),NA))
logfc_data1<-apply(logfc_data1,2,factor,levels=c("No_change","Down","Up"))
logfc_data1<-as.data.frame(logfc_data1)
colnames(logfc_data1)<-colnames(logfc_data)

survival_logfc_ids<-na.omit(match(colnames(logfc_data),survival_data$sample))
logfc_survival_ids<-survival_data$sample[survival_logfc_ids]

survival_met_ids<-na.omit(match(rownames(metabolite_status),survival_data$sample))
met_survival_ids<-survival_data$sample[survival_met_ids]

gene_met_lfc_ids<-na.omit(match(colnames(metabolite_lfc),colnames(logfc_data)))
gene_met_lfc<-logfc_data[,gene_met_lfc_ids]
gene_met_lfc_status<-as.data.frame(t(logfc_data1[,gene_met_lfc_ids]))

survival_logfc_data<-cbind(survival_data[survival_logfc_ids,c(2,4)],t(logfc_data1[,logfc_survival_ids]))
survival_met_data<-cbind(survival_data[survival_met_ids,c(2,4)],metabolite_status[met_survival_ids,])

survival_logfc_data<-data.frame(apply(survival_logfc_data,2,function(x) ifelse(is.infinite(x),NA,x)))

survival_logfc_data[,2]<-as.numeric(survival_logfc_data[,2])

row.names(survival_logfc_data)<-logfc_survival_ids
row.names(survival_met_data)<-met_survival_ids


## metabolite HR analysis independent of stage data ####
metabolite_hr_analysis<-as.data.frame(matrix(NA,ncol=19))
colnames(metabolite_hr_analysis)<-c("Metabolite","Group_1","Group_2","N_total","N_event","N_group1","N_group2","Obs(group1)","Obs(group2)","Coef","exp(coef)[HR]","SE(coef)","Z","p.value(Z)","exp(-coef)","CI(95_lower)","CI(95_upper)","Chisq","p.value(Chisq)")

for(j in 1:ncol(metabolite_status))
{
    temp_met<-survival_met_data[,c(1,2,j+2)]
    #temp_met<-temp_met[which(temp_met[,3]!= "No_change"),]
    met_name<-colnames(temp_met)[3]
    temp_met[,3]<-as.factor(temp_met[,3])
    if (length(unique(temp_met[,3]))>1)
    {
        m12<-combn(levels(temp_met[,3]),2)
        for (m in 1:ncol(m12))
        {
            temp_m1<-temp_met[grep(m12[1,m],temp_met[,3]),]
            temp_m2<-temp_met[grep(m12[2,m],temp_met[,3]),]
            temp_met12<-rbind(temp_m1,temp_m2)
            temp_met12[,3]<-factor(temp_met12[,3])
            cox_met_function<-coxph(Surv(temp_met12$OS.time,temp_met12$OS == 1)~., data = temp_met12)
            s<-summary(cox_met_function)
            colnames(temp_met12)[3]<-"met"
            group_name<-paste("Metabolite :",met_name,sep = " ")
            surv_diff <- survdiff(Surv(OS.time, OS) ~ met, data=temp_met12)
            # print(met_name)
            # print(surv_diff)
            pvalue<-( 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1));
            met_comb_name<-paste(met_name,"(",paste(m12[,m],collapse = ","),")",sep = "");
            print(paste(met_comb_name,", HR :",signif(s$coefficients[2],3),", p.value :",signif(s$coefficients[5],3)));
            
            # # making the KM plot
            # met_km_plot_title<-paste(group_name,", HR :",signif(s$coefficients[2],3),", p.value :",signif(s$coefficients[5],3))
            # met_plt<-ggsurvplot(
            #     fit = survfit(Surv(OS.time, OS) ~ met, data=temp_met12),
            #     break.time.by=500,
            #     xlab = "Days",ylab = "Overall survival probability",
            #     risk.table = T,title =  met_km_plot_title)
            # 
            # # saving the plot to file PNG format
            # sur_met_plot_file<-paste(plots_dir,met_comb_name,"_KM_plot.png",sep = "")
            # png(sur_met_plot_file,width = 8,height = 6,units = "in",type = "cairo",res = 600)
            # par(xpd=T,mar=c(5,6,4,4))
            # print(met_plt)
            # dev.off()
            # 
            # # saving the plot to file SVG format
            # sur_met_plot_file_svg<-paste(plots_dir,met_comb_name,"_KM_plot.svg",sep = "")
            # svg(sur_met_plot_file_svg,width = 8,height = 6)
            # par(xpd=T,mar=c(5,6,4,4))
            # print(met_plt)
            # dev.off()
            
            # updating the hetabolite HR file 
            temp_hr<-c(met_name,m12[1,m],m12[2,m],s$n,s$nevent,surv_diff$n,surv_diff$obs,signif(s$coefficients,6),signif(s$conf.int[-1],6),signif(surv_diff$chisq,6),signif(pvalue,6))
            metabolite_hr_analysis<-rbind(metabolite_hr_analysis,temp_hr)
        }
    }
}

metabolite_hr_analysis<-na.omit(metabolite_hr_analysis)
metabolite_hr_analysis_signif<-metabolite_hr_analysis[which(metabolite_hr_analysis$`p.value(Z)`<0.05),]

metabolite_hr_analysis_file<-paste(data_dir,"TCGA_LIHC_Metabolite_group_HR_analysis.tsv",sep = "")
write.table(metabolite_hr_analysis,metabolite_hr_analysis_file,sep = "\t",col.names = T,row.names = F)

metabolite_hr_analysis_file_signif<-paste(data_dir,"TCGA_LIHC_Metabolite_group_HR_analysis_(significant).tsv",sep = "")
write.table(metabolite_hr_analysis_signif,metabolite_hr_analysis_file_signif,sep = "\t",col.names = T,row.names = F)

# 
## NAD gene group combinations survival analysis independent of stage ####
# code editted on 02 May 2022
nad_gene_data<-cbind(phenotypic_data_final[logfc_survival_ids,c(21,20)],survival_logfc_data[,-c(1:2)])
gene_hr_analysis<-as.data.frame(matrix(NA,ncol=19))
colnames(gene_hr_analysis)<-c("Gene_comb","Group_1","Group_2","N_total","N_event","N_group1","N_group2","Obs(group1)","Obs(group2)","Coef","exp(coef)","SE(coef)","Z","p.value(Z)","exp(-coef)","CI(95_lower)","CI(95_upper)","Chisq","p.value(Chisq)")
gene_comb_hr_analysis<-gene_hr_analysis

for(g in 3:31)
{
    temp_gene<-na.omit(nad_gene_data[,c(1:2,g)])
    gene_name<-colnames(temp_gene)[3]
    temp_gene[,3]<-as.factor(temp_gene[,3])
    #temp_gene<-temp_gene[which(temp_gene[,3]!= "No_change"),]
    temp_gene2<-temp_gene
    if (length(unique(temp_gene2[,3]))>1)
    {
        g12<-combn(levels(temp_gene2[,3]),2)
        for (j in 1:ncol(g12))
        {
            temp_g1<-temp_gene2[grep(g12[1,j],temp_gene2[,3]),]
            temp_g2<-temp_gene2[grep(g12[2,j],temp_gene2[,3]),]
            comb_g<-rbind(temp_g1,temp_g2)
            comb_g[,3]<-factor(comb_g[,3])
            cox_comb_gene_function<-coxph(Surv(comb_g$OS.time,comb_g$OS == 1)~ ., data =comb_g)
            s_comb_gene<-summary(cox_comb_gene_function)
            
            gene_comb_name<-paste(gene_name,"(",paste(g12[,j],collapse = ","),")",sep = "")
            colnames(comb_g)[3]<-"gene"
            gene_comb_surv_diff <- survdiff(Surv(OS.time, OS) ~ gene, data=comb_g)
            gene_comb_pvalue<-( 1 - pchisq(gene_comb_surv_diff$chisq, length(gene_comb_surv_diff$n) - 1))
            gene_comb_km_plot_title<-paste(gene_comb_name,", HR :",signif(s_comb_gene$coefficients[2],3),", p.value :",signif(s_comb_gene$coefficients[5],3))
            print(paste(gene_comb_name,", HR :",signif(s_comb_gene$coefficients[2],3),", p.value :",signif(s_comb_gene$coefficients[5],3)))
            
            # # generating KM plots
            # gene_plt<-ggsurvplot(
            #     fit = survfit(Surv(OS.time, OS) ~ gene, data=comb_g),
            #     break.time.by=500,palette = survival_colors,#pval = T,
            #     xlab = "Days",ylab = "Overall survival probability",
            #     legend.title= "Gene status",legend.labs=g12[,j],
            #     risk.table = T,tables.height = 0.3,
            #     font.legends=16,font.tickslab=16,
            #     title = gene_comb_km_plot_title)
            # 
            # # saving the plot into PNG file 
            # gene_sur_plot_file<-paste(plots_dir,gene_comb_name,"_Group_KM_plot.png",sep = "")
            # png(gene_sur_plot_file,width = 8,height = 6,units = "in",type = "cairo",res = 600)
            # par(xpd=T,mar=c(5,6,4,4))
            # print(gene_plt)
            # dev.off()
            # 
            # # saving the plot into SVG file 
            # gene_sur_plot_file_svg<-paste(plots_dir,gene_comb_name,"_Group_KM_plot.svg",sep = "")
            # svg(gene_sur_plot_file_svg,width = 8,height = 6)
            # par(xpd=T,mar=c(5,6,4,4))
            # print(gene_plt)
            # dev.off()
            # 
            # updating the gene comb hr data
            temp_g<-c(gene_names,g12[,j],s_comb_gene$n,s_comb_gene$nevent,gene_comb_surv_diff$n,gene_comb_surv_diff$obs,signif(s_comb_gene$coefficients,6),signif(s_comb_gene$conf.int[-1],6),signif(gene_comb_surv_diff$chisq,6),signif(gene_comb_pvalue,6),signif(s_comb_gene$waldtest[3],6))
            gene_hr_analysis<-rbind(gene_hr_analysis,temp_g)
        }
    }
} 
gene_hr_analysis<-na.omit(gene_hr_analysis)
gene_hr_signif<-gene_hr_analysis[which(gene_hr_analysis$`p.value(Z)`<0.05),]

gene_hr_analysis_file<-paste(data_dir,"TCGA_LIHC_Gene_group_HR_analysis.tsv",sep = "")
write.table(gene_hr_analysis,gene_hr_analysis_file,sep = "\t",col.names = T,row.names = F)

gene_hr_signif_file<-paste(data_dir,"TCGA_LIHC_Gene_group_HR_analysis_(significant).tsv",sep = "")
write.table(gene_hr_signif,gene_hr_signif_file,sep = "\t",col.names = T,row.names = F)

# gene combination HR analysis
for(g1 in 3:31)
{
    for (g2 in 3:31)
    {
        if (!(g1==g2))
        {
            temp_gene<-na.omit(nad_gene_data[,c(1:2,g1,g2)])
            # temp_gene<-na.omit(survival_gene_met_status_data[,c(1:2,g1,g2)])
            
            temp_gene$gene_comb<-apply(temp_gene[,3:4],1, function(x) paste(x,collapse = "_"))
            gene_names<-paste(colnames(temp_gene)[3:4],collapse = "_")
            group_name<-paste("Gene :",gene_names,sep = " ")
            temp_gene$gene_comb<-as.factor(temp_gene$gene_comb)
            # rm_samples<-which(temp_gene$gene_comb=="nc_nc")
            # if(length(rm_samples)>1)
            # {
            #     temp_gene1<-temp_gene[-which(temp_gene$gene_comb=="nc_nc"),]
            # } else { temp_gene1<-temp_gene }
            temp_gene2<-temp_gene[,c(1,2,5)]
            temp_gene2[,3]<-factor(temp_gene2[,3])
            if (length(unique(temp_gene2[,3]))>1)
            {
                g12<-combn(levels(temp_gene2$gene_comb),2)
                for (j in 1:ncol(g12))
                {
                    temp_g1<-temp_gene2[grep(g12[1,j],temp_gene2$gene_comb),]
                    temp_g2<-temp_gene2[grep(g12[2,j],temp_gene2$gene_comb),]
                    comb_g<-rbind(temp_g1,temp_g2)
                    comb_g[,3]<-factor(comb_g[,3])
                    cox_comb_gene_function<-coxph(Surv(comb_g$OS.time,comb_g$OS == 1)~ ., data =comb_g)
                    s_comb_gene<-summary(cox_comb_gene_function)
                    
                    gene_comb_name<-paste(gene_names,"(",paste(g12[,j],collapse = ","),")",sep = "")
                    
                    gene_comb_surv_diff <- survdiff(Surv(OS.time, OS) ~ gene_comb, data=comb_g)
                    gene_comb_pvalue<-( 1 - pchisq(gene_comb_surv_diff$chisq, length(gene_comb_surv_diff$n) - 1))
                    gene_comb_km_plot_title<-paste(group_name,", HR :",signif(s_comb_gene$coefficients[2],3),", p.value :",signif(s_comb_gene$coefficients[5],3))
                    print(paste(gene_comb_name,", HR :",signif(s_comb_gene$coefficients[2],3),", p.value :",signif(s_comb_gene$coefficients[5],3)))
                    
                    # generating KM plots
                    gene_comb_plt<-ggsurvplot(
                        fit = survfit(Surv(OS.time, OS) ~ gene_comb, data=comb_g),
                        break.time.by=500,palette = survival_colors,#pval = T,
                        xlab = "Days",ylab = "Overall survival probability",
                        legend.title= "Gene status",legend.labs=g12[,j],
                        risk.table = T,tables.height = 0.3,
                        font.legends=16,font.tickslab=16,
                        title = gene_comb_km_plot_title)
                   
                    # saving the plot into PNG file 
                    gene_comb_sur_plot_file<-paste(plots_dir,gene_comb_name,"_Group_KM_plot.png",sep = "")
                    png(gene_comb_sur_plot_file,width = 8,height = 6,units = "in",type = "cairo",res = 600)
                    par(xpd=T,mar=c(5,6,4,4))
                    print(gene_comb_plt)
                    dev.off()
                    
                    # saving the plot into SVG file 
                    gene_comb_sur_plot_file_svg<-paste(plots_dir,gene_comb_name,"_Group_KM_plot.svg",sep = "")
                    svg(gene_comb_sur_plot_file_svg,width = 8,height = 6)
                    par(xpd=T,mar=c(5,6,4,4))
                    print(gene_comb_plt)
                    dev.off()
                    
                    # updating the gene comb hr data
                    temp_g<-c(gene_names,g12[,j],s_comb_gene$n,s_comb_gene$nevent,gene_comb_surv_diff$n,gene_comb_surv_diff$obs,signif(s_comb_gene$coefficients,6),signif(s_comb_gene$conf.int[-1],6),signif(gene_comb_surv_diff$chisq,6),signif(gene_comb_pvalue,6),signif(s_comb_gene$waldtest[3],6))
                    gene_comb_hr_analysis<-rbind(gene_comb_hr_analysis,temp_g)
                }
            }
        }
    }
} 
gene_comb_hr_analysis<-na.omit(gene_comb_hr_analysis)
gene_comb_hr_signif<-gene_comb_hr_analysis[which(gene_comb_hr_analysis$`p.value(Z)`<0.05),]

gene_comb_hr_analysis_file<-paste(data_dir,"TCGA_LIHC_Gene_combinations_HR_analysis.tsv",sep = "")
write.table(gene_comb_hr_analysis1,gene_comb_hr_analysis_file1,sep = "\t",col.names = T,row.names = F)

gene_comb_hr_signif_file<-paste(data_dir,"TCGA_LIHC_Gene_combinations_HR_analysis_(significant).tsv",sep = "")
write.table(gene_comb_hr_signif,gene_comb_hr_signif_file,sep = "\t",col.names = T,row.names = F)

## plot of NAPRT and NAD levels ####
naprt_nad_plot_data<-data.frame(NAPRT=t(gene_met_lfc[c("NAPRT"),]),NAD=t(metabolite_lfc["NAD",]),NAD_status=metabolite_status[,"NAD"])

naprt_nad_plot<-ggplot(data = naprt_nad_plot_data,aes(x=NAPRT,y=NAD,color=NAD_status))+
    geom_point(size=2)+scale_color_manual(values=rev(brewer.pal(n = 3, name = 'Set1')))+
    theme_classic()+
    scale_x_continuous(breaks = seq(-6,+3,by=1),limits =  c(-6,+3))+
    scale_y_continuous(breaks = seq(-9,+3,by=1),limits = c(-9,+3))+
    geom_hline(yintercept =  c(-1,+1),lty=2)+
    geom_vline(xintercept =  c(-1,+1),lty=2)+
    theme(axis.text = element_text(size = 16,color = "black"),axis.title = element_text(size = 18,face = "bold"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),legend.position = "bottom",
          legend.text = element_text(size = 14),legend.title = element_text(size = 16))

naprt_nad_plot_file<-paste(plots_dir,"NAPRT_NAD_correlation_plot1.png",sep = "")
png(naprt_nad_plot_file,width = 8,height = 6,units = "in",type = "cairo",res = 600)
par(xpd=T,mar=c(5,6,4,4))
print(naprt_nad_plot)
dev.off()

naprt_nad_plot_file_svg<-paste(plots_dir,"NAPRT_NAD_correlation_plot.svg",sep = "")
svg(naprt_nad_plot_file_svg,width = 5,height = 4)
par(xpd=T,mar=c(5,6,4,4))
print(naprt_nad_plot)
dev.off()

## F factor based survival analysis ####
f_factor_file<-"/home/NAD_project/TCGA_LIHC_ncount_models/TCGA-LIHC_tumour_NAD_F_Factor_data.tsv";
f_factor_data<-read.delim(f_factor_file,header = T)
colnames(f_factor_data)<-gsub("\\.","-",colnames(f_factor_data))

log_ffactor_data<-apply(f_factor_data,2,log2)
f_factor_status<-t(log_ffactor_data)
f_factor_status<-apply(f_factor_status,2, function(x) ifelse(x>1,"Up",ifelse(x<(-1),"Down","No_change")))

survival_ffactor_ids<-na.omit(match(row.names(f_factor_status),survival_data$sample))
ffactor_survival_ids<-survival_data$sample[survival_ffactor_ids]

survival_ffactor_data<-data.frame(survival_data[survival_ffactor_ids,c(3,2,4)],f_factor_status[ffactor_survival_ids,])
row.names(survival_ffactor_data)<-survival_ffactor_data[,1]
survival_ffactor_data<-survival_ffactor_data[,-1]

ffactor_hr_analysis<-as.data.frame(matrix(NA,ncol=19))
colnames(ffactor_hr_analysis)<-c("factor_combination","Group_1","Group_2","N_total","N_event","N_group1","N_group2","Obs(group1)","Obs(group2)","Coef","exp(coef)","SE(coef)","Z","p.value(Z)","exp(-coef)","CI(95_lower)","CI(95_upper)","Chisq","p.value(Chisq)")

for(f1 in 3:23)
{
    for (f2 in 3:23)
    {
        if (!(f1==f2))
        {
            temp_factor<-na.omit(survival_ffactor_data[,c(1:2,f1,f2)])
            
            temp_factor$factor_comb<-apply(temp_factor[,3:4],1, function(x) paste(x,collapse = "_"))
            factor_names<-paste(colnames(temp_factor)[3:4],collapse = "_")
            group_name<-paste("Factor :",factor_names,sep = " ")
            temp_factor$factor_comb<-as.factor(temp_factor$factor_comb)
            
            # rm_sample<-which(temp_factor$factor_comb=="No_change_No_change")
            # if(length(rm_sample)>1)
            # {
            #     temp_factor1<-temp_factor[-rm_sample,]
            # } else { temp_factor1<-temp_factor }
            temp_factor1<-temp_factor 
            temp_factor2<-temp_factor1[,c(1,2,5)]
            temp_factor2[,3]<-factor(temp_factor2[,3])
            if (length(unique(temp_factor2[,3]))>1)
            {
                f12<-combn(levels(temp_factor2$factor_comb),2)
                for (j in 1:ncol(f12))
                {
                    temp_f1<-temp_factor2[grep(f12[1,j],temp_factor2$factor_comb),]
                    temp_f2<-temp_factor2[grep(f12[2,j],temp_factor2$factor_comb),]
                    comb_f<-rbind(temp_f1,temp_f2)
                    comb_f[,3]<-factor(comb_f[,3])
                    cox_comb_factor_function<-coxph(Surv(comb_f$OS.time,comb_f$OS == 1)~ ., data =comb_f)
                    s_comb_factor<-summary(cox_comb_factor_function)
                    
                    factor_comb_name<-paste(factor_names,"(",paste(f12[,j],collapse = ","),")",sep = "")
                    factor_comb_surv_diff <- survdiff(Surv(OS.time, OS) ~ factor_comb, data=comb_f)
                    factor_comb_pvalue<-( 1 - pchisq(factor_comb_surv_diff$chisq, length(factor_comb_surv_diff$n) - 1))
                    factor_comb_km_plot_title<-paste(group_name,", HR :",signif(s_comb_factor$coefficients[2],3),", p.value :",signif(s_comb_factor$coefficients[5],3))
                    print(paste(factor_comb_name,", HR :",signif(s_comb_factor$coefficients[2],3),", p.value :",signif(s_comb_factor$coefficients[5],3)))
                    
                    # # generating KM plot
                    # factor_comb_plt<-ggsurvplot(
                    #     fit = survfit(Surv(OS.time, OS) ~ factor_comb, data=comb_f),
                    #     break.time.by=500,palette = survival_colors,#pval = T,
                    #     xlab = "Days",ylab = "Overall survival probability",
                    #     legend.title= "Factor status",legend.labs=f12[,j],
                    #     risk.table = T,tables.height = 0.3,
                    #     font.legends=16,font.tickslab=16,
                    #     title = factor_comb_km_plot_title)
                    # 
                    # factor_comb_sur_plot_file<-paste(plots_dir,factor_comb_name,"_Group_KM_plot.png",sep = "")
                    # png(factor_comb_sur_plot_file,width = 8,height = 6,units = "in",type = "cairo",res = 600)
                    # par(xpd=T,mar=c(5,6,4,4))
                    # print(factor_comb_plt)
                    # dev.off()
                    # 
                    # factor_comb_sur_plot_file_svg<-paste(plots_dir,factor_comb_name,"_Group_KM_plot.svg",sep = "")
                    # svg(factor_comb_sur_plot_file_svg,width = 8,height = 6)
                    # par(xpd=T,mar=c(5,6,4,4))
                    # print(factor_comb_plt)
                    # dev.off()
                    
                    temp_f<-c(factor_names,f12[,j],s_comb_factor$n,s_comb_factor$nevent,factor_comb_surv_diff$n,factor_comb_surv_diff$obs,signif(s_comb_factor$coefficients,6),signif(s_comb_factor$conf.int[-1],6),signif(factor_comb_surv_diff$chisq,6),signif(factor_comb_pvalue,6))
                    ffactor_hr_analysis<-rbind(ffactor_hr_analysis,temp_f)
                }
                
            }
        }
        
        else
        {
            temp_factor<-na.omit(survival_ffactor_data[,c(1:2,f1)])
            
            factor_names<-colnames(temp_factor)[3]
            group_name<-paste("Factor :",factor_names,sep = " ")
            # rm_sample<-which(temp_factor[,3]=="No_change")
            # if(length(rm_sample)>1)
            # {
            #     temp_factor1<-temp_factor[-rm_sample,]
            # } else { temp_factor1<-temp_factor }
            temp_factor1<-temp_factor 
            temp_factor1[,3]<-factor(temp_factor1[,3])
            if (length(unique(temp_factor1[,3]))>1)
            {
                f12<-combn(levels(temp_factor1[,3]),2)
                for (j in 1:ncol(f12))
                {
                    temp_f1<-temp_factor1[grep(f12[1,j],temp_factor1[,3]),]
                    temp_f2<-temp_factor1[grep(f12[2,j],temp_factor1[,3]),]
                    comb_f<-rbind(temp_f1,temp_f2)
                    comb_f[,3]<-factor(comb_f[,3])
                    cox_comb_factor_function<-coxph(Surv(comb_f$OS.time,comb_f$OS == 1)~ ., data =comb_f)
                    s_comb_factor<-summary(cox_comb_factor_function)
                
                    factor_comb_name<-paste(factor_names,"(",paste(f12[,j],collapse = ","),")",sep = "")
                    
                    colnames(comb_f)[3]<-"factor"
                    factor_comb_surv_diff <- survdiff(Surv(OS.time, OS) ~ factor, data=comb_f)
                    factor_comb_pvalue<-( 1 - pchisq(factor_comb_surv_diff$chisq, length(factor_comb_surv_diff$n) - 1))
                    factor_comb_km_plot_title<-paste(group_name,", HR :",signif(s_comb_factor$coefficients[2],3),", p.value :",signif(s_comb_factor$coefficients[5],3))
                    print(paste(factor_comb_name,", HR :",signif(s_comb_factor$coefficients[2],3),", p.value :",signif(s_comb_factor$coefficients[5],3)))
                    
                    # # generating KM plot
                    # factor_comb_plt<-ggsurvplot(
                    #     fit = survfit(Surv(OS.time, OS) ~ factor, data=comb_f),
                    #     break.time.by=500,palette = survival_colors,
                    #     xlab = "Days",ylab = "Overall survival probability",
                    #     legend.title= "Factor status",legend.labs=levels(temp_factor1$factor),
                    #     risk.table = T,tables.height = 0.3,
                    #     font.legends=16,font.tickslab=16,
                    #     title = factor_comb_km_plot_title)
                    # 
                    #factor_comb_sur_plot_file<-paste(plots_dir,factor_names,"_Group_KM_plot.png",sep = "")
                    # png(factor_comb_sur_plot_file,width = 8,height = 6,units = "in",type = "cairo",res = 600)
                    # par(xpd=T,mar=c(5,6,4,4))
                    # print(factor_comb_plt)
                    # dev.off()
                    # 
                    # factor_comb_sur_plot_file_svg<-paste(plots_dir,factor_comb_name,"_Group_KM_plot.svg",sep = "")
                    # svg(factor_comb_sur_plot_file_svg,width = 8,height = 6)
                    # par(xpd=T,mar=c(5,6,4,4))
                    # print(factor_comb_plt)
                    # dev.off()
                    
                    # updating factor HR analysis dataframe
                    temp_f<-c(factor_names,f12[,j],s_comb_factor$n,s_comb_factor$nevent,factor_comb_surv_diff$n,factor_comb_surv_diff$obs,signif(s_comb_factor$coefficients,6),signif(s_comb_factor$conf.int[-1],6),signif(factor_comb_surv_diff$chisq,6),signif(factor_comb_pvalue,6))
                    ffactor_hr_analysis<-rbind(ffactor_hr_analysis,temp_f)
                }
            }
        }  
    }    
}        
ffactor_hr_analysis<-na.omit(ffactor_hr_analysis)
ffactor_hr_analysis[,4:19]<-apply(ffactor_hr_analysis[4:19],2,as.numeric)
ffactor_hr_analysis_sig<-ffactor_hr_analysis[which(ffactor_hr_analysis$`p.value(Z)`<0.05),]

# writing the HR data into file
ffactor_hr_analysis_file<-paste(data_dir,"TCGA_LIHC_NAD_Ffactor_combination_HR_analysis.tsv",sep="")
ffactor_hr_analysis_file_sig<-paste(data_dir,"TCGA_LIHC_NAD_Ffactor_combination_HR_analysis_(significant).tsv",sep="")

write.table(ffactor_hr_analysis,ffactor_hr_analysis_file,sep = "\t",col.names = T,row.names = F)
write.table(ffactor_hr_analysis_sig,ffactor_hr_analysis_file_sig,sep = "\t",col.names = T,row.names = F)
