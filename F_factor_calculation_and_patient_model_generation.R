### the code is written to compute fold change values and f-factor values for patient specific model generation
# checking for necessary libraries
if(!("BiocManager" %in% installed.packages()))
{
    install.packages("BiocManager")
}

required_bioconductor_packages<-c("edgeR","org.Hs.eg.db")
check_bioconductor_libraries<-(required_bioconductor_packages %in% installed.packages())
bioconductor_libraries_to_install<-required_bioconductor_packages[!check_bioconductor_libraries]
if (length(bioconductor_libraries_to_install)>1)
{
    require(BiocManager)
    BiocManager::install(bioconductor_libraries_to_install)
}
# Loading necessary libraries
for( p in required_bioconductor_packages) require(p,character.only = T) 

# source file for patient model generation
# change the path of the source file based on your system
source("patient_model_generator.R")

# getting necessary output folder paths
# directory to store the tumour models
data_dir<-"/home/NAD_project/TCGA_LIHC_ncount_models/" 
# directory to store the lfc files
deg_dir<-"/home/NAD_project/TCGA_LIHC_ncount_models/edgeR_DEG_analysis/"

# reading necessary input files
count_data_file<-"/home/NAD_project/TCGA_HCC/TCGA-LIHC.htseq_counts.tsv"; # read ccounts file
f_factor_files<-"/home/NAD_project/NAD_F_factors.csv" # F_factor to gene mapping file
nad_gene_file<-"/home/NAD_project/NAD_genes.csv" # NAD net genes file

f_factor_genes<-read.delim(f_factor_files,header = T,strip.white = T,stringsAsFactors = F)
nad_genes<-read.delim(nad_gene_file,header = T,stringsAsFactors = F,strip.white = T)

#### normalized count data generation and analysis ####
# loading count data
count_data<-read.delim(count_data_file,header = F,strip.white = T,stringsAsFactors = F)
colnames(count_data)<-count_data[1,]

# extracting the gene ids and mapping them to gene symbols
ensembl_genes<-count_data[,1]
ensembl_genes<-sapply(ensembl_genes,function(x) unlist(strsplit(x,split = "\\."))[1])
rownames(count_data)<-ensembl_genes
count_data<-count_data[-1,-1]
ensembl_genes<-ensembl_genes[-1]
genesymbols<-data.frame(EnsemblID=ensembl_genes,
                        Symbols=mapIds(org.Hs.eg.db, ensembl_genes,keytype="ENSEMBL", column="SYMBOL"),
                        entrezID=mapIds(org.Hs.eg.db, ensembl_genes,keytype="ENSEMBL", column="ENTREZID"))

# selecting Normal and tumour samples data
sample_ids<-colnames(count_data)
normal_count_col<-grep("-11",sample_ids)
tumour_count_col<-grep("-01",sample_ids)

# subsetting Normal and tumour samples data
normal_data_count<-count_data[,normal_count_col]
tumour_data_count<-count_data[,tumour_count_col]
normal_ids<-colnames(normal_data_count)
tumour_ids<-colnames(tumour_data_count)

# preparing data for DEG analysis using edgeR for each patient sample
sample_type<-factor(c(rep("N",ncol(normal_data_count)),rep("T",ncol(tumour_data_count))))
data_count<-cbind(normal_data_count,tumour_data_count) # combining normla and tumour data
data_count<-apply(data_count,2,as.numeric) # converting to numeric values
data_count<-2^data_count # converting from log2 to normal counts
data_count<-data_count-1 # removing the additional normalizing factor
row.names(data_count)<-ensembl_genes
colnames(data_count)<-c(normal_ids,tumour_ids)

normal_data<-data_count[,normal_ids]
tumour_data<-data_count[,tumour_ids]

nad_genes_lfc_data<-matrix(0,nrow = nrow(nad_genes),ncol = length(tumour_ids))
row.names(nad_genes_lfc_data)<-nad_genes$Gene_symbols
colnames(nad_genes_lfc_data)<-tumour_ids
for (j in 1:length(tumour_ids))
{
    sample_types<-factor(c(rep("N",ncol(normal_data)),"T"))
    data_ncount<-cbind(normal_data,tumour_data[,j]) # combining normal and tumour data
    colnames(data_ncount)<-c(normal_ids,tumour_ids[j])
    design <- model.matrix(~ sample_types)
    # colnames(design)<-levels(sample_type)
    y<-DGEList(counts = data_ncount,group = sample_types,genes = genesymbols)
    # checking if the gene symbols are present in the data
    idfound <- y$genes$Symbols %in% mappedRkeys(org.Hs.egSYMBOL)
    y<-y[idfound,] # subsetting oly gene symbol containing count data
    # Reordering couynt data according to counts
    o <- order(rowSums(y$counts), decreasing=TRUE)
    y <- y[o,]
    d <- duplicated(y$genes$Symbols) # checking for duplicates
    y <- y[!d,] # vremoving duplicate count data
    # nrow(y)
    
    keep <- filterByExpr(y,design = design,min.count=10) # checking for samples and gene with less than 10 count values
    # summary(keep)
    y <- y[keep, , keep.lib.sizes=FALSE] # selecting genes with more than 10 count values
    # gene_sum<-apply(y,1,sum)
    y$samples$lib.size <- colSums(y$counts) # calculating the library sizes
    # setting gene symbols as rownames
    rownames(y$counts)<- rownames(y$genes) <- y$genes$Symbols
    y<-calcNormFactors(y) # calculating the normalizing factors
    norm_counts<-cpm(y)
    y <- estimateDisp(y, design, robust=TRUE,Verbose=T)
    fit_glm <- glmFit(y, design)
    lrt <- glmLRT(fit_glm)
    # topTags(lrt)
    gene_lfc<-data.frame(topTags(lrt,n=nrow(y))$table)
    # summary(decideTests(lrt))
    nad_gene_lfc<-gene_lfc[nad_genes$Gene_symbols,]
    nad_genes_lfc_data[,j]<-nad_gene_lfc$logFC
    # staring the patient specific lfc files
    tumour_lfc_file_name<-paste(deg_dir,tumour_ids[j],"_edgeR_DEG_list.txt",sep = "")
    write.table(nad_gene_lfc,tumour_lfc_file_name,sep = "\t",row.names = F,col.names = T)
}
# saving the all the NAD_net genes from tumour samples
nad_genes_lfc_data_name<-paste(data_dir,"TCGA-LIHC_tumour_NAD_LFC_data_(edgeR).tsv",sep="")
write.table(nad_genes_lfc_data,nad_genes_lfc_data_name,sep = "\t",row.names = T,col.names = T)

# selecting F factors related genes in the data
f_factor_rows<-list()
for(i in 1:nrow(f_factor_genes))
{
    a<-unlist(strsplit(f_factor_genes[i,2],split = ","))
    if(length(a)>1)
    {
        b=c()
        for (j in a)
        {
            b=c(b,grep(j,row.names(nad_genes_lfc_data)))
        }
    } else
    {
        b=grep(a,row.names(nad_genes_lfc_data))
    }
    f_factor_rows<-append(f_factor_rows,list(b))
}
names(f_factor_rows)<-f_factor_genes[,1]

# calculate F factors for each tumour sample 
ffactor_data<-matrix(0,nrow = nrow(f_factor_genes),ncol = ncol(nad_genes_lfc_data))
for (p in 1:ncol(ffactor_data))
{
    patient_data<-2^nad_genes_lfc_data[,p]
    for (i in 1:nrow(ffactor_data))
    {
        ffactor_data[i,p]<-sum(patient_data[unlist(f_factor_rows[i])],na.rm = T)
    }
}
colnames(ffactor_data)<-colnames(nad_genes_lfc_data)
row.names(ffactor_data)<-f_factor_genes[,1]
# saving the f-factor data for each tumour sample
ffactor_data_name<-paste(data_dir,"TCGA-LIHC_tumour_NAD_F_Factor_data.tsv",sep="")
write.table(ffactor_data,ffactor_data_name,sep = "\t",row.names = T,col.names = T)

### generating patient specific NAD Models ####
patient_models_dir<-"/home/NAD_project/TCGA_LIHC_ncount_models/TCGA_LIHC_patient_models/"
# reading the base model 
base_model_file<-"/home/NAD_project/TCGA_LIHC_patient_ncount_models/Normal_base_NAD_model.xml"
# creating tumur speific NAD_net models
patient_model_generator(ffactor_data,base_model_file,patient_models_dir)
