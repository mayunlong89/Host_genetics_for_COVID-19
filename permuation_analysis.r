#Rscript
#Author: Yunlong Ma
#E-mail: glb-biotech@zju.edu.cn
#100,000 times in silico permutation analysis for genes identified from MAGMA and S-MultiXcan analysis


#Set the work directory
setwd("F:\\Desktop\\")
set.seed(12345)

#Part I Read data on significant genes and background genes

#Read significant genes from MAGMA gene-based association analysis
magma_sig <- read.table("MAGMA_genes.txt", header=T)
magma_gene <- magma_sig$Gene_name

#Read significant genes from S-MultiXcan integrative genomics analysis
multixcan_sig <- read.table("S_MultiXcan.txt", header=T)
multixcan_gene <- multixcan_sig$Gene_name

#Read background genes of S-MultiXcan integrative genomics analysis
Backgroud_multixcan<- read.table("S_MultiXcan_background.txt", header=T)
Backgroud_gene <- Backgroud_multixcan$Gene_name


#Calculate the numebr of genes in each gene set
len_Sig_magma <- length(magma_gene)
len_Sig_multixcan <- length(multixcan_gene)
len_Backgroud_gene <- length(Backgroud_gene)


#Part II establish a function for permutation analysis

#Permutation Function
Permut_analysis <- function(x,y,z){
  
  random_selected_genes <- sample(x,y)
  
  temp <- match(random_selected_genes,z)
  
  random_overlaped_genes <- na.omit(temp)
  
  num<-length(random_overlaped_genes)
  
  return(num)
  
}


#100000 times permutation analysis for MAGMA and S-MultiXcan analysis
results_perumt <- replicate(100000,Permut_analysis(Backgroud_gene,len_Sig_multixcan,magma_gene))


#Ploting function
Fig_random <- function(x,y,z){
  
  hist(x, col="red",xlab="Counts of overlapped genes",xlim = c(0,320),main=NULL)
  temp1 <- match(y,z)
  Observed_overlaped_genes <- na.omit(temp1)
  Observed_gene_num <- length(Observed_overlaped_genes)
  abline(v=Observed_gene_num,col="darkblue",lty="longdash")
  P_value=length(x[x>Observed_gene_num])/length(x)
  x1= Observed_gene_num
  freq <- table(x)
  y1 = max(freq)
  text(x1,y1,P_value)
  
}



#Visulization for MAGMA and S-MultiXcan analysis
Fig_random(results_perumt,multixcan_gene,magma_gene)




#End
