## ---------------------------
##
## Script name: 
##
## Purpose of script:
##
## Author: Chathura J. Gunasekara, PhD
##
## Date Created: 2020-02-26
##
## Copyright (c) Chathura J Gunasekara, 2020
## Email: gunaseka@bcm.edu/cjgunase@mtu.edu
##
## ---------------------------
##
## Notes:
##   	Each 100bp bin should consist of at least 1 CpG site and 
#satisfy the following coverage requirement:
#If number of CpGs (nCpGs) <= 2, nCpGs sites should be covered by N reads
#If nCpGs >2 , ⌈nCpG/2 ⌉ CpG sites should be covered by N reads.

##
## ---------------------------

## set working directory for Mac and PC

setwd("~/Documents/CoRSIV_Capture")          # Chathura's working directory (mac)


## ---------------------------

options(scipen = 6, digits = 4) # I prefer to view outputs in non-scientific notation

## ---------------------------
## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(tidyr)
require(data.table)
require(gplots)
require(qqman)
require(pracma)
require(imputeTS)
require(corrplot)
## ---------------------------

READ_DEPTH_CUTOFF = 10

USC_Sample_ID_table <- read.table("./USC_pediatric_glioma_alignments/sample_id_table1.txt",header = T)
Covariate_table <- read.table("./USC_pediatric_glioma_alignments/Baylor_covariate_file.txt",header = T)

temp_table <- merge(USC_Sample_ID_table,Covariate_table,by.x = "ID2",by.y = "Blindid")
temp_table$merged_col <- paste(temp_table$ID1,temp_table$Case_status,sep = "_")
temp_table <- temp_table[order(temp_table$ID1),]


files <- "./USC_pediatric_glioma_bed_file"
file_names <- dir(files,full.names = T,pattern = "*.bed")
df <- do.call(rbind,
              lapply(file_names, 
                     function(x) cbind(read.table(x),
                                       #name=strsplit(strsplit(x,'\\.')[[1]][2],"/")[[1]][6]
                                       name=paste0("USC_",strsplit(strsplit(x,split = "/")[[1]][3],split = "_")[[1]][2])
                     )
              )
)

individuals <- as.character(unique(df$name))
corsivs <- as.character(unique(df$V10))
length(corsivs)
corsiv_methy_df <- data.frame(matrix(nrow = length(corsivs),ncol = 0))
corsiv_depth_df <- data.frame(matrix(nrow = length(corsivs),ncol = 0))

for (individual in individuals){
    #individual <- "USC_08"
    ind_df <- df[df$name==individual,]
    corsiv_methylation_values <-c()
    avg_depth_corsivs <- c()
    for(corsiv in corsivs){
        #corsiv="174aab9d"
        print(paste0(individual,"--",corsiv))
        corsiv_df <- ind_df[ind_df$V10==corsiv,]
        corsiv_df$total <- corsiv_df$V5 + corsiv_df$V6
        nCpGs <- dim(corsiv_df)[1]
        
        if(nCpGs>2){
            corsiv_df <- corsiv_df[corsiv_df$total>=READ_DEPTH_CUTOFF,]
            if(dim(corsiv_df)[1] > ceil(nCpGs/2)){
                avg_depth_corsivs <- c(avg_depth_corsivs,1)
            }else{
                avg_depth_corsivs <- c(avg_depth_corsivs,0)
            }
            
        }else if(0 < nCpGs & nCpGs <=2 ){
            corsiv_df <- corsiv_df[corsiv_df$total>=READ_DEPTH_CUTOFF,]
            if(dim(corsiv_df)[1] == nCpGs){
                avg_depth_corsivs <- c(avg_depth_corsivs,1)
            }else{
                avg_depth_corsivs <- c(avg_depth_corsivs,0)
            }
        }else{
            avg_depth_corsivs <- c(avg_depth_corsivs,0)
        }
        
        corsiv_methylation_values <- c(corsiv_methylation_values, (sum(corsiv_df$V5)/sum(corsiv_df$total))*100)

    }
    corsiv_methy_df[individual] <- corsiv_methylation_values
    corsiv_depth_df[individual] <- avg_depth_corsivs
}
rownames(corsiv_methy_df) <- corsivs
rownames(corsiv_depth_df) <- corsivs
   

########################################### Depth analysis ############

dim(corsiv_depth_df)[1]
colnames(corsiv_depth_df) <- temp_table$merged_col

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
depth_data <- as.matrix(corsiv_depth_df)
#depth_data <- depth_data>=10
#depth_data[is.na(depth_data)] <- 0
#depth_data <- 1*depth_data

heatmap.2(depth_data,trace="none",dendrogram = "row",margins=c(12,8),col=my_palette,main = "read depth cutoff >=25",Colv="NA")

corsiv_covered_by_samples <- data.frame(rowSums(depth_data))
corsiv_covered_by_samples$CoRSIV_Id <- rownames(corsiv_covered_by_samples)
rownames(corsiv_covered_by_samples) <- NULL
colnames(corsiv_covered_by_samples) <- c("nSamples","CoRSIV_Id")

barplot((table(corsiv_covered_by_samples$nSamples)/dim(corsiv_depth_df)[1])*100,ylim = c(0,100),
        ylab = "% ",xlab="Number of Samples (n)", main = " % of CoRSIVs satisfy coverage criteria in 24 samples")


selected_corsivs <- corsiv_covered_by_samples[corsiv_covered_by_samples$nSamples>19,]$CoRSIV_Id


####################################### Methylation analysis ##########

corsiv_methy_df_subset <- subset(corsiv_methy_df,rownames(corsiv_methy_df) %in% selected_corsivs)

colnames(corsiv_methy_df_subset) <-  temp_table$merged_col

corsiv_methy_df_subset <- na_mean(corsiv_methy_df_subset)
heatmap.2(as.matrix(corsiv_methy_df_subset),trace="none",
          dendrogram = "row",margins=c(12,8),col=my_palette,main = "hierarchical clustering of Capture Methylation Data",Colv = NA)




# range(corsiv_methy_df[1,])
# temp <- data.frame(t(apply(corsiv_methy_df,1,range)))
# temp$IIR <- abs(temp$X1 - temp$X2)
# 
# temp$CoRSIV_ID <- rownames(temp)
# 
# df <- df[!duplicated(df$V10),]
# df <- df[c("V10","V12")]
# 
# final <- merge(temp,df,by.x = "CoRSIV_ID",by.y = "V10")
# cor(final$IIR,final$V12,method = "pearson")
# plot(final$IIR,final$V12,xlim = c(0,100),ylim = c(0,100),ylab = "CoRSIV Methy. Range in GTeX 10 indi.",xlab="CoRSIV Methy. Range in USC_pediatric 24 individuals",col="blue")
# 
# 
# 
# colnames(corsiv_methy_df) <- table_1_2$merged_col
# 

# 
# png("./heatmaps_in_r.png",    # create PNG for the heat map        
#     width = 5*300,        # 5 x 300 pixels
#     height = 8*300,
#     res = 300,            # 300 pixels per inch
#     pointsize = 8)
# my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
# 
# depth_data <- as.matrix(corsiv_methy_df)
# depth_data <- depth_data>=50
# depth_data[is.na(depth_data)] <- 0
# depth_data <- 1*depth_data
# 
# heatmap.2(depth_data,trace="none",Colv = "NA",dendrogram = "row",margins=c(12,8),col=my_palette)
# 
# dev.off()

corsiv_methy_df_subset <- corsiv_methy_df_subset[c(
     sort(grep("control",colnames(corsiv_methy_df_subset))),sort(grep("case",colnames(corsiv_methy_df_subset)))
     )]

########################### t-test comparison ################
t.result <- apply(corsiv_methy_df_subset, 1, function (x) t.test(x[1:12],x[13:24],paired=TRUE))
corsiv_methy_df_subset$p_value <- unlist(lapply(t.result, function(x) x$p.value))
corsiv_methy_df_subset$p_value <- corsiv_methy_df_subset$p_value
corsiv_methy_df_subset$fdr <- p.adjust(corsiv_methy_df_subset$p_value, method = "holm")
corsiv_methy_df_subset$CoRSIV_ID <- rownames(corsiv_methy_df_subset)

complete_bed_file <- read.table("./combined_corsivs.bed",sep = "\t")
dim(complete_bed_file)



to_manhattan <- corsiv_methy_df_subset[c("CoRSIV_ID","p_value")]
to_manhattan <- merge(to_manhattan,complete_bed_file,by.x = "CoRSIV_ID",by.y = "V4")
to_manhattan <- to_manhattan[c("CoRSIV_ID","V1","V2","p_value")]
colnames(to_manhattan) <- c("SNP","CHR","BP","P")
#strtoi(substring(to_manhattan$CHR, 4),10)
to_manhattan$CHR <- strtoi(substring(to_manhattan$CHR, 4),10)
to_manhattan$BP <- as.numeric(to_manhattan$BP)
to_manhattan <- na.omit(to_manhattan)
options(repr.plot.width=18, repr.plot.height=16)
manhattan(to_manhattan, chr="CHR", bp="BP", snp="SNP", p="P",
          suggestiveline = FALSE,col = c("red","blue"),
          annotatePval = 0.001,annotateTop = F)


##############################################################
complete_bed_file <- read.table("./combined_corsivs.bed",sep = "\t")
corsiv_methy_df_subset$CoRSIV_ID <- rownames(corsiv_methy_df_subset)
corsiv_methy_df_subset <- merge(corsiv_methy_df_subset,complete_bed_file,by.x ="CoRSIV_ID",by.y = "V4",no.dups = T,all.x = T)
corsiv_methy_df_subset$V1 <- strtoi(substring(corsiv_methy_df_subset$V1, 4),10)
corsiv_methy_df_subset <- corsiv_methy_df_subset[
    with(corsiv_methy_df_subset, order(V1,V3)),
    ]
for(i in 1:22){
    i=20
    corsiv_methy_df_subset_chr <- corsiv_methy_df_subset[corsiv_methy_df_subset$V1==i,]
    corsiv_methy_df_subset_chr <- na.omit(corsiv_methy_df_subset_chr)
    
    corsiv_methy_df_subset_chr_control <- corsiv_methy_df_subset_chr[c(2:13,31)]
    corsiv_methy_df_subset_chr_case <- corsiv_methy_df_subset_chr[c(14:25,31)]
    
    rownames(corsiv_methy_df_subset_chr_control) <- corsiv_methy_df_subset_chr_control$V5
    corsiv_methy_df_subset_chr_control$V5 <- NULL
    M_control<-cor(t(corsiv_methy_df_subset_chr_control),use="complete.obs",method = "spearman") 
    
    rownames(corsiv_methy_df_subset_chr_case) <- corsiv_methy_df_subset_chr_case$V5
    corsiv_methy_df_subset_chr_case$V5 <- NULL
    M_case<-cor(t(corsiv_methy_df_subset_chr_case),use="complete.obs",method = "spearman") 
    
    M_case[lower.tri(M_case,diag = T)] <-0
    M_control[upper.tri(M_control,diag = T)] <-0
    M <- M_control + M_case
    
    #pdf(file=paste("./pearson_corrplots/USC_Pearson_corrplot_Controls_glioma_Cases_",i,".pdf",sep = ""),family = "Courier",width = 50,height = 50)
    pdf(file=paste("./spearman_corrplots/USC_Spearman_corrplot_Controls_glioma_Cases_",i,".pdf",sep = ""),family = "Courier",width = 50,height = 50)
    
    corrplot(M,
             method = "square", tl.cex = 1, tl.col = 'black',
             diag = FALSE,tl.srt=45,addgrid.col = NA,is.corr = FALSE,cl.lim = c(-1,1))#addCoef.col = "black",number.cex = 0.2,
    dev.off()
}




library(psych)
library("PerformanceAnalytics")
case_focus_region <- t(corsiv_methy_df_subset_chr_case[89:133,])
M_case<-cor(case_focus_region,use="complete.obs",method = "pearson")

control_focus_region <- t(corsiv_methy_df_subset_chr_control[89:133,])
M_control<-cor(control_focus_region,use="complete.obs",method = "pearson")

M_case[lower.tri(M_case,diag = T)] <-0
M_control[upper.tri(M_control,diag = T)] <-0
M <- M_control + M_case

pdf(file=paste("./pearson_corrplots/Pearson_correlation_plot_control_vs_case_focus_region_in_chr",i,".pdf",sep = ""),family = "Courier",width = 15,height = 15)
corrplot(M,
         method = "square", tl.cex = 1, tl.col = 'black',
         diag = FALSE,tl.srt=45,addgrid.col = NA,is.corr = FALSE,cl.lim = c(-1,1),addCoef.col = "black",
         number.cex = 0.3,
         title = "Pearson : chr20:29306401-29307400 /n to chr20:30747301-30747500")
dev.off()



pdf(file=paste("./pearson_corrplots/Pearson_Case_Scatter_plots_Focus_region_in_chr",i,".pdf",sep = ""),family = "Courier",width = 50,height = 50)

pairs.panels(case_focus_region,
             smooth = F, scale = T,lm=T,cor = T,stars = T,
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = F,  # show density plots
             ellipses = F # show correlation ellipses,
)


dev.off()


# ?chart.Correlation
# 
# v1 <- as.numeric(corsiv_methy_df_subset_chr_control[rownames(corsiv_methy_df_subset_chr_control)=="chr20:29752801-29753300",])
# v2 <- as.numeric(corsiv_methy_df_subset_chr_control[rownames(corsiv_methy_df_subset_chr_control)=="chr20:29753301-29754200",])
# cor.test(v1,v2,method = "pearson")
# plot(v1,v2,xlab = "chr20:29752801-29753300",ylab="chr20:29753301-29754200",main = "Control samples",col="black",pch=16)
# 
# 
# 
# 
# 
# CoRSIV1_all <- df[df$V11=="chr20:29752801-29753300",]
# CoRSIV1_all <- merge(CoRSIV1_all,temp_table,by.x="name",by.y = "ID1")
# CoRSIV1_all$tatal_read_count <- CoRSIV1_all$V5 + CoRSIV1_all$V6
# CoRSIV1_all_cases <- CoRSIV1_all[CoRSIV1_all$Case_status=="case",]
# CoRSIV1_all_control <- CoRSIV1_all[CoRSIV1_all$Case_status=="control",]
# 
# 
# 
# control_depth_R1 <- data.frame( 
#     CoRSIV1_all_control  %>%
#     group_by(name) %>%
#     summarise(mean_read_depth = mean(tatal_read_count))
# )
# 
# ggplot(control_depth_R1, aes(x=name, y=mean_read_depth, fill=name)) +
#     geom_bar(stat="identity")+theme_minimal()+ylim(0,75)
# 
# 
# case_depth_R1 <- data.frame(
#     CoRSIV1_all_cases  %>%
#     group_by(name) %>%
#     summarise(mean_read_depth = mean(tatal_read_count))
# )
# 
# ggplot(case_depth_R1, aes(x=name, y=mean_read_depth, fill=name)) +
#     geom_bar(stat="identity")+theme_minimal()+ylim(0,75)
# 
# 
# 
# 
# CoRSIV1_all <- df[df$V11=="chr20:29753301-29754200",]
# CoRSIV1_all <- merge(CoRSIV1_all,temp_table,by.x="name",by.y = "ID1")
# CoRSIV1_all$tatal_read_count <- CoRSIV1_all$V5 + CoRSIV1_all$V6
# CoRSIV1_all_cases <- CoRSIV1_all[CoRSIV1_all$Case_status=="case",]
# CoRSIV1_all_control <- CoRSIV1_all[CoRSIV1_all$Case_status=="control",]
# 
# 
# 
# control_depth_R2 <- data.frame( 
#     CoRSIV1_all_control  %>%
#         group_by(name) %>%
#         summarise(mean_read_depth = mean(tatal_read_count))
# )
# 
# ggplot(control_depth_R2, aes(x=name, y=mean_read_depth, fill=name)) +
#     geom_bar(stat="identity")+theme_minimal()+ylim(0,75)
# 
# 
# case_depth_R2 <- data.frame(
#     CoRSIV1_all_cases  %>%
#         group_by(name) %>%
#         summarise(mean_read_depth = mean(tatal_read_count))
# )
# 
# ggplot(case_depth_R2, aes(x=name, y=mean_read_depth, fill=name)) +
#     geom_bar(stat="identity")+theme_minimal()+ylim(0,75)
# 
# 
# plot(control_depth_R1$mean_read_depth,control_depth_R2$mean_read_depth,xlab = "chr20:29752801-29753300",ylab="chr20:29753301-29754200",main = "control depth",col="black",pch=16)
# 
# plot(case_depth_R1$mean_read_depth,case_depth_R2$mean_read_depth,xlab = "chr20:29752801-29753300",ylab="chr20:29753301-29754200",main = "case depth",col="black",pch=16)
# 
# 
# 
# v1 <- as.numeric(corsiv_methy_df_subset_chr_case[rownames(corsiv_methy_df_subset_chr_case)=="chr20:29752801-29753300",])
# v2 <- as.numeric(corsiv_methy_df_subset_chr_case[rownames(corsiv_methy_df_subset_chr_case)=="chr20:29753301-29754200",])
# 
# cor.test(v1,v2,method = "pearson")
# plot(v1,v2,xlab = "chr20:29752801-29753300",ylab="chr20:29753301-29754200",main = "Case samples",col="black",pch=16)



