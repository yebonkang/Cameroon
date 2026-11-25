library(dplyr)
library(stringr)
library(xlsx)

file_paths <- list.files(path = "results-clade2-new", 
                         pattern = "\\.tsv$", 
                         full.names = TRUE)

#check if NA values 
na_vector=c()
for (file_path in file_paths){
  result = read.delim(file_path, as.is=T, check.names=F)
  na_vector <- c(na_vector, any(is.na(result)))
}
print(ifelse(any(na_vector)==TRUE, "NA included", "No NAs")) # check if any NA value -> return TRUE (meaning 

result_list <- list()

for (file_path in file_paths){
  file_name=strsplit(file_path, split="/")[[1]][2]
  result = read.delim(file_path, as.is=T, check.names=F)
  
  result  <- result  %>%
    mutate(strain = case_when(
      str_detect(sseqid, "^N16961") ~ "N16961",
      str_detect(sseqid,"Cam23") ~ "Cam23",
      str_detect(sseqid,"TUC_T2734") ~ "TUC_T2734",
      str_detect(sseqid,"2010V-1116") ~ "2010V-1116",
      str_detect(sseqid, "_") ~ str_extract(sseqid, ".*?(?=_)")
    )) 

  data = result %>% 
    group_by(strain) %>% 
    slice_max(score, n = 1, with_ties = FALSE) %>%
    arrange(desc(score))
  
  reference=max(data$score) # reference score 
  
  # Calculate BSR
  data = data %>% 
    mutate(bsr=score/reference) 
  
  result_list[[file_name]] <- data
}

#final result
final_result <- bind_rows(result_list)

# write at once 
write.xlsx(as.data.frame(final_result), "clade2-new-bsr_result.xlsx", row.names=FALSE)

# Create a matrix for BSR scores
n=length(unique(final_result$qseqid)) # locus tags = 155: same as no. of original aa sequence 
n_strain=length(unique(final_result$strain))# strains = 8 

bsr_table=matrix(NA, nrow = n_strain, ncol = n)

#strain order
strain_labels<-read.table("./tree/strain_order_clade2-new.txt", header =FALSE)
strain_labels=strain_labels$V1
all(strain_labels %in% unique(final_result$strain))
all(unique(final_result$strain) %in% strain_labels)

rownames(bsr_table)=strain_labels # edit row names
colnames(bsr_table)=unique(final_result$qseqid) #edit col names
all(final_result$qseqid %in% colnames(bsr_table))
all(final_result$strain %in% rownames(bsr_table))

# fill BSR score to bsr_table matrix
for (row_name in rownames(bsr_table)) {
  for (col_name in colnames(bsr_table)) {
    # Find matching location 
    df_index <- which(final_result$strain == row_name & final_result$qseqid == col_name)
    
    if (length(df_index)==1){
      bsr_table[row_name, col_name] <- final_result[df_index, "bsr"][[1]]
    }
  }
}

bsr_table[is.na(bsr_table)]=0 # fill NA values with "0"



##############################################################
################HEAT MAP for each genomic island##############
library(ComplexHeatmap)

library(pheatmap)
#extract locus tags for each genomic island
VSP1=colnames(bsr_table)[1:11] # VSP1
VSP2=colnames(bsr_table)[13:36] # VSP2
VPI1=colnames(bsr_table)[38:74] # VPI1
ctx=colnames(bsr_table)[77:103] # CTX-Phi cluster 
VPI2=colnames(bsr_table)[104:(length(colnames(bsr_table))-2)] # VPI2
virulence_genes = c("VC_0409", "VC_0633","VC_0984","VC_1450","VC_A0219","VC_A0865")
idx=colnames(bsr_table)%in%virulence_genes
Others=colnames(bsr_table)[idx] # Other virulence genes 

#color palette
color_pal=colorRampPalette(c("white","red3"))
color_pal=color_pal(100)

### CTX-Phi cluster 
ctx_map <- Heatmap(bsr_table[,ctx],
                   col=color_pal,
                   cluster_rows=FALSE,cluster_columns = FALSE,
                   row_names_side = "left",
                   column_names_side = "top",
                   column_names_gp = gpar(fontsize = 10),
                   border=FALSE,
                   rect_gp = gpar(col = "grey1"),
                   name="BSR",
                   heatmap_height = unit(1, "npc")*nrow(bsr_table[,ctx]),
                   heatmap_width=unit(0.7,"cm")*ncol(bsr_table[,ctx]),
                   height=50
)
plot(ctx_map)

### VSP1 
VSP1_map <- Heatmap(bsr_table[,VSP1],
                    col=color_pal,
                    cluster_rows=FALSE,cluster_columns = FALSE,
                    row_names_side = "left",
                    column_names_side = "top",
                    column_names_gp = gpar(fontsize = 10),
                    border=FALSE,
                    rect_gp = gpar(col = "grey1"),
                    name="BSR",
                    heatmap_height = unit(1, "npc")*nrow(bsr_table[,VSP1]),
                    heatmap_width=unit(0.9,"cm")*ncol(bsr_table[,VSP1]),
                    height=50
)
plot(VSP1_map)

### VSP2
VSP2_map <- Heatmap(bsr_table[,VSP2],
                    col=color_pal,
                    cluster_rows=FALSE,cluster_columns = FALSE,
                    row_names_side = "left",
                    column_names_side = "top",
                    column_names_gp = gpar(fontsize = 10),
                    border=FALSE,
                    rect_gp = gpar(col = "grey1"),
                    name="BSR",
                    heatmap_height = unit(1, "npc")*nrow(bsr_table[,VSP2]),
                    heatmap_width=unit(0.7,"cm")*ncol(bsr_table[,VSP2]),
                    height=50
)
plot(VSP2_map)

### VPI1
VPI1_map <- Heatmap(bsr_table[,VPI1],
                    col=color_pal,
                    cluster_rows=FALSE,cluster_columns = FALSE,
                    row_names_side = "left",
                    column_names_side = "top",
                    column_names_gp = gpar(fontsize = 10),
                    border=FALSE,
                    rect_gp = gpar(col = "grey1"),
                    name="BSR",
                    heatmap_height = unit(1, "npc")*nrow(bsr_table[,VPI1]),
                    heatmap_width=unit(0.5,"cm")*ncol(bsr_table[,VPI1]),
                    height=50
)
plot(VPI1_map)

### VPI2
VPI2_map <- Heatmap(bsr_table[,VPI2],
                    col=color_pal,
                    cluster_rows=FALSE,cluster_columns = FALSE,
                    row_names_side = "left",
                    column_names_side = "top",
                    column_names_gp = gpar(fontsize = 10),
                    border=FALSE,
                    rect_gp = gpar(col = "grey1"),
                    name="BSR",
                    heatmap_height = unit(1, "npc")*nrow(bsr_table[,VPI2]),
                    heatmap_width=unit(0.5,"cm")*ncol(bsr_table[,VPI2]),
                    height=50
)
plot(VPI2_map)

### Others
Others_map <- Heatmap(bsr_table[,Others],
                    col=color_pal,
                    cluster_rows=FALSE,cluster_columns = FALSE,
                    row_names_side = "left",
                    column_names_side = "top",
                    column_names_gp = gpar(fontsize = 10),
                    border=FALSE,
                    rect_gp = gpar(col = "grey1"),
                    name="BSR",
                    heatmap_height = unit(1, "npc")*nrow(bsr_table[,Others]),
                    heatmap_width=unit(0.9,"cm")*nrow(bsr_table[,Others]),
                    height=50
)
plot(Others_map)
