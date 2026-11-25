library(igraph)
library(xlsx)
library(tidyverse)

# allelic profile
file_paths <- list.files(path = ".\\data", 
                         pattern = "\\_cgmlst.xlsx", 
                         full.names = TRUE)
  
# combine all the data 
allelicProfile_list <- list()

for (file_path in file_paths){
  file_name=gsub("_cgmlst.xlsx", "",strsplit(file_path, split="/")[[1]][2])
  result = read.xlsx(file_path,sheetName="output")
  
  result <- result %>% # filter out contigs without any match
    filter(Match!="-") %>% 
    arrange(Locus)
  result$Contig <- file_name # file name as contig name
  result$Locus <- sub(" \\(.*", "", result$Locus) # tidy up the locus name
  
  # for duplicated loci, select exact match 
  result <- result %>%
    group_by(Locus) %>%
    arrange(Match != "exact") %>%  # This ensures that "exact" matches come first within each group
    slice_head(n = 1) %>%  # Select the first row within each group after arranging
    ungroup()
  
  result$Allele <- as.numeric(result$Allele)  # if character, convert it into numeric values
  
  result <- result %>% 
    select(c(1,3,4)) %>% 
    pivot_wider(names_from = Locus, values_from = Allele)
 
  allelicProfile_list[[file_name]] <- result
}

#final result
final_result <- bind_rows(allelicProfile_list)

# also import all the environmental strains from pubMLST cgMLST database
library(readxl)
database <- read_excel("data/pubMLST_cgMLST_info.xlsx", sheet = "cgmlst")
colnames(database)[1] <- "Contig"
database <- database %>% 
  select(-c(2,ncol(database)))
database <- as.data.frame(database)


# bind with the cgMLST result above
all(colnames(database) %in% colnames(final_result)) 
all(colnames(final_result) %in% colnames(database))
cgMLST_result <- rbind(database, final_result)
cgMLST_result <- as.data.frame(cgMLST_result)
table(is.na(cgMLST_result)) # check NA values
# Replace NA with the most common allele in each column
for(i in 1:ncol(cgMLST_result)) {
  if(any(is.na(cgMLST_result[, i]))) {
    # Find the most common value, excluding NA
    most_common <- names(sort(table(cgMLST_result[,i], useNA = "no"),decreasing=TRUE))[1]
    # Replace NA values with the most common allele
    cgMLST_result[is.na(cgMLST_result[,i]), i] <- as.numeric(most_common)
  }
}
table(is.na(cgMLST_result)) # no NA values!

rownames(cgMLST_result) <- cgMLST_result$Contig
cgMLST_result$Contig <- NULL
cgMLST_result <- as.matrix(cgMLST_result)
dim(cgMLST_result)

# Pairwise allelic differences 
calculate_difference <- function(profile1, profile2) {
  sum(profile1 != profile2)
}

n <- nrow(cgMLST_result)
difference_matrix <- matrix(0, n, n)

# Fill the matrix with allelic differences
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    difference_matrix[i, j] <- calculate_difference(cgMLST_result[i, ], cgMLST_result[j, ])
    difference_matrix[j, i] <- difference_matrix[i, j]  # Symmetric
  }
}

rownames(difference_matrix) <- rownames(cgMLST_result)
colnames(difference_matrix) <- rownames(cgMLST_result)

write.xlsx(as.data.frame(difference_matrix),"allele_matrix_new.xlsx")


## Create adjacency matrix
adjacency_matrix <- 1 * (difference_matrix <= 133)
rownames(adjacency_matrix)<-rownames(cgMLST_result)
colnames(adjacency_matrix)<-rownames(cgMLST_result)
diag(adjacency_matrix) <- 0 # change diagonal values into 0 


## convert them into edge & node lists (i already have node info, so just create edge list)
#edge list
edge_list <- as.data.frame(as.table(adjacency_matrix))
names(edge_list) <- c("source", "target", "weight")
edge_list <- subset(edge_list, weight > 0)
edge_list$weight <- NULL # don't need weight info
write.csv(edge_list, "edge_list_new.csv", row.names = FALSE)

