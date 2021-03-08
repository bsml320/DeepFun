setwd("Heatmap_Plot/")

# Load SAD summary for input variants.
SAD_summary = read.delim("JOB3984_screen_example/screen_HD.DeepFun_screen.SAD_summary.txt", head = F)
input_list <- SAD_summary[-1,(which(SAD_summary[1,]=="Otherinfo1")+1):ncol(SAD_summary)]
input_list_tag <- input_list[,3]

# Load global SAD profiles.
SAD = read.delim("JOB3984_screen_example/screen_HD.DeepFun.sad.txt", head = T, colClasses = "character")
SAD_matrix = SAD[-c(1:4), -c(1:4)]
SAD_matrix_df <- matrix(as.numeric(unlist(SAD_matrix)), ncol = ncol(SAD_matrix))
rownames(SAD_matrix_df) <- rownames(SAD_matrix)
#nrow(SAD_matrix) == length(input_list_tag)

# Renamed the variants.
rownames(SAD_matrix_df) = input_list_tag

SAD_matrix_tag = SAD[-c(1:4), 1:4]
colnames(SAD_matrix_tag) <- c("Chr", "Pos", "Ref", "Alt")
SAD_header = SAD[c(1:4),-c(1:4)]
SAD_label <- as.data.frame(t(SAD_header)[,1:3])
colnames(SAD_label) <- c("Experiment_target", "Biosample_type", "Biosample_term_name")

capitalize <- function (string) {
    capped <- grep("^[A-Z]", string, invert = TRUE)
    substr(string[capped], 1, 1) <- toupper(substr(string[capped], 
        1, 1))
    return(string)
}

# Select Experiment target in tissue or cell type.
#table(SAD_label$Experiment_target)
#table(SAD_label$Biosample_type)
target <- "DNase-seq"
type <- "tissue"		
id_target <- which(SAD_label$Experiment_target == target)	
id_type <- which(SAD_label$Biosample_type == type)	
id_overlap <- intersect(id_target, id_type)

SAD_matrix_df_subset <- SAD_matrix_df[, id_overlap]

# Select top 50 variants with the highest SAD value.
#top_variants <- names(head(sort(apply((abs(SAD_matrix_df_subset)), 1, max), decreasing = T), 50))				
#SAD_matrix_df_subset <- SAD_matrix_df_subset[which(rownames(SAD_matrix_df_subset) %in% top_variants), ]

label_subset <- SAD_label[id_overlap, -1]
label_subset <- matrix(as.vector(unlist(label_subset)), ncol = ncol(label_subset))
colnames(label_subset) <- c("Biosample_type", "Biosample_term_name")
Group <- capitalize(label_subset[, 2])
		
color1 <- colorRampPalette(c("white", "red"))(floor(max(SAD_matrix_df_subset) * 10000))
color2 <- colorRampPalette(c("blue", "white"))(floor(min(SAD_matrix_df_subset) * -10000))
color <- c(color2, color1)
set.seed(3)

library(ComplexHeatmap)
ha_subset = HeatmapAnnotation(df = data.frame(Group))

pdf("DNase-seq_tissue.pdf", 20, 11)	
print(Heatmap(SAD_matrix_df_subset, name = "SAD", top_annotation = ha_subset, cluster_columns = FALSE, show_column_names = FALSE, cluster_rows = FALSE,
	show_row_names = TRUE, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8), col = color) )
dev.off()

# Merge SAD profile from duplicate profiles.

SAD_matrix_merged = matrix(0, nrow = nrow(SAD_matrix_df_subset), ncol = length(table(Group)))
rownames(SAD_matrix_merged) = rownames(SAD_matrix_df_subset)
colnames(SAD_matrix_merged) = names(table(Group))

for (i in 1:nrow(SAD_matrix_merged)){
	SAD_matrix_merged[i,] <- tapply(SAD_matrix_df_subset[i,], Group, median)
}

write.table(SAD_matrix_merged, "SAD_matrix_merged.txt", sep = "\t", quote = F)

# Find the most relevant tissue or cell type.

relevant_tissue <- matrix(NA, nrow = nrow(SAD_matrix_df_subset), ncol = 5)
rownames(relevant_tissue) = rownames(SAD_matrix_df_subset)

for (i in 1:5) {
        relevant_tissue[, i] <- colnames(SAD_matrix_merged)[apply(abs(SAD_matrix_merged), 1, order)[ncol(SAD_matrix_merged) - i + 1, ]]
}
colnames(relevant_tissue) = paste0("Top", 1:ncol(relevant_tissue))

write.table(relevant_tissue, "relevant_tissue.txt", sep = "\t", quote = F)
