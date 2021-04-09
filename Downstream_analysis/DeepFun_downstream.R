setwd("Downstream_analysis")

# 1. Load DeepFun predition results.
# 1.1 Load the SAD scores for all chromatin features.
SAD = read.delim("JOB6929_3682845047/screen_HD.DeepFun.sad.txt", head = T, colClasses = "character")
SAD_matrix = SAD[-c(1:4), -c(1:4)]
SAD_matrix_df <- matrix(as.numeric(unlist(SAD_matrix)), ncol = ncol(SAD_matrix))
rownames(SAD_matrix_df) <- rownames(SAD_matrix)

# 1.2 Load accesson quality information, including the AUROC and AUPRC result, for downstream QC.
Accession_quality <- read.delim("Accession_quality.txt", head = T, row.names = 1)

# 1.3 Load SAD scores summary, as it saved the raw input information. 
SAD_summary = read.delim("JOB6929_3682845047/screen_HD.DeepFun_screen.SAD_summary.txt", head = F, colClasses = "character")
# Renamed the input variant if it is necessary. Here we use the variant id as example (the 3th column in raw input file).
input_list <- SAD_summary[-1,(which(SAD_summary[1,] == "Otherinfo1") + 1):ncol(SAD_summary)]
input_list_tag <- input_list[,3]
rownames(SAD_matrix_df) = input_list_tag

# Prepare necessary file for downstream heatmap plot.
SAD_matrix_tag = SAD[-c(1:4), 1:4]
colnames(SAD_matrix_tag) <- c("Chr", "Pos", "Ref", "Alt")
SAD_header = SAD[c(1:4),-c(1:4)]
SAD_label <- as.data.frame(t(SAD_header)[,1:3])
SAD_label[,4] <- Accession_quality$AUROC[match(SAD_header[4,], rownames(Accession_quality))]
SAD_label[,5] <- Accession_quality$AUPRC[match(SAD_header[4,], rownames(Accession_quality))]
colnames(SAD_label) <- c("Experiment_target", "Biosample_type", "Biosample_term_name", "AUROC", "AUPRC")

capitalize <- function (string) {
    capped <- grep("^[A-Z]", string, invert = TRUE)
    substr(string[capped], 1, 1) <- toupper(substr(string[capped], 
        1, 1))
    return(string)
}

# 2. DeepFun results for heatmap plot.
# 2.1 Select Experiment target in tissue or cell type.
# Specifiy the experiment target, e.g. DNase-seq, H3K27ac-human. More details: table(SAD_label$Experiment_target)	
target <- "DNase-seq"
# Specifiy the biosample type, e.g. tissue, cell_line. More details: table(SAD_label$Biosample_type)
type <- "tissue"
# Find the overlapped profile.
id_interested <- intersect(which(SAD_label$Experiment_target == target), which(SAD_label$Biosample_type == type))

# 2.2 QC step to remove out low quality profile.
id_HQ <- intersect(which(SAD_label$AUROC > 0.8), which(SAD_label$AUPRC > 0.4))
id_overlap <- intersect(id_interested, id_HQ)

# 2.3 Get the final interested target passed the quality control, for downstream heatmap plot
SAD_score <- SAD_matrix_df[, id_overlap]
colnames(SAD_score) <- as.character(SAD_header[4, id_overlap])

# 2.4 Variants filter, e.g. only use the top 50 variants with the highest SAD value.
#top_variants <- names(head(sort(apply((abs(SAD_score)), 1, max), decreasing = T), 50))				
#SAD_score <- SAD_score[which(rownames(SAD_score) %in% top_variants), ]

label_subset <- SAD_label[id_overlap, -1]
label_subset <- matrix(as.vector(unlist(label_subset)), ncol = ncol(label_subset))
colnames(label_subset) <- c("Biosample_type", "Biosample_term_name", "AUROC", "AUPRC")
Group <- capitalize(label_subset[, 2])

# Load ComplexHeatmap library			
library(ComplexHeatmap)
ha_subset = HeatmapAnnotation(df = data.frame(Group))
color1 <- colorRampPalette(c("white", "red"))(floor(max(SAD_score) * 10000))
color2 <- colorRampPalette(c("blue", "white"))(floor(min(SAD_score) * -10000))
color <- c(color2, color1)
set.seed(3)

# Heatmap plot
pdf("DNase_tissue.pdf", 20, 11)	
print(Heatmap(SAD_score, name = "SAD", top_annotation = ha_subset, cluster_columns = FALSE, show_column_names = FALSE, cluster_rows = FALSE,
	show_row_names = TRUE, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8), col = color) )
dev.off()

# 3. Others downstream analysis.
# 3.1 Evaluate the significant changes of variants base on SAD score.
# To better interpret and take advantage of the SAD scores, user can transform SAD value to standard scale, e.g. Z-score.
SAD_z_score = SAD_p_value = SAD_score
for (i in 1:ncol(SAD_score)){
	SAD_z_score[,i] <- (SAD_score[,i] - mean(SAD_score[,i]))/sd(SAD_score[,i])
	for (j in 1:nrow(SAD_score)){
		SAD_p_value[j, i] <- (1 - pnorm(abs(SAD_z_score[j, i]))) * 2
	}
}

# It is note that DeepFun only accepts up to 3000 variants per job, compare to million level variants in most GWAS study, this is a small proportion.
# We hope user combine all interested and background variant to togethoer before variant significance evaluation. 

# 3.2 Merge SAD scores from duplicate profiles.
SAD_score_merged = matrix(0, nrow = nrow(SAD_score), ncol = length(table(Group)))
rownames(SAD_score_merged) = rownames(SAD_score)
colnames(SAD_score_merged) = names(table(Group))

# Calculate the median SAD score.
for (i in 1:nrow(SAD_score_merged)){
	SAD_score_merged[i,] <- tapply(SAD_score[i,], Group, median)
}

# Save the median SAD scores for interested profiles.
write.table(SAD_score_merged, "SAD_score_merged.txt", sep = "\t", quote = F)

# 3.3 Find the top relevant tissue or cell type for interested experiment target.
relevant_tissue <- matrix(NA, nrow = nrow(SAD_score), ncol = 5)
rownames(relevant_tissue) = rownames(SAD_score)

for (i in 1:5) {
        relevant_tissue[, i] <- colnames(SAD_score_merged)[apply(abs(SAD_score_merged), 1, order)[ncol(SAD_score_merged) - i + 1, ]]
}
colnames(relevant_tissue) = paste0("Top", 1:ncol(relevant_tissue))

write.table(relevant_tissue, "relevant_tissue.txt", sep = "\t", quote = F)
