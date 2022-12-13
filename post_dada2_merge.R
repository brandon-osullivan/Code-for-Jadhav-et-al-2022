#
# Filename: post_dada2_merge.R
#
# Purpose: Merge dada2 frequency tables.
#		For local windows cmd or terminal execution with Rscript, interactive cluster session with Rscript, or cluster batch submission with Rscript
#
# Created: 28-Jan-2021
#
# Revision: A1
#
# History:
#

#Merge dada2 frequency tables
#Run sbsearch.exe to get .tax files
#gather the csv, fasta, and tax for each run you want to merge into their own directory within another directory
#example C:/runs_to_merge/run1/run1.csv C:/runs_to_merge/run2/run2.csv C:/runs_to_merge/run3/run3.csv
#give the directory C:/runs_to_merge/ as an argument, with NO OTHER DIRECTORIES inside it, besides the runs themselves

biocm_pkg_list <- c("Biostrings","Rhtslib","Rsamtools","GenomicAlignments","hwriter","ShortRead","RcppParallel","plyr","reshape2","dada2","Rhdf5lib","pixmap","rhdf5","sp","ade4","permute","vegan","multtest","igraph","iterators","foreach","biomformat","ape","phyloseq","ggplot2","gridExtra","RColorBrewer","docopt")
for (i in 1:length(biocm_pkg_list)) {
	library(biocm_pkg_list[i],quietly=TRUE,character.only=TRUE)
}

"Usage:
	script.R (-n RUNNAME) (-i INDIR)

Options:
	-n RUNNAME	appends to filenames
	-i INDIR	directory with csvs and asv fastas" -> doc

opts <- docopt(doc,help=TRUE)
run_name <- opts$n
setwd(opts$i)
if (.Platform$OS.type == "windows") {
	path <- gsub("/", "\\\\", getwd())
} else if (.Platform$OS.type == "unix") {
	path <- getwd()
}
#check user supplied arguments are correct format etc
if (!(dir.exists(opts$i))) {
	stop("Directory does not exist")
}
if(grepl(pattern="\\s", run_name, perl=TRUE)) {
	stop("Do not use spaces in the run name")
}

#makes a list of all the directories in the given directory
data_dir_list <- list.dirs(path, full.names=TRUE)
print(data_dir_list)

#get the filepaths for all the inputs
csv_in_list <- vector()
fasta_in_list <- vector()
tax_in_list <- vector()
print(data_dir_list[-1])
for (indir in data_dir_list[-1]) {
	csv_in_list <- append(csv_in_list, list.files(indir, pattern=".csv", full.names=TRUE))
	fasta_in_list <- append(fasta_in_list, list.files(indir, pattern=".fasta", full.names=TRUE))
	tax_in_list <- append(tax_in_list, list.files(indir, pattern=".tax", full.names=TRUE))
}
print(csv_in_list)
print(fasta_in_list)
print(tax_in_list)

#load the inputs
csv_load <- list()
fasta_load <- list()
tax_load <- list()
for (i in 1:(length(data_dir_list)-1)) {
	csv_load[[i]] <- read.csv(csv_in_list[i], check.names=FALSE, row.names=1)
	fasta_load[[i]] <- getSequences(fasta_in_list[i])
	tax_load[[i]] <- read.table(tax_in_list[i])[,2]
	names(tax_load[[i]]) <- read.table(tax_in_list[i])[,1]
}

#harmonize csv columns for merging
csv_premerge <- list()
for (i in 1:length(csv_load)) {
	csv_premerge[[i]] <- csv_load[[i]]
	for (j in 1:length(csv_load)) {
		if(i != j) {
			for (colnm in colnames(csv_load[[j]][!(colnames(csv_load[[j]]) %in% colnames(csv_premerge[[i]]))])) {
				csv_premerge[[i]][[colnm]] <- 0
			}
		}
	}
}

csv_merge <- csv_premerge[[1]]
#merge csvs
for (i in 1:(length(csv_premerge)-1)) {
	csv_merge <- rbind(csv_merge, csv_premerge[[i+1]])
}

#concatenate fastas and tax
fasta_premerge <- fasta_load[[1]]
for (i in 1:(length(fasta_load)-1)) {
	fasta_premerge <- c(fasta_premerge, fasta_load[[i+1]])
}

tax_premerge <- tax_load[[1]]
for (i in 1:(length(tax_load)-1)) {
	tax_premerge <- c(tax_premerge, tax_load[[i+1]])
}

#merge into one table
#make heatmap to get order
st_tot_heatmap <- t(as.matrix(csv_merge))
png(file.path(path, paste(run_name, "_asv_heatmap_merged.png", sep="")), width = 14, height = 8, units = "in", res = 450)
st_tot_heatmap_dendro <- heatmap(st_tot_heatmap[], Colv = NA, main="ASV heatmap", distfun=function (y) dist(y,method = "canb"), col = hcl.colors(64, palette="Blue-Red 3"), scale = "column", cexRow = 0.5, cexCol = 0.5, keep.dendro=TRUE)
dev.off()
#put table in order
csv_to_write <- t(st_tot_heatmap)[,order.dendrogram(st_tot_heatmap_dendro$Rowv)]

new_mat_1 <- matrix(nrow=dim(csv_to_write)[1]+1, ncol=dim(csv_to_write)[2])
new_mat_1[-1,] <- csv_to_write

seq_vector <- vector()
for (colnm in colnames(csv_to_write)) {
	seq_vector <- append(seq_vector, unname(fasta_premerge[colnm]))
}

new_mat_1[1,] <- seq_vector

new_mat_2 <- matrix(nrow=dim(new_mat_1)[1]+1, ncol=dim(new_mat_1)[2])
new_mat_2[-1,] <- new_mat_1

new_mat_2[1,] <- tax_premerge[colnames(csv_to_write)]
rownames(new_mat_2) <- append(c("taxonomy","sequence"),rownames(csv_to_write))
colnames(new_mat_2) <- colnames(csv_to_write)

write.csv(t(new_mat_2), file.path(path, "RosettaReport_merged.csv"))

write.csv(t(csv_merge), file.path(path, "sequence_table_merged.csv"))
csv_merge_biom <- make_biom(t(csv_merge))
write_biom(csv_merge_biom, file.path(path, "sequence_table_merged.biom"))

fasta_merge <- fasta_premerge[unique(names(fasta_premerge))]
jal_ready_out <- file(description = file.path(path, "ASVs_merged.fasta"), open = "w")
for (i in 1:length(fasta_merge)) {
	writeLines(paste(">", names(fasta_merge[i]), sep=""), jal_ready_out)
	writeLines(unname(fasta_merge[i]), jal_ready_out)
}
close(jal_ready_out)

tax_merge <- tax_premerge[unique(names(tax_premerge))]
jal_ready_out <- file(description = file.path(path, "taxonomy_merged.tax"), open = "w")
for (i in 1:length(tax_merge)) {
	writeLines(paste(names(tax_merge[i]), unname(tax_merge[i]), sep="\t"), jal_ready_out)
}
close(jal_ready_out)

metadata_starter <- file(description = file.path(path, "metadata_starter.tsv"), open = "w")
writeLines("sample-id\tmetadata_name", metadata_starter)
writeLines("#q2:types\tcategorical", metadata_starter)
for (i in rownames(csv_merge)) {
	writeLines(i, metadata_starter)
}
close(metadata_starter)

save.image(file.path(path, paste(run_name, ".RData", sep="")))
