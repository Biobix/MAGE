library("MAGE")

############
#PARAMETERS#
############
chromosomes <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22)

wd_samples_control <- "" #working directory with control sample file
wd_res_control <- "" #working directory with imprinting results files

wd_samples_tumor <- "" #working directory with case sample file
wd_seq_tumor <- "" #working directory with case count files
wd_res_tumor <- ""; if(!dir.exists(wd_res_tumor)) dir.create(wd_res_tumor) #working directory for differential imprinting results

file_samples_control <- paste(wd_samples_control, "kidney_samples.txt", sep = "")
file_samples_tumor <- paste(wd_seq_tumor, "samples_KIRC.txt", sep = "")
file_sequences_tumor <- paste(wd_seq_tumor, "counts_SNP_chr", sep = "")

#######
#START#
#######
#READ SAMPLE INFO OF CONTROL DATA AND TUMOUR DATA
sample_info_control <- read.table(file_samples_control, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
names(sample_info_control) <- c("sample_nr", "sample")
samples_control <- sample_info_control$sample
nr_samples_control <- length(samples_control)

sample_info_tumor <- read.table(file_samples_tumor, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
names(sample_info_tumor) <- c("sample_nr", "sample")
samples_tumor <- sample_info_tumor$sample
nr_samples_tumor <- length(samples_tumor)

for (chr in chromosomes) {
  results <- read.table(paste(wd_res_control, "chr", chr, "_imprinted_genes.txt", sep = ""), header = T, stringsAsFactors = F)
  allelecounts <- read.table(paste(wd_res_control, "chr", chr, "_allelecounts_imprinted.txt", sep = ""), header = T, stringsAsFactors = F)

  if (nrow(results) > 0) {
    positions_impr <- as.character(results$position)
    data <- hash::hash()
    for (z in positions_impr) {
      ref_allele <- results[which(results$position == as.numeric(z)), "reference"]
      var_allele <- results[which(results$position == as.numeric(z)), "variant"]

      if(ref_allele == TRUE) {ref_allele <- "T"; results[which(results$position == as.numeric(z)), "reference"] <- "T"}
      if(var_allele == TRUE) {var_allele <- "T"; results[which(results$position == as.numeric(z)), "variant"] <- "T"}

      data[[z]] <- data.frame("chromosome" = chr,  "position" = z, "ref_alleles" = paste(ref_allele, var_allele, sep = "/"),
                              "variationtype" = "SNV", "dbSNP_ref" = results[which(results$position == z), "dbSNP"],
                              "gene" = results[which(results$position == z), "gene"],
                              "sample" = allelecounts[which(allelecounts$position == as.numeric(z)), "sample_id"],
                              "sample_nr" = allelecounts[which(allelecounts$position == as.numeric(z)), "sample_nr"],
                              "ref_count" = allelecounts[which(allelecounts$position == as.numeric(z)), "ref"],
                              "var_count" = allelecounts[which(allelecounts$position == as.numeric(z)), "var"],
                              "ref" = ref_allele,
                              "var" = var_allele, stringsAsFactors = F)
    }

    ####################
    #ANALYSE TUMOR DATA#
    ####################
    if (length(positions_impr) > 0) {
      #READ COUNT FILES
      #sample_data_t <- data.table::fread(paste(file_sequences_tumor, chr, "_tumor.txt", sep=""), header = FALSE, stringsAsFactors = F, sep = "\t", data.table = getOption("datatable.fread.datatable", FALSE))
      sample_data_t <- data.table::fread(paste(file_sequences_tumor, chr, ".txt", sep=""), header = FALSE, stringsAsFactors = F, sep = "\t", data.table = getOption("datatable.fread.datatable", FALSE))
      names(sample_data_t) <- c("chromosome", "position", "ref_alleles", "variationtype", "dbSNP_ref", "gene", "A", "T", "C", "G", "sample")
      sample_data_t$sample_nr <- unlist(lapply(sample_data_t$sample, function(x) sample_info_tumor[which(sample_info_tumor$sample == x), "sample_nr"]))
      sample_data_t$position <- as.character(sample_data_t$position)

      #ONLY RETAIN IMPRINTED POSITIONS ALSO PRESENT IN TUMOUR DATA
      positions_impr_t <- positions_impr[which(positions_impr %in% unique(sample_data_t$position))]
      #USE SAME REFERENCE AND VARIANT ALLELE AS IN CONTROL DATA AND DETERMINE REFERENCE/VARIANT COUNTS
      data_t <- hash::hash()
      for (z in positions_impr_t) {
        data_t[[z]] <- sample_data_t[which(sample_data_t$position == z), ]
        data_t[[z]]$ref_alleles <- data[[z]]$ref_alleles[1]
        data_t[[z]] <- MAGE::standard_alleles(data_t[[z]])

        if (data_t[[z]]$ref[1] == data[[z]]$var[1] && data_t[[z]]$var[1] == data[[z]]$ref[1]) {
          data_t[[z]]$var <- data[[z]]$var[1]; data_t[[z]]$ref <- data[[z]]$ref[1]
          ref_c <- data_t[[z]]$var_count; var_c <- data_t[[z]]$ref_count
          data_t[[z]]$ref_count <- ref_c; data_t[[z]]$var_count <- var_c
        } else if (data_t[[z]]$ref[1] != data[[z]]$ref[1]) {
          print(z)
        }

        if (any(data_t[[z]]$ref_count + data_t[[z]]$var_count == 0)) data_t[[z]] <- data_t[[z]][- which(data_t[[z]]$ref_count + data_t[[z]]$var_count == 0), ]

      }

      #CHECK FOR LOI BETWEEN CASES AND CONTROLS WITH THE LOGISTIC REGRESSION FUNCTION
      p_LOI <- hash::hash(); p_LOI_df <- data.frame()
      for (z in positions_impr_t) {
        p_LOI[[z]] <- MAGE::binomial_logistic(data[[z]]$ref_count, data[[z]]$var_count, data_t[[z]]$ref_count, data_t[[z]]$var_count)$p.value
        p_LOI_df <- rbind(p_LOI_df, data.frame("position" = z, "p_LOI" = p_LOI[[z]]))
      }

      impr_genes <- unique(unlist(strsplit(results$gene, ",")))

      results <- merge(results, p_LOI_df, by = "position")

      #USE HARMONIC MEAN TO COMBINE P-VALUES OF SNPS PER GENE
      p_total <- hash::hash(); p_total_df <- data.frame()
      for (gene in impr_genes) {
        pos_gene <- as.character(results[grep(gene, results$gene), "position"])
        pos_gene <- pos_gene[which(pos_gene %in% positions_impr_t)]

        if (length(pos_gene) > 0) {
          ##COMBINE P-VALUES LOI PER GENE
          p_values_LOI <- hash::values(p_LOI[pos_gene])
          p_total[[gene]] <- MAGE::combine_p_gene(p_values_LOI)

          p_total_df <- rbind(p_total_df, data.frame("position" = gene, "p_LOI" = p_total[[gene]]))
        }
      }

      f2 <- file(paste(wd_res_tumor, "LOI_SNPs_chr", chr, ".txt", sep = ""), "w")
      out <- write.table(results, f2, quote = F, sep = "\t")
      close(f2)

      f2 <- file(paste(wd_res_tumor, "LOI_gene_chr", chr, ".txt", sep = ""), "w")
      out <- write.table(p_total_df, f2, quote = F, sep = "\t")
      close(f2)
      
      for (z in positions_impr_t) {
        MAGE::plot_histo(data[[z]]$ref_count, data[[z]]$var_count, wd_res = wd_res_tumor, chr = chr, position = z, gene  = paste(data[[z]]$gene[1], "control", sep = "_"))
        MAGE::plot_histo(data_t[[z]]$ref_count, data_t[[z]]$var_count, wd_res = wd_res_tumor, chr = chr, position = z, gene  = paste(data[[z]]$gene[1], "tumour", sep = "_"))
      }
    }
  }
}