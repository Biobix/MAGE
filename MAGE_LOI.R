library("MAGE")

############
#PARAMETERS#
############
chromosomes <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22)
wd_samples <- ""
wd_seq <-  ""
wd_seq_tumor <-  ""
wd_res_tumor <- ""; if(!dir.exists(wd_res)) dir.create(wd_res)

file_samples <- paste(wd_samples, "kidney_samples.txt", sep = "")
file_samples_t <- paste(wd_seq_tumor, "samples_KIRC.txt", sep = "")
file_sequences_t <- paste(wd_seq_tumor, "counts_SNP_chr", sep = "")
# file_sequences_t <- paste(wd_seq, "counts_SNP_tumor_chr", sep = "")

#######
#START#
#######
#READ SAMPLE INFO OF CONTROL DATA AND TUMOUR DATA
sample_info <- read.table(file_samples, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
names(sample_info) <- c("sample_nr", "sample")
samples_all <- sample_info$sample
nr_samples_all <- length(samples_all)

sample_info_t <- read.table(file_samples_t, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
names(sample_info_t) <- c("sample_nr", "sample")
samples_t <- sample_info_t$sample
nr_samples_t <- length(samples_t)

for (chr in chromosomes) {
  results <- read.table(paste(wd_seq, "chr", chr, "_imprinted_genes.txt", sep = ""), header = T, stringsAsFactors = F)
  allelecounts <- read.table(paste(wd_seq, "chr", chr, "_allelecounts_imprinted.txt", sep = ""), header = T, stringsAsFactors = F)

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
      #sample_data_t <- data.table::fread(paste(file_sequences_t, chr, "_tumor.txt", sep=""), header = FALSE, stringsAsFactors = F, sep = "\t", data.table = getOption("datatable.fread.datatable", FALSE))
      sample_data_t <- data.table::fread(paste(file_sequences_t, chr, ".txt", sep=""), header = FALSE, stringsAsFactors = F, sep = "\t", data.table = getOption("datatable.fread.datatable", FALSE))
      names(sample_data_t) <- c("chromosome", "position", "ref_alleles", "variationtype", "dbSNP_ref", "gene", "A", "T", "C", "G", "sample")
      sample_data_t$sample_nr <- unlist(lapply(sample_data_t$sample, function(x) sample_info_t[which(sample_info_t$sample == x), "sample_nr"]))
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


      #CALCULATE MEAN COUNT OF CONTROL AND TUMOUR DATA AND ANALYSE DIFFERENTIAL EXPRESSION
      DE_df <- data.frame()
      for (z in positions_impr_t) {
        logcount_c <- log(data[[z]]$ref_count + data[[z]]$var_count + 0.5)
        logcount_t <- log(data_t[[z]]$ref_count + data_t[[z]]$var_count + 0.5)
        meancount_c <- mean(data[[z]]$ref_count + data[[z]]$var_count)
        meancount_t <- mean(data_t[[z]]$ref_count + data_t[[z]]$var_count)

        DE_z <- t.test(logcount_c, logcount_t)$p.value
        log_fc <- log2(meancount_t / meancount_c)

        DE_df_x <- data.frame("position" = z, "DE_p" = DE_z, "logFC" = log_fc, "mean_c" = meancount_c, "mean_t" = meancount_t)
        DE_df <- rbind(DE_df, DE_df_x)
      }

      results_DE <- merge(results, merge(p_LOI_df, DE_df, by = "position"), by = "position")

      #USE HARMONIC MEAN TO COMBINE P-VALUES OF SNPS PER GENE
      p_total <- hash::hash(); p_total_df <- data.frame()
      for (gene in impr_genes) {
        pos_gene <- as.character(results[grep(gene, results$gene), "position"])
        pos_gene <- pos_gene[which(pos_gene %in% positions_impr_t)]

        if (length(pos_gene) > 0) {
          ##COMBINE P-VALUES LOI PER GENE
          p_values_LOI <- hash::values(p_LOI[pos_gene])
          p_total[[gene]] <- MAGE::combine_p_gene(p_values_LOI)

          p_values_DE <- results_DE[which(results_DE$position %in% pos_gene), "DE_p"]
          p_total_DE <- MAGE::combine_p_gene(p_values_DE)
          logFC_total_DE <- mean(results_DE[which(results_DE$position %in% pos_gene), "logFC"])

          p_total_df <- rbind(p_total_df, data.frame("position" = gene, "p_LOI" = p_total[[gene]], "p_DE" = p_total_DE,"logFC" = logFC_total_DE))
        }
      }

      f2 <- file(paste(wd_res_tumor, "LOI_DE_SNPs_chr", chr, ".txt", sep = ""), "w")
      out <- write.table(results_DE, f2, quote = F, sep = "\t")
      close(f2)

      f2 <- file(paste(wd_res_tumor, "LOI_DE_gene_chr", chr, ".txt", sep = ""), "w")
      out <- write.table(p_total_df, f2, quote = F, sep = "\t")
      close(f2)
      
      for (z in positions_impr_t) {
        MAGE::plot_histo(data[[z]]$ref_count, data[[z]]$var_count, wd_res = wd_res_tumor, chr = chr, position = z, gene  = paste(data[[z]]$gene[1], "control", sep = "_"))
        MAGE::plot_histo(data_t[[z]]$ref_count, data_t[[z]]$var_count, wd_res = wd_res_tumor, chr = chr, position = z, gene  = paste(data[[z]]$gene[1], "tumour", sep = "_"))
      }
    }
  }
}


