library("MAGE")

############
#PARAMETERS#
############
chromosomes <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22)

wd_samples_control <- "" #working directory with control sample file
wd_seq_control <- "" #working directory with control count files
wd_res_control <- "" #working directory with imprinting results files

wd_samples_tumor <- "" #working directory with case sample file
wd_seq_tumor <- "" #working directory with case count files
wd_res_tumor <- ""; if(!dir.exists(wd_res_tumor)) dir.create(wd_res_tumor) #working directory for differential imprinting results

file_samples_control <- paste(wd_samples_control, "kidney_samples.txt", sep = "")
file_sequences_control <- paste(wd_seq_control, "counts_SNP_chr", sep = "")
file_samples_tumor <- paste(wd_seq_tumor, "samples_KIRC.txt", sep = "")
file_sequences_tumor <- paste(wd_seq_tumor, "counts_SNP_chr", sep = "")
# file_sequences_t <- paste(wd_seq, "counts_SNP_tumor_chr", sep = "")

pA_estimate_filt <- 0.1
pA_filt <- 0.15
cov_filt <- 4
nr_samples_filt <- 30
SE_filt <- 0.035
sym_filt <- 0.05
GOF_filt <- 0.8
impr_filt <- 0.6
median_filt <- 0.8
p_filt <- 0.05
plot_coverage <- 40
plot_HWE <- FALSE

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
  #READ COUNT FILES
  sample_data <- data.table::fread(paste(file_sequences_control, chr, ".txt", sep=""), header = FALSE, stringsAsFactors = F, sep = "\t", data.table = getOption("datatable.fread.datatable", FALSE))
  names(sample_data) <- c("chromosome", "position", "ref_alleles", "variationtype", "dbSNP_ref", "gene", "A", "T", "C", "G", "sample")
  sample_data$sample_nr <- unlist(lapply(sample_data$sample, function(x) sample_info_control[which(sample_info_control$sample == x), "sample_nr"]))
  sample_data$position <- as.character(sample_data$position)

  #DETERMINE AMOUNT OF LOCI
  positions <- unique(sample_data$position)
  loci <- length(positions)

  #CREATE HASH WITH DATA FRAME FOR EVERY POSITION
  data <- hash::hash()
  for (i in positions) {
    data[[ i ]] <- sample_data[which(sample_data$position==i), ]
  }

  #DELETE INDELS AND SUBSTITUTIONS
  for (z in positions) {
    if (data[[z]]$variationtype[1] != "SNV" || grepl('\\.', strsplit(data[[z]]$ref_alleles[1], "/"))) {
      data[[z]] <- NULL
    }
  }
  positions <- hash::keys(data); loci <- length(positions)

  #DETERMINE STANDARD ALLELES AND REFERENCE/VARIANT COUNTS
  for (z in positions) {
    data[[z]] <- MAGE::standard_alleles(data[[z]])
  }

  #PRIOR FILTERING OF LOCI ON ALLELE FREQUENCY, COVERAGE, NUMBER OF SAMPLES AND BAYES FILTERING. ALSO NON-STANDARD ALLELIC COUNTS ARE REMOVED.
  for (z in positions) {
    data[[z]] <- MAGE::prior_filter(data[[z]], prior_allelefreq = pA_estimate_filt, coverage_filter = cov_filt, samples_filter = nr_samples_filt)
  }
  positions <- hash::keys(data); loci <- length(positions)

  #ESIMATE PARAMETERS WITH SeqEM (SEQUENCING ERROR RATE, ALLELE FREQUENCY AND INBREEDING), CALCUTE SYMMETRY GOF PER SNP AND FILTER SNPS ON ALLELE FREQUENCY, SEQUENCING ERROR RATE AND SYMMETRY
  SE_all <- NULL; F_all <- NULL
  for (z in positions) {
    parameters_results <- MAGE::estimate_parameters(data[[z]]$ref_count, data[[z]]$var_count)
    data[[z]]$allelefreq <- parameters_results$allelefreq
    data[[z]]$est_SE <- parameters_results$SE
    data[[z]]$est_inbr <- parameters_results$F_inbr
    
    data[[z]]$sym <- MAGE::symmetry_gof(data[[z]]$ref_count, data[[z]]$var_count, parameters_results$allelefreq)

    if (parameters_results$allelefreq <= pA_filt || parameters_results$allelefreq >= (1 - pA_filt) || parameters_results$SE > SE_filt) {
      data[[z]] <- NULL
    } else {
      SE_all <- c(SE_all, parameters_results$SE)
      F_all <- c(F_all, parameters_results$F_inbr)
      
      if (data[[z]]$sym[1] <= sym_filt) {
        data[[z]] <- NULL
      }
    }
  }
  f_inb_chr <- median(F_all)
  SE <- median(SE_all)
  positions <- hash::keys(data); loci <- length(positions)

  #PERFORM IMPRINTING ANALYSIS PER SNP AND CREATE RESULTS DATA FRAME
  results <- data.frame()
  for (z in positions) {
    lrt_results <- MAGE::lrt_i(data[[z]]$ref_count, data[[z]]$var_count, allelefreq = data[[z]]$allelefreq[1], SE = SE, inbr = f_inb_chr)
    med_imp <- MAGE::median_imprinting(data[[z]]$ref_count, data[[z]]$var_count, allelefreq = data[[z]]$allelefreq[1], inbr = f_inb_chr)
    results_z <- data.frame("position" = z, "gene" = data[[z]]$gene[1], "LRT" = lrt_results$LRT, "p" = lrt_results$p_value, "estimated.i" = lrt_results$est_i, "allele.frequency" = data[[z]]$allelefreq[1], "dbSNP" = data[[z]]$dbSNP_ref[1], "reference" = data[[z]]$ref[1], "variant" = data[[z]]$var[1], "est_SE" = data[[z]]$est_SE[1], "coverage" = data[[z]]$coverage[1], "nr_samples" = nrow(data[[z]]), "GOF" = lrt_results$GOF_likelihood, "symmetry" = data[[z]]$sym[1], "med_impr" = med_imp, "est_inbreeding" = data[[z]]$est_inbr[1], "tot_inbreeding" = f_inb_chr, stringsAsFactors = FALSE)
    results <- rbind(results, results_z)
  }
  
  f2 <- file(paste(wd_res, "SE_inbr_chr", chr, ".txt", sep = ""), "w")
  out <- write.table(data.frame("parameter" = c("SE", "inbr"), "value" = c(SE, f_inb_chr)), f2, quote = F, sep = "\t")
  close(f2)
  
  #FINAL FILTERING OF INTERESTING AND IMPRINTED SNPS AND WRITE RESULTS FILES
  results <- MAGE::final_filter(chr, data, results, wd_res_control, gof_filt = GOF_filt, med_impr_filt = median_filt, i_filt = impr_filt, file_all = TRUE, file_impr = TRUE, file_all_counts = FALSE, file_impr_counts = TRUE)
  
  positions_impr <- as.character(results$position)
  
  #PLOT SIGNIFICANTLY IMPRINTED GENES#
  if (length(positions_impr) > 0) {
    for (z in positions_impr) {
      MAGE::plot_imprinting(data[[z]]$ref_count, data[[z]]$var_count, allelefreq = results[which(results$position==z), "allele.frequency"], impr = results[which(results$position==z), "estimated.i"], SE = SE, wd_res = wd_res_control, chr = chr, position = z, gene  = results[which(results$position==z), "gene"], inbr = f_inb_chr, coverage = plot_coverage, plot_hwe = plot_HWE)
    }
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


