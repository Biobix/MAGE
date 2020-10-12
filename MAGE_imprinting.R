library("MAGE")

chromosomes <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22)
wd_samples <- "" #working directory with sample file
wd_seq <-  "" #working directory with  count files
wd_res <- ""; if(!dir.exists(wd_res)) dir.create(wd_res) #working directory for  imprinting results

file_samples <- paste(wd_samples, "kidney_samples.txt", sep = "")
file_sequences <- paste(wd_seq, "counts_SNP_chr", sep = "")

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
sample_info <- read.table(file_samples, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
names(sample_info) <- c("sample_nr", "sample")
samples_all <- sample_info$sample
nr_samples_all <- length(samples_all)

for (chr in chromosomes) {
  #READ COUNT FILES
  sample_data <- data.table::fread(paste(file_sequences, chr, ".txt", sep=""), header = FALSE, stringsAsFactors = F, sep = "\t", data.table = getOption("datatable.fread.datatable", FALSE))
  names(sample_data) <- c("chromosome", "position", "ref_alleles", "variationtype", "dbSNP_ref", "gene", "A", "T", "C", "G", "sample")
  sample_data$sample_nr <- unlist(lapply(sample_data$sample, function(x) sample_info[which(sample_info$sample == x), "sample_nr"]))
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
  
  #PRIOR FILTER OF LOCI
  for (z in positions) {
    data[[z]] <- MAGE::prior_filter(data[[z]], prior_allelefreq = pA_estimate_filt, coverage_filter = cov_filt, samples_filter = nr_samples_filt)
  }
  positions <- hash::keys(data); loci <- length(positions)
  
  #ESIMATE PARAMETERS WITH SeqEM (SEQUENCING ERROR RATE, ALLELE FREQUENCY AND INBREEDING), CALCUTE SYMMETRY GOF PER SNP AND FILTER SNPS
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
  results <- MAGE::final_filter(chr, data, results, wd_res, gof_filt = GOF_filt, med_impr_filt = median_filt, i_filt = impr_filt, file_all = TRUE, file_impr = TRUE, file_all_counts = FALSE, file_impr_counts = TRUE)
  
  positions_impr <- as.character(results$position)
  
  #PLOT SIGNIFICANTLY IMPRINTED GENES#
  if (length(positions_impr) > 0) {
    for (z in positions_impr) {
      MAGE::plot_imprinting(data[[z]]$ref_count, data[[z]]$var_count, allelefreq = results[which(results$position==z), "allele.frequency"], impr = results[which(results$position==z), "estimated.i"], SE = SE, wd_res = wd_res, chr = chr, position = z, gene  = results[which(results$position==z), "gene"], inbr = f_inb_chr, coverage = plot_coverage, plot_hwe = plot_HWE)
    }
  }
}


