### SV_comparison_program analysis 
### Author: Rodrigo López
### Date: 14/02/2024


rm(list=ls())


########################################################################################################################################################
# PACKAGE REQUIREMENTS #
########################


required_packages <- c("optparse","data.table", "bedr", "dplyr", "tidyverse")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
  install.packages(missing_packages)
}

library(optparse)
library(data.table)
library(bedr)
library(dplyr)
library(tidyverse)


########################################################################################################################################################
# Arguments #
#############


option_list=list(
  make_option(c('-i','--inputdir'),type="character", help="Input directory."),
  make_option(c('-o','--outputdir'),type="character", help="Output directory.")
)
  
opt_parser=OptionParser(option_list = option_list)
opt=parse_args(opt_parser) #list of the args


# opt=list()
# opt$inputdir <- "/home/rodrigo/tblab/rodrigo/data/VCF"
# opt$outputdir <- "/mnt/tblab/LongSeq/data/"

inputdir <- opt$inputdir
outputdir <- opt$outputdir

setwd(inputdir)


########################################################################################################################################################
#   Input loading  #
####################


toAnnotate_files = list.files(inputdir, pattern = ".vcf")
all_progrmas_df <- data.frame(stringsAsFactors = FALSE)

for (toAnnotate_file in toAnnotate_files) {
  print(paste("Loading file:", toAnnotate_file))
  sample_name <- sub("^(.*?)\\..*", "\\1", tools::file_path_sans_ext(basename(toAnnotate_file)))
  toAnnotate_df <- fread(toAnnotate_file, sep = "\t", header = TRUE, skip = "#CHROM")
  colnames(toAnnotate_df)[10] <- "SAMPLE"
  toAnnotate_df$SAMPLE_ID <- sample_name
  toAnnotate_df$N_SAMPLE <- 1
  toAnnotate_df$SAMPLE_SIMILARITY_PERCENTAGE <- 0
  if (!any(grepl("[a-zA-Z]", toAnnotate_df$ID))) {
    toAnnotate_df$PROGRAM <- sub(".*SVMETHOD=([^;]+).*", "\\1", toAnnotate_df$INFO)
  } else {
    toAnnotate_df$PROGRAM <- sub("^(\\w+).*", "\\1", toAnnotate_df$ID)
  }
  if (any(grepl("AF=", toAnnotate_df$INFO))) {
    af_values <- gsub(".*AF=([0-9.]+).*", "\\1", toAnnotate_df$INFO)
    toAnnotate_df <- toAnnotate_df[af_values >= 0.3 | is.na(af_values), ]
  }
  toAnnotate_df$N_PROGRAM <- 1
  toAnnotate_df$PROGRAM_SIMILARITY_PERCENTAGE <- 0
  toAnnotate_df$N_NUCLEOTIDE_COINCIDENCE <- 0
  toAnnotate_df <- toAnnotate_df[toAnnotate_df$QUAL > 10 | toAnnotate_df$QUAL == ".", ]
  all_progrmas_df <- rbind(all_progrmas_df, toAnnotate_df)
}

colnames(all_progrmas_df)[colnames(all_progrmas_df) == "#CHROM"] <- "CHROM"


########################################################################################################################################################
# DEL #
#######


SV_DEL <- data.frame()
print("Merging DEL")
for (i in 1:nrow(all_progrmas_df)) {
  current_row <- all_progrmas_df[i, ]
  if (grepl("SVTYPE=.*DEL", current_row$INFO)) {
    end_value <- as.numeric(sub(".*END=([0-9]+).*", "\\1", current_row$INFO, perl = TRUE))
    current_row$END <- end_value
    if (grepl("DR:DV", current_row$FORMAT)) {
      dr_position <- grep("DR", unlist(strsplit(as.character(current_row$FORMAT), ":")))
      dv_position <- grep("DV", unlist(strsplit(as.character(current_row$FORMAT), ":")))
      sample_values <- unlist(strsplit(as.character(current_row$SAMPLE), ":"))
      sample_values[dr_position] <- paste(sample_values[dr_position], sample_values[dv_position], sep = ",")
      sample_values <- sample_values[-dv_position]
      current_row$SAMPLE <- paste(sample_values, collapse = ":")
      current_row$FORMAT <- gsub("DR:DV", "AD", current_row$FORMAT)
    }
    current_row$GT <- unlist(strsplit(as.character(current_row$SAMPLE[1]), ":"))[1]
    SV_DEL <- rbind(SV_DEL, current_row)
  }
}


#In case of similar POS values, it retains maximum END values
DEL_pos <- data.frame()

for (cromosoma in unique(SV_DEL$CHROM)) {
  subset_data <- SV_DEL[SV_DEL$CHROM == cromosoma, ]
  for (pos_value in unique(subset_data$POS)) {
    matching_rows <- subset_data[subset_data$POS == pos_value, ]
    if (nrow(matching_rows) > 1 & length(unique(matching_rows$GT)) == 1) {
      max_similarity <- max(matching_rows$END)
      min_coherence <- min(matching_rows$END)
      if (length(unique(matching_rows$PROGRAM)) > 1) {
        programs_concat <- paste(unique(matching_rows$PROGRAM), collapse = ",")
        program_count <- length(unique(unlist(strsplit(programs_concat, ","))))
        matching_rows$PROGRAM <- programs_concat
        matching_rows$N_PROGRAM <- program_count
        matching_rows$PROGRAM_SIMILARITY_PERCENTAGE <- round((min_coherence / max_similarity) * 100, 3)
      }
      if (length(unique(matching_rows$SAMPLE_ID)) > 1) {
        sample_concat <- paste(unique(matching_rows$SAMPLE_ID), collapse = ",")
        sample_count <- length(unique(unlist(strsplit(sample_concat, ","))))
        matching_rows$SAMPLE_ID <- sample_concat
        matching_rows$N_SAMPLE <- sample_count
        matching_rows$SAMPLE_SIMILARITY_PERCENTAGE <- round((min_coherence / max_similarity) * 100, 3)
      }
      matching_rows$N_NUCLEOTIDE_COINCIDENCE <- min_coherence - matching_rows$POS
      matching_rows$ID <- paste(matching_rows$ID, collapse = ",")
      max_info_row <- matching_rows[which.max(nchar(matching_rows$INFO)), ]
      matching_rows$INFO <- max_info_row$INFO
      max_format_row <- matching_rows[which.max(nchar(matching_rows$FORMAT)), ]
      matching_rows$FORMAT <- max_format_row$FORMAT
      max_sample_row <- matching_rows[which.max(nchar(matching_rows$SAMPLE)), ]
      matching_rows$SAMPLE <- max_sample_row$SAMPLE
      matching_rows$QUAL <- matching_rows$QUAL[which.max(nchar(matching_rows$QUAL))]
      matching_rows$REF <- matching_rows$REF[which.max(nchar(matching_rows$REF))]
      matching_rows$ALT <- matching_rows$ALT[which.max(nchar(matching_rows$ALT))]
      row_max_end <- matching_rows[which.max(matching_rows$END), ]
      DEL_pos <- rbind(DEL_pos, row_max_end)
    } else {
      DEL_pos <- rbind(DEL_pos, matching_rows)
    }
  }
}


#In case of similar END values, it retains minimum POS values
DEL_pos_end <- data.frame()

for (cromosoma in unique(DEL_pos$CHROM)) {
  subset_data <- DEL_pos[DEL_pos$CHROM == cromosoma, ]
  for (end_value in unique(subset_data$END)) {
    matching_rows <- subset_data[subset_data$END == end_value, ]
    if (nrow(matching_rows) > 1 & length(unique(matching_rows$GT)) == 1) {
      max_similarity <- max(matching_rows$POS)
      min_coherence <- min(matching_rows$POS)
      if (length(unique(matching_rows$PROGRAM)) > 1) {
        programs_concat <- paste(unique(matching_rows$PROGRAM), collapse = ",")
        program_count <- length(unique(unlist(strsplit(programs_concat, ","))))
        matching_rows$PROGRAM <- programs_concat
        matching_rows$N_PROGRAM <- program_count
        matching_rows$PROGRAM_SIMILARITY_PERCENTAGE <- round((min_coherence / max_similarity) * 100, 3)
      }
      if (length(unique(matching_rows$SAMPLE_ID)) > 1) {
        sample_concat <- paste(unique(matching_rows$SAMPLE_ID), collapse = ",")
        sample_count <- length(unique(unlist(strsplit(sample_concat, ","))))
        matching_rows$SAMPLE_ID <- sample_concat
        matching_rows$N_SAMPLE <- sample_count
        matching_rows$SAMPLE_SIMILARITY_PERCENTAGE <- round((min_coherence / max_similarity) * 100, 3)
      }
      matching_rows$N_NUCLEOTIDE_COINCIDENCE <- matching_rows$END - max_similarity
      matching_rows$ID <- paste(matching_rows$ID, collapse = ",")
      max_info_row <- matching_rows[which.max(nchar(matching_rows$INFO)), ]
      matching_rows$INFO <- max_info_row$INFO
      max_format_row <- matching_rows[which.max(nchar(matching_rows$FORMAT)), ]
      matching_rows$FORMAT <- max_format_row$FORMAT
      max_sample_row <- matching_rows[which.max(nchar(matching_rows$SAMPLE)), ]
      matching_rows$SAMPLE <- max_sample_row$SAMPLE
      matching_rows$QUAL <- matching_rows$QUAL[which.max(nchar(matching_rows$QUAL))]
      matching_rows$REF <- matching_rows$REF[which.max(nchar(matching_rows$REF))]
      matching_rows$ALT <- matching_rows$ALT[which.max(nchar(matching_rows$ALT))]
      row_min_pos <- matching_rows[which.min(matching_rows$POS), ]
      DEL_pos_end <- rbind(DEL_pos_end, row_min_pos)
    } else {
      DEL_pos_end <- rbind(DEL_pos_end, matching_rows)
    }
  }
}
SV_DEL <- DEL_pos_end


#Filtered coincident regions (new POS value = minimum, new END value = maximum)
DEL_coincidence <- data.frame()

for (cromosoma in unique(SV_DEL$CHROM)) {
  subset_data <- SV_DEL[SV_DEL$CHROM == cromosoma, ]
  overlaps <- data.frame()
  for (i in 1:nrow(subset_data)) {
    current_row <- subset_data[i, ]
    overlapping_rows <- subset_data[subset_data$POS <= current_row$END & subset_data$END >= current_row$POS, ]
    if (nrow(overlapping_rows) > 1 & length(unique(overlapping_rows$GT)) == 1) {
      min_pos <- min(overlapping_rows$POS)
      max_end <- max(overlapping_rows$END)
      total_length <- max_end - min_pos + 1 #+1 evita que sea 0
      overlap_length <- min(current_row$END, max_end) - max(current_row$POS, min_pos) + 1
      overlap_percentage <- round((overlap_length / total_length) * 100, 3)
      current_row$POS <- min_pos
      current_row$END <- max_end
      current_row$N_NUCLEOTIDE_COINCIDENCE <- min(overlapping_rows$END) - max(overlapping_rows$POS)
      current_row$ID <- paste(overlapping_rows$ID, collapse = ",")
      max_info_row <- overlapping_rows[which.max(nchar(overlapping_rows$INFO)), ]
      current_row$INFO <- max_info_row$INFO
      max_format_row <- overlapping_rows[which.max(nchar(overlapping_rows$FORMAT)), ]
      current_row$FORMAT <- max_format_row$FORMAT
      max_sample_row <- overlapping_rows[which.max(nchar(overlapping_rows$SAMPLE)), ]
      current_row$SAMPLE <- max_sample_row$SAMPLE
      current_row$QUAL <- overlapping_rows$QUAL[which.max(nchar(overlapping_rows$QUAL))]
      current_row$REF <- overlapping_rows$REF[which.max(nchar(overlapping_rows$REF))]
      current_row$ALT <- overlapping_rows$ALT[which.max(nchar(overlapping_rows$ALT))]
      svlen_index <- grep("SVLEN=", current_row$INFO)
      if (length(svlen_index) > 0) {
        svlen_new <- max_end - min_pos
        current_row$INFO[svlen_index] <- gsub("SVLEN=[-0-9]+", paste("SVLEN=-", svlen_new, sep = ""), current_row$INFO[svlen_index])
      }
      end_index <- grep("END=", current_row$INFO)
      if (length(end_index) > 0) {
        current_row$INFO[end_index] <- gsub("END=[0-9]+", paste("END=", max_end, sep = ""), current_row$INFO[end_index])
      }
      if (length(unique(overlapping_rows$PROGRAM)) > 1) {
        programs_concat <- paste(unique(overlapping_rows$PROGRAM), collapse = ",")
        program_count <- length(unique(unlist(strsplit(programs_concat, ","))))
        current_row$PROGRAM <- programs_concat
        current_row$N_PROGRAM <- program_count
        current_row$PROGRAM_SIMILARITY_PERCENTAGE <- overlap_percentage
      }
      if (length(unique(overlapping_rows$SAMPLE_ID)) > 1) {
        sample_concat <- paste(unique(overlapping_rows$SAMPLE_ID), collapse = ",")
        sample_count <- length(unique(unlist(strsplit(sample_concat, ","))))
        current_row$SAMPLE_ID <- sample_concat
        current_row$N_SAMPLE <- sample_count
        current_row$SAMPLE_SIMILARITY_PERCENTAGE <- overlap_percentage
      }
      non_missing_qual <- overlapping_rows$QUAL[overlapping_rows$QUAL != "."]
      if (length(non_missing_qual) > 0) {
        max_qual <- max(as.numeric(non_missing_qual), na.rm = TRUE)
        current_row$QUAL <- ifelse(is.na(max_qual), ".", as.character(max_qual))
      }
      overlaps <- rbind(overlaps, current_row)
    } else {
      overlaps <- rbind(overlaps, current_row)
    }
  }
  overlaps <- overlaps[order(-as.numeric(sapply(overlaps$INFO, function(x) gsub(".*END=([0-9]+).*", "\\1", x)))), ]
  overlaps <- overlaps[!duplicated(overlaps$POS), ]
  DEL_coincidence <- rbind(DEL_coincidence, overlaps)
}
DEL_coincidence$GT <- NULL
SV_DEL <- DEL_coincidence


########################################################################################################################################################
# INV #
#######


SV_INV <- data.frame()
print("Merging INV")
for (i in 1:nrow(all_progrmas_df)) {
  current_row <- all_progrmas_df[i, ]
  if (grepl("SVTYPE=.*INV", current_row$INFO)) {
    end_value <- as.numeric(sub(".*END=([0-9]+).*", "\\1", current_row$INFO, perl = TRUE))
    current_row$END <- end_value
    if (grepl("DR:DV", current_row$FORMAT)) {
      dr_position <- grep("DR", unlist(strsplit(as.character(current_row$FORMAT), ":")))
      dv_position <- grep("DV", unlist(strsplit(as.character(current_row$FORMAT), ":")))
      sample_values <- unlist(strsplit(as.character(current_row$SAMPLE), ":"))
      sample_values[dr_position] <- paste(sample_values[dr_position], sample_values[dv_position], sep = ",")
      sample_values <- sample_values[-dv_position]
      current_row$SAMPLE <- paste(sample_values, collapse = ":")
      current_row$FORMAT <- gsub("DR:DV", "AD", current_row$FORMAT)
    }
    current_row$GT <- unlist(strsplit(as.character(current_row$SAMPLE[1]), ":"))[1]
    SV_INV <- rbind(SV_INV, current_row)
  }
}

#In case of similar POS values, it retains maximum END values
INV_pos <- data.frame()

for (cromosoma in unique(SV_INV$CHROM)) {
  subset_data <- SV_INV[SV_INV$CHROM == cromosoma, ]
  for (pos_value in unique(subset_data$POS)) {
    matching_rows <- subset_data[subset_data$POS == pos_value, ]
    if (nrow(matching_rows) > 1 & length(unique(matching_rows$GT)) == 1) {
      max_similarity <- max(matching_rows$END)
      min_coherence <- min(matching_rows$END)
      if (length(unique(matching_rows$PROGRAM)) > 1) {
        programs_concat <- paste(unique(matching_rows$PROGRAM), collapse = ",")
        program_count <- length(unique(unlist(strsplit(programs_concat, ","))))
        matching_rows$PROGRAM <- programs_concat
        matching_rows$N_PROGRAM <- program_count
        matching_rows$PROGRAM_SIMILARITY_PERCENTAGE <- round((min_coherence / max_similarity) * 100, 3)
      }
      if (length(unique(matching_rows$SAMPLE_ID)) > 1) {
        sample_concat <- paste(unique(matching_rows$SAMPLE_ID), collapse = ",")
        sample_count <- length(unique(unlist(strsplit(sample_concat, ","))))
        matching_rows$SAMPLE_ID <- sample_concat
        matching_rows$N_SAMPLE <- sample_count
        matching_rows$SAMPLE_SIMILARITY_PERCENTAGE <- round((min_coherence / max_similarity) * 100, 3)
      }
      matching_rows$N_NUCLEOTIDE_COINCIDENCE <- min_coherence - matching_rows$POS
      matching_rows$ID <- paste(matching_rows$ID, collapse = ",")
      max_info_row <- matching_rows[which.max(nchar(matching_rows$INFO)), ]
      matching_rows$INFO <- max_info_row$INFO
      max_format_row <- matching_rows[which.max(nchar(matching_rows$FORMAT)), ]
      matching_rows$FORMAT <- max_format_row$FORMAT
      max_sample_row <- matching_rows[which.max(nchar(matching_rows$SAMPLE)), ]
      matching_rows$SAMPLE <- max_sample_row$SAMPLE
      matching_rows$QUAL <- matching_rows$QUAL[which.max(nchar(matching_rows$QUAL))]
      matching_rows$REF <- matching_rows$REF[which.max(nchar(matching_rows$REF))]
      matching_rows$ALT <- matching_rows$ALT[which.max(nchar(matching_rows$ALT))]
      row_max_end <- matching_rows[which.max(matching_rows$END), ]
      INV_pos <- rbind(INV_pos, row_max_end)
    } else {
      INV_pos <- rbind(INV_pos, matching_rows)
    }
  }
}


#In case of similar END values, it retains minimum POS values
INV_pos_end <- data.frame()

for (cromosoma in unique(INV_pos$CHROM)) {
  subset_data <- INV_pos[INV_pos$CHROM == cromosoma, ]
  for (end_value in unique(subset_data$END)) {
    matching_rows <- subset_data[subset_data$END == end_value, ]
    if (nrow(matching_rows) > 1 & length(unique(matching_rows$GT)) == 1) {
      max_similarity <- max(matching_rows$POS)
      min_coherence <- min(matching_rows$POS)
      if (length(unique(matching_rows$PROGRAM)) > 1) {
        programs_concat <- paste(unique(matching_rows$PROGRAM), collapse = ",")
        program_count <- length(unique(unlist(strsplit(programs_concat, ","))))
        matching_rows$PROGRAM <- programs_concat
        matching_rows$N_PROGRAM <- program_count
        matching_rows$PROGRAM_SIMILARITY_PERCENTAGE <- round((min_coherence / max_similarity) * 100, 3)
      }
      if (length(unique(matching_rows$SAMPLE_ID)) > 1) {
        sample_concat <- paste(unique(matching_rows$SAMPLE_ID), collapse = ",")
        sample_count <- length(unique(unlist(strsplit(sample_concat, ","))))
        matching_rows$SAMPLE_ID <- sample_concat
        matching_rows$N_SAMPLE <- sample_count
        matching_rows$SAMPLE_SIMILARITY_PERCENTAGE <- round((min_coherence / max_similarity) * 100, 3)
      }
      matching_rows$N_NUCLEOTIDE_COINCIDENCE <- matching_rows$END - max_similarity
      matching_rows$ID <- paste(matching_rows$ID, collapse = ",")
      max_info_row <- matching_rows[which.max(nchar(matching_rows$INFO)), ]
      matching_rows$INFO <- max_info_row$INFO
      max_format_row <- matching_rows[which.max(nchar(matching_rows$FORMAT)), ]
      matching_rows$FORMAT <- max_format_row$FORMAT
      max_sample_row <- matching_rows[which.max(nchar(matching_rows$SAMPLE)), ]
      matching_rows$SAMPLE <- max_sample_row$SAMPLE
      matching_rows$QUAL <- matching_rows$QUAL[which.max(nchar(matching_rows$QUAL))]
      matching_rows$REF <- matching_rows$REF[which.max(nchar(matching_rows$REF))]
      matching_rows$ALT <- matching_rows$ALT[which.max(nchar(matching_rows$ALT))]
      row_min_pos <- matching_rows[which.min(matching_rows$POS), ]
      INV_pos_end <- rbind(INV_pos_end, row_min_pos)
    } else {
      INV_pos_end <- rbind(INV_pos_end, matching_rows)
    }
  }
}
SV_INV <- INV_pos_end


#Filtered coincident regions (new POS value = minimum, new END value = maximum)
INV_coincidence <- data.frame()

for (cromosoma in unique(SV_INV$CHROM)) {
  subset_data <- SV_INV[SV_INV$CHROM == cromosoma, ]
  overlaps <- data.frame()
  for (i in 1:nrow(subset_data)) {
    current_row <- subset_data[i, ]
    overlapping_rows <- subset_data[subset_data$POS <= current_row$END & subset_data$END >= current_row$POS, ]
    if (nrow(overlapping_rows) > 1 & length(unique(overlapping_rows$GT)) == 1) {
      min_pos <- min(overlapping_rows$POS)
      max_end <- max(overlapping_rows$END)
      total_length <- max_end - min_pos + 1 #en el caso que la coincidencia sea igual, para que no sea 0 se usa +1
      overlap_length <- min(current_row$END, max_end) - max(current_row$POS, min_pos) + 1
      overlap_percentage <- round((overlap_length / total_length) * 100, 3)
      current_row$POS <- min_pos
      current_row$END <- max_end
      current_row$N_NUCLEOTIDE_COINCIDENCE <- min(overlapping_rows$END) - max(overlapping_rows$POS)
      current_row$ID <- paste(overlapping_rows$ID, collapse = ",")
      max_info_row <- overlapping_rows[which.max(nchar(overlapping_rows$INFO)), ]
      current_row$INFO <- max_info_row$INFO
      max_format_row <- overlapping_rows[which.max(nchar(overlapping_rows$FORMAT)), ]
      current_row$FORMAT <- max_format_row$FORMAT
      max_sample_row <- overlapping_rows[which.max(nchar(overlapping_rows$SAMPLE)), ]
      current_row$SAMPLE <- max_sample_row$SAMPLE
      current_row$QUAL <- overlapping_rows$QUAL[which.max(nchar(overlapping_rows$QUAL))]
      current_row$REF <- overlapping_rows$REF[which.max(nchar(overlapping_rows$REF))]
      current_row$ALT <- overlapping_rows$ALT[which.max(nchar(overlapping_rows$ALT))]
      svlen_index <- grep("SVLEN=", current_row$INFO)
      if (length(svlen_index) > 0) {
        svlen_new <- max_end - min_pos
        current_row$INFO[svlen_index] <- gsub("SVLEN=[-0-9]+", paste("SVLEN=", svlen_new, sep = ""), current_row$INFO[svlen_index])
      }
      end_index <- grep("END=", current_row$INFO)
      if (length(end_index) > 0) {
        current_row$INFO[end_index] <- gsub("END=[0-9]+", paste("END=", max_end, sep = ""), current_row$INFO[end_index])
      }
      if (length(unique(overlapping_rows$PROGRAM)) > 1) {
        programs_concat <- paste(unique(overlapping_rows$PROGRAM), collapse = ",")
        program_count <- length(unique(unlist(strsplit(programs_concat, ","))))
        current_row$PROGRAM <- programs_concat
        current_row$N_PROGRAM <- program_count
        current_row$PROGRAM_SIMILARITY_PERCENTAGE <- overlap_percentage
      }
      if (length(unique(overlapping_rows$SAMPLE_ID)) > 1) {
        sample_concat <- paste(unique(overlapping_rows$SAMPLE_ID), collapse = ",")
        sample_count <- length(unique(unlist(strsplit(sample_concat, ","))))
        current_row$SAMPLE_ID <- sample_concat
        current_row$N_SAMPLE <- sample_count
        current_row$SAMPLE_SIMILARITY_PERCENTAGE <- overlap_percentage
      }
      non_missing_qual <- overlapping_rows$QUAL[overlapping_rows$QUAL != "."]
      if (length(non_missing_qual) > 0) {
        max_qual <- max(as.numeric(non_missing_qual), na.rm = TRUE)
        current_row$QUAL <- ifelse(is.na(max_qual), ".", as.character(max_qual))
      }
      overlaps <- rbind(overlaps, current_row)
    } else {
      overlaps <- rbind(overlaps, current_row)
    }
  }
  overlaps <- overlaps[order(-as.numeric(sapply(overlaps$INFO, function(x) gsub(".*END=([0-9]+).*", "\\1", x)))), ]
  overlaps <- overlaps[!duplicated(overlaps$POS), ]
  INV_coincidence <- rbind(INV_coincidence, overlaps)
}
INV_coincidence$GT <- NULL
SV_INV <- INV_coincidence


########################################################################################################################################################
# DUP #
#######


SV_DUP <- data.frame()
print("Merging DUP")
for (i in 1:nrow(all_progrmas_df)) {
  current_row <- all_progrmas_df[i, ]
  if (grepl("SVTYPE=.*DUP", current_row$INFO)) {
    end_value <- as.numeric(sub(".*END=([0-9]+).*", "\\1", current_row$INFO, perl = TRUE))
    current_row$END <- end_value
    if (grepl("DR:DV", current_row$FORMAT)) {
      dr_position <- grep("DR", unlist(strsplit(as.character(current_row$FORMAT), ":")))
      dv_position <- grep("DV", unlist(strsplit(as.character(current_row$FORMAT), ":")))
      sample_values <- unlist(strsplit(as.character(current_row$SAMPLE), ":"))
      sample_values[dr_position] <- paste(sample_values[dr_position], sample_values[dv_position], sep = ",")
      sample_values <- sample_values[-dv_position]
      current_row$SAMPLE <- paste(sample_values, collapse = ":")
      current_row$FORMAT <- gsub("DR:DV", "AD", current_row$FORMAT)
    }
    current_row$GT <- unlist(strsplit(as.character(current_row$SAMPLE[1]), ":"))[1]
    SV_DUP <- rbind(SV_DUP, current_row)
  }
}


#In case of similar POS values, it retains maximum END values
DUP_pos <- data.frame()

for (cromosoma in unique(SV_DUP$CHROM)) {
  subset_data <- SV_DUP[SV_DUP$CHROM == cromosoma, ]
  for (pos_value in unique(subset_data$POS)) {
    matching_rows <- subset_data[subset_data$POS == pos_value, ]
    if (nrow(matching_rows) > 1 & length(unique(matching_rows$GT)) == 1) {
      max_similarity <- max(matching_rows$END)
      min_coherence <- min(matching_rows$END)
      if (length(unique(matching_rows$PROGRAM)) > 1) {
        programs_concat <- paste(unique(matching_rows$PROGRAM), collapse = ",")
        program_count <- length(unique(unlist(strsplit(programs_concat, ","))))
        matching_rows$PROGRAM <- programs_concat
        matching_rows$N_PROGRAM <- program_count
        matching_rows$PROGRAM_SIMILARITY_PERCENTAGE <- round((min_coherence / max_similarity) * 100, 3)
      }
      if (length(unique(matching_rows$SAMPLE_ID)) > 1) {
        sample_concat <- paste(unique(matching_rows$SAMPLE_ID), collapse = ",")
        sample_count <- length(unique(unlist(strsplit(sample_concat, ","))))
        matching_rows$SAMPLE_ID <- sample_concat
        matching_rows$N_SAMPLE <- sample_count
        matching_rows$SAMPLE_SIMILARITY_PERCENTAGE <- round((min_coherence / max_similarity) * 100, 3)
      }
      matching_rows$N_NUCLEOTIDE_COINCIDENCE <- min_coherence - matching_rows$POS
      matching_rows$ID <- paste(matching_rows$ID, collapse = ",")
      max_info_row <- matching_rows[which.max(nchar(matching_rows$INFO)), ]
      matching_rows$INFO <- max_info_row$INFO
      max_format_row <- matching_rows[which.max(nchar(matching_rows$FORMAT)), ]
      matching_rows$FORMAT <- max_format_row$FORMAT
      max_sample_row <- matching_rows[which.max(nchar(matching_rows$SAMPLE)), ]
      matching_rows$SAMPLE <- max_sample_row$SAMPLE
      matching_rows$QUAL <- matching_rows$QUAL[which.max(nchar(matching_rows$QUAL))]
      matching_rows$REF <- matching_rows$REF[which.max(nchar(matching_rows$REF))]
      matching_rows$ALT <- matching_rows$ALT[which.max(nchar(matching_rows$ALT))]
      row_max_end <- matching_rows[which.max(matching_rows$END), ]
      DUP_pos <- rbind(DUP_pos, row_max_end)
    } else {
      DUP_pos <- rbind(DUP_pos, matching_rows)
    }
  }
}


#In case of similar END values, it retains minimum POS values
DUP_pos_end <- data.frame()

for (cromosoma in unique(DUP_pos$CHROM)) {
  subset_data <- DUP_pos[DUP_pos$CHROM == cromosoma, ]
  for (end_value in unique(subset_data$END)) {
    matching_rows <- subset_data[subset_data$END == end_value, ]
    if (nrow(matching_rows) > 1 & length(unique(matching_rows$GT)) == 1) {
      max_similarity <- max(matching_rows$POS)
      min_coherence <- min(matching_rows$POS)
      if (length(unique(matching_rows$PROGRAM)) > 1) {
        programs_concat <- paste(unique(matching_rows$PROGRAM), collapse = ",")
        program_count <- length(unique(unlist(strsplit(programs_concat, ","))))
        matching_rows$PROGRAM <- programs_concat
        matching_rows$N_PROGRAM <- program_count
        matching_rows$PROGRAM_SIMILARITY_PERCENTAGE <- round((min_coherence / max_similarity) * 100, 3)
      }
      if (length(unique(matching_rows$SAMPLE_ID)) > 1) {
        sample_concat <- paste(unique(matching_rows$SAMPLE_ID), collapse = ",")
        sample_count <- length(unique(unlist(strsplit(sample_concat, ","))))
        matching_rows$SAMPLE_ID <- sample_concat
        matching_rows$N_SAMPLE <- sample_count
        matching_rows$SAMPLE_SIMILARITY_PERCENTAGE <- round((min_coherence / max_similarity) * 100, 3)
      }
      matching_rows$N_NUCLEOTIDE_COINCIDENCE <- matching_rows$END - max_similarity
      matching_rows$ID <- paste(matching_rows$ID, collapse = ",")
      max_info_row <- matching_rows[which.max(nchar(matching_rows$INFO)), ]
      matching_rows$INFO <- max_info_row$INFO
      max_format_row <- matching_rows[which.max(nchar(matching_rows$FORMAT)), ]
      matching_rows$FORMAT <- max_format_row$FORMAT
      max_sample_row <- matching_rows[which.max(nchar(matching_rows$SAMPLE)), ]
      matching_rows$SAMPLE <- max_sample_row$SAMPLE
      matching_rows$QUAL <- matching_rows$QUAL[which.max(nchar(matching_rows$QUAL))]
      matching_rows$REF <- matching_rows$REF[which.max(nchar(matching_rows$REF))]
      matching_rows$ALT <- matching_rows$ALT[which.max(nchar(matching_rows$ALT))]
      row_min_pos <- matching_rows[which.min(matching_rows$POS), ]
      DUP_pos_end <- rbind(DUP_pos_end, row_min_pos)
    } else {
      DUP_pos_end <- rbind(DUP_pos_end, matching_rows)
    }
  }
}
SV_DUP <- DUP_pos_end


#Filtered coincident regions (new POS value = minimum, new END value = maximum)
DUP_coincidence <- data.frame()

for (cromosoma in unique(SV_DUP$CHROM)) {
  subset_data <- SV_DUP[SV_DUP$CHROM == cromosoma, ]
  overlaps <- data.frame()
  for (i in 1:nrow(subset_data)) {
    current_row <- subset_data[i, ]
    overlapping_rows <- subset_data[subset_data$POS <= current_row$END & subset_data$END >= current_row$POS, ]
    if (nrow(overlapping_rows) > 1 & length(unique(overlapping_rows$GT)) == 1) {
      min_pos <- min(overlapping_rows$POS)
      max_end <- max(overlapping_rows$END)
      total_length <- max_end - min_pos + 1 #en el caso que la coincidencia sea igual, para que no sea 0 se usa +1
      overlap_length <- min(current_row$END, max_end) - max(current_row$POS, min_pos) + 1
      overlap_percentage <- round((overlap_length / total_length) * 100, 3)
      current_row$POS <- min_pos
      current_row$END <- max_end
      current_row$N_NUCLEOTIDE_COINCIDENCE <- min(overlapping_rows$END) - max(overlapping_rows$POS)
      current_row$ID <- paste(overlapping_rows$ID, collapse = ",")
      max_info_row <- overlapping_rows[which.max(nchar(overlapping_rows$INFO)), ]
      current_row$INFO <- max_info_row$INFO
      max_format_row <- overlapping_rows[which.max(nchar(overlapping_rows$FORMAT)), ]
      current_row$FORMAT <- max_format_row$FORMAT
      max_sample_row <- overlapping_rows[which.max(nchar(overlapping_rows$SAMPLE)), ]
      current_row$SAMPLE <- max_sample_row$SAMPLE
      current_row$QUAL <- overlapping_rows$QUAL[which.max(nchar(overlapping_rows$QUAL))]
      current_row$REF <- overlapping_rows$REF[which.max(nchar(overlapping_rows$REF))]
      current_row$ALT <- overlapping_rows$ALT[which.max(nchar(overlapping_rows$ALT))]
      svlen_index <- grep("SVLEN=", current_row$INFO)
      if (length(svlen_index) > 0) {
        svlen_new <- max_end - min_pos
        current_row$INFO[svlen_index] <- gsub("SVLEN=[-0-9]+", paste("SVLEN=", svlen_new, sep = ""), current_row$INFO[svlen_index])
      }
      end_index <- grep("END=", current_row$INFO)
      if (length(end_index) > 0) {
        current_row$INFO[end_index] <- gsub("END=[0-9]+", paste("END=", max_end, sep = ""), current_row$INFO[end_index])
      }
      if (length(unique(overlapping_rows$PROGRAM)) > 1) {
        programs_concat <- paste(unique(overlapping_rows$PROGRAM), collapse = ",")
        program_count <- length(unique(unlist(strsplit(programs_concat, ","))))
        current_row$PROGRAM <- programs_concat
        current_row$N_PROGRAM <- program_count
        current_row$PROGRAM_SIMILARITY_PERCENTAGE <- overlap_percentage
      }
      if (length(unique(overlapping_rows$SAMPLE_ID)) > 1) {
        sample_concat <- paste(unique(overlapping_rows$SAMPLE_ID), collapse = ",")
        sample_count <- length(unique(unlist(strsplit(sample_concat, ","))))
        current_row$SAMPLE_ID <- sample_concat
        current_row$N_SAMPLE <- sample_count
        current_row$SAMPLE_SIMILARITY_PERCENTAGE <- overlap_percentage
      }
      non_missing_qual <- overlapping_rows$QUAL[overlapping_rows$QUAL != "."]
      if (length(non_missing_qual) > 0) {
        max_qual <- max(as.numeric(non_missing_qual), na.rm = TRUE)
        current_row$QUAL <- ifelse(is.na(max_qual), ".", as.character(max_qual))
      }
      overlaps <- rbind(overlaps, current_row)
    } else {
      overlaps <- rbind(overlaps, current_row)
    }
  }
  overlaps <- overlaps[order(-as.numeric(sapply(overlaps$INFO, function(x) gsub(".*END=([0-9]+).*", "\\1", x)))), ]
  overlaps <- overlaps[!duplicated(overlaps$POS), ]
  DUP_coincidence <- rbind(DUP_coincidence, overlaps)
}

DUP_coincidence$GT <- NULL
SV_DUP <- DUP_coincidence


########################################################################################################################################################
# INS #
#######


SV_INS <- data.frame()
print("Merging INS")
for (i in 1:nrow(all_progrmas_df)) {
  current_row <- all_progrmas_df[i, ]
  if (grepl("SVTYPE=.*INS", current_row$INFO)) {
    end_value <- as.numeric(sub(".*SVLEN=([0-9]+).*", "\\1", current_row$INFO, perl = TRUE))
    current_row$END <- end_value + current_row$POS
    if (grepl("DR:DV", current_row$FORMAT)) {
      dr_position <- grep("DR", unlist(strsplit(as.character(current_row$FORMAT), ":")))
      dv_position <- grep("DV", unlist(strsplit(as.character(current_row$FORMAT), ":")))
      sample_values <- unlist(strsplit(as.character(current_row$SAMPLE), ":"))
      sample_values[dr_position] <- paste(sample_values[dr_position], sample_values[dv_position], sep = ",")
      sample_values <- sample_values[-dv_position]
      current_row$SAMPLE <- paste(sample_values, collapse = ":")
      current_row$FORMAT <- gsub("DR:DV", "AD", current_row$FORMAT)
    }
    current_row$GT <- unlist(strsplit(as.character(current_row$SAMPLE[1]), ":"))[1]
    SV_INS <- rbind(SV_INS, current_row)
  }
}


#In case of similar POS values, it retains maximum END values
INS_pos <- data.frame()

for (cromosoma in unique(SV_INS$CHROM)) {
  subset_data <- SV_INS[SV_INS$CHROM == cromosoma, ]
  for (pos_value in unique(subset_data$POS)) {
    for (end_value in unique(subset_data$END)) {
      matching_rows <- subset_data[subset_data$POS == pos_value & subset_data$END == end_value, ]
      if (nrow(matching_rows) > 1 & length(unique(matching_rows$GT)) == 1 & length(unique(matching_rows$ALT)) == 1) { 
         if (length(unique(matching_rows$PROGRAM)) > 1) {
            programs_concat <- paste(unique(matching_rows$PROGRAM), collapse = ",")
            program_count <- length(unique(unlist(strsplit(programs_concat, ","))))
            matching_rows$PROGRAM <- programs_concat
            matching_rows$N_PROGRAM <- program_count
            matching_rows$PROGRAM_SIMILARITY_PERCENTAGE <- round((min_coherence / max_similarity) * 100, 3)
         }
         if (length(unique(matching_rows$SAMPLE_ID)) > 1) {
            sample_concat <- paste(unique(matching_rows$SAMPLE_ID), collapse = ",")
            sample_count <- length(unique(unlist(strsplit(sample_concat, ","))))
            matching_rows$SAMPLE_ID <- sample_concat
            matching_rows$N_SAMPLE <- sample_count
            matching_rows$SAMPLE_SIMILARITY_PERCENTAGE <- round((min_coherence / max_similarity) * 100, 3)
         }
         matching_rows$N_NUCLEOTIDE_COINCIDENCE <- matching_rows$END - matching_rows$POS
         matching_rows$ID <- paste(matching_rows$ID, collapse = ",")
         max_info_row <- matching_rows[which.max(nchar(matching_rows$INFO)), ]
         matching_rows$INFO <- max_info_row$INFO
         max_format_row <- matching_rows[which.max(nchar(matching_rows$FORMAT)), ]
         matching_rows$FORMAT <- max_format_row$FORMAT
         max_sample_row <- matching_rows[which.max(nchar(matching_rows$SAMPLE)), ]
         matching_rows$SAMPLE <- max_sample_row$SAMPLE
         matching_rows$QUAL <- matching_rows$QUAL[which.max(nchar(matching_rows$QUAL))]
         matching_rows$REF <- matching_rows$REF[which.max(nchar(matching_rows$REF))]
         matching_rows$ALT <- matching_rows$ALT[which.max(nchar(matching_rows$ALT))]
         row_max_end <- matching_rows[which.max(matching_rows$END), ]
         INS_pos <- unique(rbind(INS_pos, matching_rows))
         } else {
           INS_pos <- unique(rbind(INS_pos, matching_rows))
         }
    }
  }
}
INS_pos$GT <- NULL
SV_INS <- INS_pos


########################################################################################################################################################
# BND #
#######

SV_BND <- data.frame()
print("Merging BND")
for (i in 1:nrow(all_progrmas_df)) {
  current_row <- all_progrmas_df[i, ]
  if (grepl("SVTYPE=.*BND", current_row$INFO)) {
    current_row$END <- sub("\\[([^\\[\\]]*)\\]", "\\1", current_row$ALT, perl = TRUE)
    if (grepl("DR:DV", current_row$FORMAT)) {
      dr_position <- grep("DR", unlist(strsplit(as.character(current_row$FORMAT), ":")))
      dv_position <- grep("DV", unlist(strsplit(as.character(current_row$FORMAT), ":")))
      sample_values <- unlist(strsplit(as.character(current_row$SAMPLE), ":"))
      sample_values[dr_position] <- paste(sample_values[dr_position], sample_values[dv_position], sep = ",")
      sample_values <- sample_values[-dv_position]
      current_row$SAMPLE <- paste(sample_values, collapse = ":")
      current_row$FORMAT <- gsub("DR:DV", "AD", current_row$FORMAT)
    }
    current_row$GT <- unlist(strsplit(as.character(current_row$SAMPLE[1]), ":"))[1]
    SV_BND <- rbind(SV_BND, current_row)
  }
}


#In case of similar POS values, it retains maximum END values
BND_pos <- data.frame()

for (cromosoma in unique(SV_BND$CHROM)) {
  subset_data <- SV_BND[SV_BND$CHROM == cromosoma, ]
  for (pos_value in unique(subset_data$POS)) {
    for (end_value in unique(subset_data$END)) {
     matching_rows_end <- subset_data[subset_data$POS == pos_value & subset_data$END == end_value, ]
      if (nrow(matching_rows_end) > 1 & length(unique(matching_rows_end$GT)) == 1) {
        if (length(unique(matching_rows_end$PROGRAM)) > 1) {
          programs_concat <- paste(unique(matching_rows_end$PROGRAM), collapse = ",")
          program_count <- length(unique(unlist(strsplit(programs_concat, ","))))
          matching_rows_end$PROGRAM <- programs_concat
          matching_rows_end$N_PROGRAM <- program_count
          matching_rows_end$PROGRAM_SIMILARITY_PERCENTAGE <- 100
        }
        if (length(unique(matching_rows_end$SAMPLE_ID)) > 1) {
          sample_concat <- paste(unique(matching_rows_end$SAMPLE_ID), collapse = ",")
          sample_count <- length(unique(unlist(strsplit(sample_concat, ","))))
          matching_rows_end$SAMPLE_ID <- sample_concat
          matching_rows_end$N_SAMPLE <- sample_count
          matching_rows_end$SAMPLE_SIMILARITY_PERCENTAGE <- 100
        }
        matching_rows_end$ID <- paste(matching_rows_end$ID, collapse = ",")
        max_info_row <- matching_rows_end[which.max(nchar(matching_rows_end$INFO)), ]
        matching_rows_end$INFO <- max_info_row$INFO
        max_format_row <- matching_rows_end[which.max(nchar(matching_rows_end$FORMAT)), ]
        matching_rows_end$FORMAT <- max_format_row$FORMAT
        max_sample_row <- matching_rows_end[which.max(nchar(matching_rows_end$SAMPLE)), ]
        matching_rows_end$SAMPLE <- max_sample_row$SAMPLE
        matching_rows_end$QUAL <- matching_rows_end$QUAL[which.max(nchar(matching_rows_end$QUAL))]
        matching_rows_end$REF <- matching_rows_end$REF[which.max(nchar(matching_rows_end$REF))]
        matching_rows_end$ALT <- matching_rows_end$ALT[which.max(nchar(matching_rows_end$ALT))]
        BND_pos <- unique(rbind(BND_pos, matching_rows_end))
      } else {
        BND_pos <- unique(rbind(BND_pos, matching_rows_end))
      }
    }
  }
}
BND_pos$GT <- NULL
SV_BND <- BND_pos


#aportar PRECISE o IMPRECISE a la columna INFO
reorder_info <- function(info) {
  info_elements <- strsplit(info, ";")[[1]]
  imprecise <- grep("^IMPRECISE", info_elements, value = TRUE)
  precise <- grep("^PRECISE", info_elements, value = TRUE)
  extra <- setdiff(info_elements, c(imprecise, precise))
  if (length(precise) > 0) {
    new_info <- paste(c(precise, extra), collapse = ";")
  } else if (length(imprecise) > 0) {
    new_info <- paste(c(imprecise, extra), collapse = ";")
  } else {
    new_info <- paste(c("PRECISE", extra), collapse = ";")
  }
  return(new_info)
}
SV_BND <- SV_BND %>% mutate(INFO = purrr::map_chr(INFO, reorder_info))


########################################################################################################################################################
# OUTPUT #
##########


#unir los data.frame sin BND
SV_merge <- rbind(SV_DEL, SV_DUP, SV_INS, SV_INV)


#ordenar la tabla INFO
reorder_info <- function(info) {
  info_elements <- strsplit(info, ";")[[1]]
  svtype <- grep("^SVTYPE=", info_elements, value = TRUE)
  svlen <- grep("^SVLEN=", info_elements, value = TRUE)
  end <- grep("^END=", info_elements, value = TRUE)
  imprecise <- grep("^IMPRECISE", info_elements, value = TRUE)
  precise <- grep("^PRECISE", info_elements, value = TRUE)
  extra <- setdiff(info_elements, c(svtype, svlen, end, imprecise, precise))
  if (length(precise) > 0) {
    new_info <- paste(c(precise, svtype, svlen, end, extra), collapse = ";")
  } else if (length(imprecise) > 0) {
    new_info <- paste(c(imprecise, svtype, svlen, end, extra), collapse = ";")
  } else {
    new_info <- paste(c("PRECISE", svtype, svlen, end, extra), collapse = ";")
  }
  return(new_info)
}
SV_merge <- SV_merge %>% mutate(INFO = purrr::map_chr(INFO, reorder_info))


#unir todos BND a todos los SV
SV_merge <- rbind(SV_merge, SV_BND)

##ADECUAR LAS COLUMNAS##
#Eliminar repeticiones PROGRAM
unique_programs <- sapply(strsplit(SV_merge$PROGRAM, ","), function(x) {
  if(length(x) > 1) {
    return(paste(unique(x), collapse = ","))
  } else {
    return(x)
  }
})
SV_merge$PROGRAM <- unique_programs

#Eliminar repeticiones SAMPLE_ID
unique_sample_ids <- sapply(strsplit(SV_merge$SAMPLE_ID, ","), function(x) {
  if(length(x) > 1) {
    return(paste(unique(x), collapse = ","))
  } else {
    return(x)
  }
})
SV_merge$SAMPLE_ID <- unique_sample_ids

#nº samples en cada SV
SV_merge$N_SAMPLE <- apply(SV_merge, 1, function(row) { length(unique(unlist(strsplit(as.character(row["SAMPLE_ID"]), ","))))})

#Ajustar la columna SIMILARITY_PERCENTAGE
SV_merge <- SV_merge %>%
  mutate(
    SAMPLE_SIMILARITY_PERCENTAGE = ifelse(SAMPLE_SIMILARITY_PERCENTAGE == 0 & N_SAMPLE == 1, NA, SAMPLE_SIMILARITY_PERCENTAGE),
    PROGRAM_SIMILARITY_PERCENTAGE = ifelse(PROGRAM_SIMILARITY_PERCENTAGE == 0 & N_PROGRAM == 1, NA, PROGRAM_SIMILARITY_PERCENTAGE),
    N_NUCLEOTIDE_COINCIDENCE = ifelse(N_NUCLEOTIDE_COINCIDENCE == 0, NA, N_NUCLEOTIDE_COINCIDENCE)
  )


SV_merge <- SV_merge %>% rename("#CHROM" = 1)
SV_merge$END <- NULL

#generar el output en la pathway indicada
outputfile <- file.path(opt$outputdir, ("result_comparison.vcf"))

print(paste("Generating output:", outputfile))
write.table(SV_merge, outputfile, col.names = T, row.names = F, sep = "\t", quote = F)


