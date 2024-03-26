### SV_comparison_program analysis 
### Author: Rodrigo López
### Date: 14/02/2024

rm(list=ls())


#**********************#
# package requirements #
#**********************#

required_packages <- ("optparse")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
  install.packages(missing_packages)
}

library(optparse)

#==============#
#   Arguments  #
#==============#

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="\t\t input TSV file (output from VEP)", metavar="character"),
  make_option(c('-o','--outputdir'),type="character", help="Output directory.")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# opt=list()
# opt$input = "/home/rodrigo/tblab/rodrigo/data/TSV/LongSeq/result_comparison.SV.annotated.tsv"
# opt$outputdir <- "/home/rodrigo/tblab/rodrigo/data/TSV/LongSeq"


####################
#   Input loading  #
####################


tsv_data = read.delim(opt$input, stringsAsFactors = FALSE)
print(paste("Loading file:",opt$input))


###################
# Input filtering #
###################


max_program_n <- max(tsv_data$PROGRAM_N)
tsv_data_filtered <- data.frame()
print("Filtering: N_PROGRAM == max; N_SAMPLE < 2; LOEUF_bin < 6 or !=NA; AnnotSV_ranking_score >-0.89; ACMG_class >2 or !=NA")
# LOEUF -> loss-of-function observed/expected upper bound fraction: Low LOEUF scores indicate loss-of-function variation in a given gene; High LOEUF scores suggest a relatively higher tolerance to inactivation
# AnnotSV_ranking_score -> SV ranking score:  pathogenic ≥0.99, likely pathogenic [0.90;0.98], variant of uncertain significance [0.89;-0.89], likely benign [-0.90;-0.98], benign ≤-0.99.
# ACMG_class -> SV ranking class into 1 of 5: class 1 (benign); class 2 (likely benign); class 3 (variant of unknown significance); class 4 (likely pathogenic); class 5 (pathogenic); class NA (Non Attributed)

for (row in 1:nrow(tsv_data)) {
  current_row <- tsv_data[row, ]
  if (current_row$N_PROGRAM == max_program_n) {
    if (current_row$N_SAMPLE < 2) {
      if (!is.na(current_row$LOEUF_bin) && current_row$LOEUF_bin < 6) {
        if (is.na(current_row$AnnotSV_ranking_score) || current_row$AnnotSV_ranking_score > -0.89){
          full_valor <- sub("^full=(\\d+)$", "\\1", current_row$ACMG_class)
          if (!is.na(current_row$ACMG_class) && current_row$ACMG_class != "full=NA" && full_valor > 2) {
            tsv_data_filtered <- rbind(tsv_data_filtered, current_row)
          }
        }
      }
    }
  }
}


##########
# Output #
##########

outputfile <- file.path(opt$outputdir, ("result_comparison.SV.annotated.filtered.tsv"))
print(paste("Generating output:", outputfile))
write.table(tsv_data_filtered, outputfile, col.names = T, row.names = F, sep = "\t", quote = F)

