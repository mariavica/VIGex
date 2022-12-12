#!/usr/bin/env Rscript

### load your normalized data (in log2-values)
library(openxlsx)
library(argparser)

argp = arg_parser(description = "Calculate VIGex Score");


argp = add_argument(
  parser  = argp                 ,
  arg     = "--expression"           ,
  help    = "Expression data table"
);

argp = add_argument(
  parser  = argp                 ,
  arg     = "--output"           ,
  help    = "Output file"
);

argv = parse_args(parser = argp);

data <- read.xlsx(argv$expression, sheet=1)

expression <- t(data[,c(-1)])
colnames(expression) <- data$samples

### make sure your data is centered around the mean
centered <- t(scale(t(expression),center = TRUE, scale = FALSE))


vigex_genes <- c("CTLA4","IL7R","GZMA","PRF1","IFNg","PDCD1","FOXP3","CXCL9","CXCL10","CXCL11","GZMB","CD274")


vigex_output <- data.frame(samples=colnames(centered), vigex_score=NA, vigex_cluster=NA)

vigex_output$vigex_score <- apply(centered[vigex_genes,],2,mean)

######################################################################################################################
### Caution: these cutoffs have been validated in the training dataset presented in this paper.                   ####
### If you are using vigex in other samples/techniques, you'll probably have to correct for batch effects.        ####                                                                                        ####
### Moreover, this is a score/sample grouping developed in a pan-cancer cohort, if your cohort is tumor-specific, ####
### you may have to apply further corrections or other techniques to obtain the Hot/iCold/Cold categories.        ####                                                                              ####
######################################################################################################################

vigex_output$vigex_cluster[ which(vigex_output$vigex_score>=0.75) ] <- "Hot"
vigex_output$vigex_cluster[ which(vigex_output$vigex_score<0.75 & vigex_output$vigex_score>c(-0.75)) ] <- "iCold"
vigex_output$vigex_cluster[ which(vigex_output$vigex_score<=c(-0.75)) ] <- "Cold"

write.xlsx(vigex_output,argv$output)




