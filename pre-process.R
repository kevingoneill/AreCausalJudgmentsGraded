#!/usr/bin/Rscript
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly=T)
if (length(args) != 2) {
    writeLines('Usage: ./pre-process.r <input-file-name> <output-file-name>')
    quit()
}
in_file <- args[1]
out_file <- args[2]

data_wide <- read.csv(in_file, header=TRUE, stringsAsFactors=FALSE)

## remove unneeded rows/cols
data_wide <- data_wide[-c(1, 2),]
data_wide <- data_wide[, -c(1:5, 7:21, 169:249, 251)]

## save demographic info
#write.csv(subset(data_wide, select=c('Age', 'Sex', 'AttnCheck')),
#          'demographics.csv', row.names=FALSE)

## filter out by attention check
data_wide <- subset(subset(data_wide, AttnCheck=='Yes.'), select = -c(AttnCheck))

df <- gather(data_wide, question, response, T_JC_N_MC_Resp:Cof_OD_I_VAS_Conf_1)
colnames(df) <- c('duration', 'age', 'sex', 'id', 'question', 'response')
df <- df[order(df$id),]
df <- df[-which(df$response == ''),]   # filter out missing responses

## split question identifier into conditions
df <- separate(df, 'question', c('vignette', 'structure', 'condition', 'scale', 'measure'),
               sep='_', extra='drop')
df$measure <- tolower(df$measure)

## convert back to wide format- one row per vignette, columns for MCresp, MCconf, VASresp, VASconf
df <- df %>% unite(temp, scale, measure) %>% spread(temp, response)
df <- df[, c('id', 'age', 'sex', 'duration', 'vignette', 'structure', 'condition',
             'MC_resp', 'MC_conf', 'VAS_resp', 'VAS_conf')]

write.csv(df, out_file, row.names=FALSE)
