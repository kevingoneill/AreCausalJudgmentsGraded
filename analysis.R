#!/usr/bin/Rscript
library(ggplot2)
library(brms)

args <- commandArgs(trailingOnly=T)
if (length(args) != 1) {
    writeLines('Usage: ./analisis.r <file-name>')
    quit()
}

judgments <- read.csv(args[1], header=TRUE)
judgments$MC_conf <- judgments$MC_conf / 100
judgments$VAS_resp <- judgments$VAS_resp / 100
judgments$VAS_conf <- judgments$VAS_conf / 100
judgments$MC_resp <- factor(judgments$MC_resp, #ordered=TRUE,
                            levels=c('did not cause', 'partially caused', 'totally caused'))

time <- aggregate(duration ~ id, judgments, function(x) mean(x) / 60.0)
mean(time$duration)
sd(time$duration) / sqrt(nrow(time))



## binned_pp_check(model, judgments)
##
##  Compute/display pp_checks for a model for different levels of confidence.
##  Uses type='hist' to avoid problems with bandwidth
##
binned_pp_check <- function(model, judgments, binwidth=1.0) {
    for (l in levels(judgments$MC_resp)) {
        for (upper in seq(min(judgments$MC_conf)+binwidth, max(judgments$MC_conf), by=binwidth)) {
            upper = upper
            lower = upper - binwidth
            
            mask <- (judgments$MC_resp == l) &
                (judgments$MC_conf >= lower) &
                (judgments$MC_conf <= upper)
            j <- judgments[mask, ]
            
            if (!is.null(nrow(j)) && nrow(j) > 1) {
                print(pp_check(model, type='hist', binwidth=0.01, newdata=j) + theme_bw() +
                      ggtitle(sprintf("%s Ratings (%1.2f - %1.2f Confidence)",
                                      as.character(l), lower, upper)))
            }
        }
    }
}



mNormal <- brm(VAS_resp ~ MC_resp * MC_conf + (1 + MC_resp*MC_conf || id) +
                   (1 || vignette:structure:condition),
               data=judgments, file='mNormal', cores=1, init_r = 0.99)
summary(mNormal)
pdf("Normal.pdf")
plot(mNormal)

pp_check(mNormal, resp='VASresp', bw=0.05) + theme_bw()
marginal_effects(mNormal, resp='VASresp', probs=c(0.05, 0.95), points=TRUE)
dev.off()


mZOIB <- brm(bf(VAS_resp ~ MC_resp * MC_conf + (1 + MC_resp*MC_conf | id) +
                    (1 || vignette:structure:condition),
                phi ~ MC_resp * MC_conf + (1 + MC_resp*MC_conf || id) +
                    (1 || vignette:structure:condition),
                zoi ~ MC_resp * MC_conf + (1 + MC_resp*MC_conf || id) +
                    (1 || vignette:structure:condition),
                coi ~ MC_resp * MC_conf + (1 + MC_resp*MC_conf || id) +
                    (1 || vignette:structure:condition)),
             family=zero_one_inflated_beta(), data=judgments, file='mZOIB',
             cores=4, thin=2, iter=2500, init_r=0.99, control=list(adapt_delta=0.95))
summary(mZOIB)
pdf("ZOIB.pdf")
plot(mZOIB)
pp_check(mZOIB, bw=0.005) + theme_bw()
marginal_effects(mZOIB, probs=c(0.05, 0.95), points=TRUE)


loo(mNormal, mZOIB)
