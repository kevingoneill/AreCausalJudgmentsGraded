#!/usr/bin/Rscript
library(ggplot2)
library(brms)
library(plyr)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggisoband)

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
length(unique(judgments$id))
time <- aggregate(duration ~ id, judgments, function(x) mean(x) / 60.0)
mean(time$duration)
sd(time$duration) / sqrt(nrow(time))

marginal_density <- function(judgments, var, group='MC_resp') {
    ggplot(judgments) +
        aes_string(x=var, color=group, fill=group) +
        geom_density(bw=0.01, alpha=0.33) +
        theme_classic() +
        theme_void() + theme(legend.position='none')
}


density_2d <- function(center, top, right) {
    ## align marginal plots
    aligned_x <- align_plots(top, center, align='v', axis='lr')[[1]]
    aligned_y <- align_plots(right, center, align='h', axis='tb')[[1]]
    
    ## Arrange plots
    plot_grid(aligned_x, NULL, center, aligned_y,
              ncol = 2, nrow = 2,
              rel_heights = c(0.25, 1),
              rel_widths = c(1, 0.25))
}

pdf('hist2d.pdf')
density_2d(ggplot(judgments) +
           aes(x=VAS_resp, y=MC_conf, group=MC_resp, color=MC_resp, fill=MC_resp) +
           geom_density_bands(aes(alpha=stat(ndensity)), ##h=c(0.25, 0.25),
                              show.legend=FALSE) +
           ##geom_point() +
           xlab('Causal Rating') + ylab('Confidence') +
           labs(color='Categorical Rating', fill='Categorical Rating') +
           scale_alpha_continuous(range = c(0, 1)) +
           theme_classic() +
           theme(legend.position='bottom'),
           marginal_density(judgments, 'VAS_resp'),
           marginal_density(judgments, 'MC_conf') + coord_flip())
density_2d(ggplot(judgments) +
           aes(x=VAS_resp, y=MC_conf, group=MC_resp, fill=MC_resp) +
           stat_density2d(aes(alpha=stat(ndensity)),
                          geom='tile', contour=FALSE, show.legend=FALSE) +
           xlab('Causal Rating') + ylab('Confidence') +
           labs(color='Categorical Rating', fill='Categorical Rating') +
           scale_alpha_continuous(range = c(0, 1)) +
           theme_classic() +
           theme(legend.position='bottom'),
           marginal_density(judgments, 'VAS_resp'),
           marginal_density(judgments, 'MC_conf') + coord_flip())
dev.off()

## binned_pp_check(model, judgments)
##
##  Compute/display pp_checks for a model for different levels of confidence.
##  Uses type='hist' to avoid problems with bandwidth
##
binned_pp_check <- function(model, judgments, binwidth=1.0) {
    for (l in levels(judgments$MC_resp)) {
        for (upper in seq(binwidth, 1.0, by=binwidth)) {
            upper = upper
            lower = upper - binwidth
            
            mask <- (judgments$MC_resp == l) &
                (judgments$MC_conf >= lower) &
                (judgments$MC_conf <= upper)
            j <- judgments[mask, ]
            
            if (!is.null(nrow(j)) && nrow(j) > 1) {
                print(pp_check(model, bw=0.005, newdata=j) + theme_bw() +
                      ggtitle(sprintf("%s Ratings (%1.2f - %1.2f Confidence)",
                                      as.character(l), lower, upper)))
            }
        }
    }
}

mNormal <- brm(bf(VAS_resp ~ MC_resp * MC_conf + (1 |i| id) +
                      (1 |v| vignette:structure:condition),
                  sigma ~ MC_resp * MC_conf + (1 |i| id) +
                      (1 |v| vignette:structure:condition)),
               prior=c(set_prior('normal(0, 10.0)', class='b'),
                       set_prior('normal(0, 10.0)', class='b', dpar='sigma')),
               data=judgments, file='mNormal2', cores=4, inits="0")
summary(mNormal)


conditions <- data.frame(MC_resp=rep(levels(judgments$MC_resp), each=100),
                         MC_conf=rep((1:100)/100), 3)

effects <- conditions %>%
    cbind(fitted(mNormal, newdata=conditions, re_formula=NA)) %>%
    rename(VAS_resp=Estimate)

#samples <- mNormal %>% fitted(re_formula=NA, summary=FALSE) %>% t %>%
#    cbind(select(judgments, c(MC_resp, MC_conf))) %>%
#    pivot_longer(-c('MC_resp', 'MC_conf'),
#                 names_to='sample', values_to='VAS_resp')

samples <- mNormal %>% fitted() %>%
    cbind(select(judgments, c(MC_resp, MC_conf))) %>%
    rename(VAS_resp=Estimate)

density_2d(effects %>%
           ggplot + aes(x=MC_conf, y=VAS_resp, group=MC_resp) +
           ##geom_density_bands(aes(alpha=stat(ndensity)), h=c(0.2, 0.2),
           ##                   bins=20, show.legend=FALSE) +
           geom_line(aes(color=MC_resp), size=2) +
           geom_ribbon(aes(ymin=Q2.5, ymax=Q97.5, fill=MC_resp), alpha=0.3) +
           coord_flip() +
           xlab('Causal Rating') + ylab('Confidence') +
           labs(color='Categorical Rating', fill='Categorical Rating') +
           ##scale_alpha_continuous(range = c(0, 1)) +
           theme_classic() +
           theme(legend.position='bottom'),
           marginal_density(samples, 'VAS_resp'),
           marginal_density(samples, 'MC_conf') + coord_flip())

density_2d(ggplot(samples) +
           aes(x=VAS_resp, y=MC_conf, group=MC_resp, fill=MC_resp) +
           stat_density2d(aes(alpha=stat(ndensity)),
                          geom='tile', contour=FALSE, show.legend=FALSE) +
           xlab('Causal Rating') + ylab('Confidence') +
           labs(color='Categorical Rating', fill='Categorical Rating') +
           scale_alpha_continuous(range = c(0, 1)) +
           theme_classic() +
           theme(legend.position='bottom'),
           marginal_density(samples, 'VAS_resp'),
           marginal_density(samples, 'MC_conf') + coord_flip())

##pdf("Normal2.pdf")
##plot(mNormal)
##pp_check(mNormal, bw=0.005) + theme_bw()
##binned_pp_check(mNormal, judgments, binwidth=0.5)
##marginal_effects(mNormal)
##dev.off()


mZOIB <- brm(bf(VAS_resp ~ MC_resp * MC_conf + (1 |i| id) +
                    (1 |v| vignette:structure:condition),
                phi ~ MC_resp * MC_conf + (1 |i| id) +
                    (1 |v| vignette:structure:condition),
                zoi ~ MC_resp * MC_conf + (1 |i| id) +
                    (1 |v| vignette:structure:condition),
                coi ~ MC_resp * MC_conf + (1 |i| id) +
                    (1 |v| vignette:structure:condition)),
             family=zero_one_inflated_beta(), data=judgments, file='mZOIB2',
             prior=c(set_prior('normal(0, 10.0)', class='b'),
                     set_prior('normal(0, 10.0)', class='b', dpar='phi'),
                     set_prior('normal(0, 10.0)', class='b', dpar='zoi'),
                     set_prior('normal(0, 10.0)', class='b', dpar='coi')),
             cores=4, iter=2500, inits="0", control=list(adapt_delta=0.95))
summary(mZOIB)
##pdf("ZOIB2.pdf")
##plot(mZOIB)
##pp_check(mZOIB, bw=0.005) + theme_bw()
##binned_pp_check(mZOIB, judgments, binwidth=0.5)
##marginal_effects(mZOIB, probs=c(0.05, 0.95), points=TRUE)
##dev.off()

##loo(mNormal, mZOIB)

effects <- conditions %>%
    cbind(fitted(mZOIB, re_formula=NA, newdata=conditions)) %>%
    rename(VAS_resp=Estimate)

#samples <- mZOIB %>% fitted(re_formula=NA, summary=FALSE) %>% t %>%
#    cbind(select(judgments, c(MC_resp, MC_conf))) %>%
#    pivot_longer(-c('MC_resp', 'MC_conf'),
#                 names_to='sample', values_to='VAS_resp')

samples <- mZOIB %>% fitted() %>%
    cbind(select(judgments, c(MC_resp, MC_conf))) %>%
    rename(VAS_resp=Estimate)

density_2d(effects %>%
           ggplot + aes(x=MC_conf, y=VAS_resp, group=MC_resp) +
           ##geom_density_bands(aes(alpha=stat(ndensity)), h=c(0.2, 0.2),
           ##                   bins=20, show.legend=FALSE) +
           geom_line(aes(color=MC_resp), size=2) +
           geom_ribbon(aes(ymin=Q2.5, ymax=Q97.5, fill=MC_resp), alpha=0.3) +
           coord_flip() +
           xlab('Causal Rating') + ylab('Confidence') +
           labs(color='Categorical Rating', fill='Categorical Rating') +
           ##scale_alpha_continuous(range = c(0, 1)) +
           theme_classic() +
           theme(legend.position='bottom'),
           marginal_density(samples, 'VAS_resp'),
           marginal_density(samples, 'MC_conf') + coord_flip())

density_2d(ggplot(samples) +
           aes(x=VAS_resp, y=MC_conf, group=MC_resp, fill=MC_resp) +
           stat_density2d(aes(alpha=stat(ndensity)),
                          geom='tile', contour=FALSE, show.legend=FALSE) +
           xlab('Causal Rating') + ylab('Confidence') +
           labs(color='Categorical Rating', fill='Categorical Rating') +
           scale_alpha_continuous(range = c(0, 1)) +
           theme_classic() +
           theme(legend.position='bottom'),
           marginal_density(samples, 'VAS_resp'),
           marginal_density(samples, 'MC_conf') + coord_flip())

quit()

kfold(mNormal, mZOIB, save_fits=TRUE)
model_weights(mNormal, mZOIB, weights='loo')
