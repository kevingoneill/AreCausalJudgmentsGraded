#!/usr/bin/Rscript
library(ggplot2)
library(brms)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(tidybayes)
library(emmeans)
library(bayestestR)
library(modelr)
library(wesanderson)
library(patchwork)

## Read and normalize data
judgments <- read.csv('data/processed_data.csv', header=TRUE) %>%
    mutate(MC_conf=MC_conf/100,
           VAS_resp=VAS_resp/100,
           VAS_conf=VAS_conf/100,
           MC_resp=factor(MC_resp, levels=c('partially caused',
                                            'did not cause',
                                            'totally caused')),
           normality=factor(ifelse(condition %in% c('N', 'I'), 'normal', 'abnormal'),
                            levels=c('normal', 'abnormal')))

## load custom plotting options/functions
source('plot.R')

### Print out descriptives
writeLines(sprintf('# of participants: %d', length(unique(judgments$id))))
subjData <- judgments %>% group_by(id) %>%
    summarize(duration=duration[1]/60,
              sex=sex[1], age=age[1])
writeLines(sprintf('Duration: %2.2f (%2.2f)',
                   mean(subjData$duration), sd(subjData$duration)))
writeLines(sprintf('Age: %2.2f (%2.2f)', mean(subjData$age), sd(subjData$age)))
subjData %>% select(sex) %>% table

writeLines(sprintf('Number of low-confidence discrete judgments: %d (%f percent)',
                   sum(judgments$MC_conf < 0.5),
                   sum(judgments$MC_conf < 0.5) / nrow(judgments) * 100))
writeLines(sprintf('Number of low-confidence continuous judgments: %d (%f percent)',
                   sum(judgments$VAS_conf < 0.5),
                   sum(judgments$VAS_conf < 0.5) / nrow(judgments) * 100))


## Fit heteroscedastic Gaussian model
mNormal <- brm(bf(VAS_resp ~ MC_resp * MC_conf + (1 |i| id) +
                      (1 |v| vignette:structure:condition),
                  sigma ~ MC_resp * MC_conf + (1 |i| id) +
                      (1 |v| vignette:structure:condition)),
               prior=c(set_prior('normal(0, 1.0)'),
                       set_prior('normal(0, 5.0)', dpar='sigma')),
               save_pars=save_pars(all=TRUE), data=judgments,
               file='mNormal', cores=4, iter=5000, inits="0")
summary(mNormal)


## Gather posterior draws on the linear scale for contrasts
em <- emmeans(mNormal, ~ MC_resp*MC_conf, at=list(MC_conf=0:1))
em.sd <- emmeans(mNormal, ~ MC_resp*MC_conf, at=list(MC_conf=0:1), dpar='sigma')

#########################################################################################
## Descriptives:
##

## plot marginals
figConfMarginal <- marginal_density(judgments, 'MC_conf')
figRatingMarginal <- marginal_density(judgments, 'VAS_resp') + coord_flip()

## plot raw data
fig2A <- dplot(judgments, legend=FALSE)

## plot model predictions
predictions <- add_predicted_draws(judgments, mNormal) %>% ungroup

fig2B <- predictions %>%
    filter(.draw %in% sample(min(.draw):max(.draw), 50)) %>%
    dplot(y='.prediction', ylab='Predicted Causal Judgment', bw.x=0.15, legend=FALSE)

figPredictionsMarginal <- marginal_density(predictions, '.prediction') +
    coord_flip()

fig2C <- mNormal %>%
    emmeans(~ MC_resp*MC_conf, at=list(MC_conf=seq(0, 1, 0.01))) %>%
    as.data.frame %>%
    mplot(y='emmean', ymin='lower.HPD', ymax='upper.HPD',
          ylab='Estimated Mean Causal Judgment')

fig2D <- mNormal %>%
    emmeans(~ MC_resp*MC_conf, at=list(MC_conf=seq(0, 1, 0.01)),
            dpar='sigma', type='response') %>%
    as.data.frame %>%
    mplot(y='response', ymin='lower.HPD', ymax='upper.HPD',
          ylab='Estimated Standard Deviation\nof Causal Judgments')


figConfMarginal + plot_spacer() + figConfMarginal + plot_spacer() +
    fig2A + figRatingMarginal + fig2B + figPredictionsMarginal +
        fig2C + plot_spacer() + fig2D + plot_spacer() +
        plot_layout(nrow=3, heights=c(0.3, 1, 1),
                    ncol=4, widths=c(1, 0.3, 1, 0.3),
                    guides='collect') +
        plot_annotation(tag_levels=list(c('', '', 'A', '', 'B', '', 'C', 'D'))) &
        theme(plot.tag=element_text(size=28),
              plot.tag.position=c(0.05, 1.1),
              legend.position='bottom')

ggsave('Figure2_j.png', width=12.5, height=11)



(marginal_density(judgments, 'MC_conf')) + plot_spacer() +
    (dplot(judgments, legend=FALSE) +
     scale_x_continuous(name='Confidence Rating', labels=c('0', '.25', '.5', '.75', '1')) +
     scale_y_continuous(labels=c('0', '.25', '.5', '.75', '1'))) +
    (marginal_density(judgments, 'VAS_resp') + coord_flip()) +
        plot_layout(nrow=2, ncol=2, heights=c(.3, 1), widths=c(1, .3))
ggsave('data_marginals.png', height=5, width=5)

plot_spacer() + plot_spacer() +
    (fig2D +
     theme(axis.text=element_text(size=12),
           axis.title=element_text(size=24),
           legend.position='none') +
     scale_x_continuous(name='Confidence Rating', labels=c('0', '.25', '.5', '.75', '1')) +
     scale_y_continuous(name='SD of Causal Ratings', labels=c('0', '.25', '.5', '.75', '1'))) +
    plot_spacer() +
    plot_layout(nrow=2, ncol=2, heights=c(.3, 1), widths=c(1, .3))
ggsave('sd_marginals.png', height=5, width=5)


for (i in list(3, 2:3, 1:3)) {
    l <- levels(judgments$MC_resp)[i]
    j <- filter(judgments, MC_resp %in% l)
    pal <- PALETTE[i]
    p <- (marginal_density(j, 'MC_conf', palette=pal)) + plot_spacer() +
        (dplot(j, legend=FALSE, palette=pal) +
         theme(axis.text=element_text(size=12),
               axis.title=element_text(size=24),
               legend.position='none') +
         scale_x_continuous(name='Confidence Rating', labels=c('0', '.25', '.5', '.75', '1')) +
         scale_y_continuous(name='Causal Rating', labels=c('0', '.25', '.5', '.75', '1'))) +
        (marginal_density(j, 'VAS_resp', palette=pal) + coord_flip()) +
        plot_layout(nrow=2, ncol=2, heights=c(.3, 1), widths=c(1, .3))
    ggsave(paste0('data_marginals_', length(l), '.png'), p, height=5, width=5)

    p <- (marginal_density(j, 'MC_conf', palette=pal)) + plot_spacer() +
        (mNormal %>%
         emmeans(~ MC_resp*MC_conf, at=list(MC_conf=seq(0, 1, 0.01))) %>%
         as.data.frame %>%
         filter(MC_resp %in% l) %>%
         mplot(y='emmean', ymin='lower.HPD', ymax='upper.HPD',
               ylab='Estimated Mean Causal Rating', palette=pal) +
         theme(axis.text=element_text(size=12),
               axis.title=element_text(size=24),
               legend.position='none') +
         scale_x_continuous(name='Confidence Rating', labels=c('0', '.25', '.5', '.75', '1')) +
         scale_y_continuous(name='Predicted Causal Rating', labels=c('0', '.25', '.5', '.75', '1'))) +
        (marginal_density(filter(predictions, MC_resp %in% l), '.prediction', palette=pal) + coord_flip()) +
        plot_layout(nrow=2, ncol=2, heights=c(.3, 1), widths=c(1, .3))
    ggsave(paste0('predictions_marginals_', length(l), '.png'), p, height=5, width=5)
}



##    Perform model testing on alternate models
##
mReduced <- brm(bf(VAS_resp ~ MC_resp * MC_conf + (1 |i| id) +
                       (1 |v| vignette:structure:condition)),
                prior=c(set_prior('normal(0, 1.0)')),
                save_all_pars=TRUE, data=judgments,
                file='mReduced', cores=4, iter=5000, inits="0")

mZOIB <- brm(bf(VAS_resp ~ MC_resp * MC_conf + (1 |i| id) +
                    (1 |v| vignette:structure:condition),
                phi ~ MC_resp * MC_conf + (1 |i| id) +
                    (1 |v| vignette:structure:condition),
                zoi ~ MC_resp * MC_conf + (1 |i| id) +
                    (1 |v| vignette:structure:condition),
                coi ~ MC_resp * MC_conf + (1 |i| id) +
                    (1 |v| vignette:structure:condition)),
             family=zero_one_inflated_beta(), data=judgments,
             file='mZOIB', save_all_pars=TRUE, cores=4,
             prior=c(set_prior('normal(0, 10.0)', class='b'),
                     set_prior('normal(0, 10.0)', class='b', dpar='phi'),
                     set_prior('normal(0, 10.0)', class='b', dpar='zoi'),
                     set_prior('normal(0, 10.0)', class='b', dpar='coi')),
             iter=5000, inits="0", control=list(adapt_delta=0.99))

mGAM <- brm(VAS_resp ~ s(MC_conf, by=MC_resp) + MC_resp + (1 |i| id) +
                (1 |v| vignette:structure:condition),
            save_all_pars=TRUE, data=judgments, file='mGAM',
            iter=5000, cores=4, control=list(adapt_delta=0.999,
                                             max_treedepth=25))

LOO(mNormal, mReduced, mZOIB, mGAM)
model_weights(mNormal, mReduced, mZOIB, mGAM)









#########################################################################################
##
## Mean Causal Judgment:
## 
##
## Set ROPE width
ROPE <- sd(judgments$VAS_resp) * 0.1

## Perform contrasts
c1.mean <- em %>% contrast('trt.vs.ctrl', simple='MC_resp') %>%
    describe_posterior(ci=0.95, rope_ci=.95, rope_range=c(-ROPE, ROPE)) %>%
    mutate(Parameter=paste0('Mean, ', Parameter),
           P=pd_to_p(pd))
c1.mean

em %>% contrast('trt.vs.ctrl', simple='MC_resp') %>% gather_emmeans_draws %>%
    ggplot(aes(x=contrast, y=.value, group=contrast,
               fill=stat(ifelse(abs(y) < ROPE, '0', group)))) +
    stat_eye(position=position_dodge(width=1), n=10000, show.legend=TRUE) +
    facet_wrap(~ MC_conf,
               labeller=labeller(MC_conf=c('0'='Low Confidence',
                                           '1'='High Confidence'))) +
    geom_hline(yintercept=c(-ROPE, ROPE), linetype='dashed') +
    scale_fill_manual(values=c('gray80', PALETTE[-1])) +
    scale_x_discrete(labels=c('DNC - PC', 'TC - PC'), name='Contrast') +
    ylab('Mean Causal Rating Contrasts:\nDiscrete Causal Judgmnet') +
    theme_classic() + theme(legend.position='none')
ggsave('contrast-discrete.png', width=6, height=4)

c2.mean <- em %>% contrast('trt.vs.ctrl', simple='MC_conf') %>%
    describe_posterior(ci=0.95, rope_ci=0.95, rope_range=c(-ROPE, ROPE)) %>%
    mutate(Parameter=paste0('Mean, ', Parameter),
           P=pd_to_p(pd))
c2.mean

em %>% contrast('trt.vs.ctrl', simple='MC_conf') %>% gather_emmeans_draws %>%
    ggplot(aes(x=MC_resp, y=.value,
               fill=stat(ifelse(abs(y) < ROPE, '0', group)))) +
    stat_eye(position=position_dodge(width=1), n=10000) +
    geom_hline(yintercept=c(-ROPE, ROPE), linetype='dashed') +
    scale_fill_manual(values=c('gray80', PALETTE)) +
    ylab('Mean Causal Judgment Contrasts: Confidence') +
    xlab('Discrete Causal Judgment') +
    theme_classic() + theme(legend.position='none')
ggsave('contrast-confidence.png', width=6, height=4)

c3.mean <- em %>% contrast('trt.vs.ctrl', interaction=TRUE) %>%
    describe_posterior(ci=0.95, rope_ci=0.95, rope_range=c(-ROPE, ROPE)) %>%
    mutate(Parameter=paste0('Mean, ', Parameter),
           P=pd_to_p(pd))
c3.mean



em %>% contrast('trt.vs.ctrl', interaction=TRUE) %>% gather_emmeans_draws %>%
    rename(MC_resp=MC_resp_trt.vs.ctrl,
           MC_conf=MC_conf_trt.vs.ctrl) %>%
    ggplot(aes(x=MC_resp, y=.value,
               fill=stat(ifelse(abs(y) < ROPE, '0', group)))) +
    stat_eye(position=position_dodge(width=1), n=10000) +
    geom_hline(yintercept=c(-ROPE, ROPE), linetype='dashed') +
    scale_fill_manual(values=c('gray80', PALETTE[-1])) +
    ylab('Mean Causal Rating Contrasts:\nConfidence x Discrete Causal Judgment') +
    xlab('Discrete Causal Judgment') +
    theme_classic() + theme(legend.position='none')
ggsave('contrast-discreteXconfidence.png', width=6, height=4)



#########################################################################################
##
## Standard Deviation of Causal Judgments:
## 
##
## Set ROPE width
sdROPE <- sd(fitted(mNormal, dpar='sigma', scale='linear', summary=FALSE)) * 0.1

## Contrasts for sd(causal rating)
c1.sd <- em.sd %>% contrast('trt.vs.ctrl', simple='MC_resp') %>%
    describe_posterior(ci=0.95, rope_ci=0.95, rope_range=c(-sdROPE, sdROPE)) %>%
    mutate(Parameter=paste0('SD, ', Parameter),
           P=pd_to_p(pd))
c1.sd

em.sd %>% contrast('trt.vs.ctrl', simple='MC_resp') %>%
    gather_emmeans_draws('sigma') %>%
    ggplot(aes(x=contrast, y=sigma, group=contrast,
               fill=stat(ifelse(abs(y) < ROPE, '0', group)))) +
    stat_eye(position=position_dodge(width=1), n=10000) +
    facet_wrap(~ MC_conf, labeller=labeller(MC_conf=c('0'='Low Confidence', '1'='High Confidence'))) +
    geom_hline(yintercept=c(-ROPE, ROPE), linetype='dashed') +
    scale_fill_manual(values=c('gray80', PALETTE[-1])) +
    scale_x_discrete(labels=c('DNC - PC', 'TC - PC'), name='Contrast') +
    ylab('Standard Deviation of Causal Judgment Contrasts:\nDiscrete Causal Judgment') +
    theme_classic() + theme(legend.position='none')
ggsave('contrast-discrete-sd.png', width=6, height=4)



c2.sd <- em.sd %>% contrast('trt.vs.ctrl', simple='MC_conf') %>%
    describe_posterior(ci=0.95, rope_ci=0.95, rope_range=c(-sdROPE, sdROPE)) %>%
    mutate(Parameter=paste0('SD, ', Parameter),
           P=pd_to_p(pd))
c2.sd

em.sd %>% contrast('trt.vs.ctrl', simple='MC_conf') %>%
    gather_emmeans_draws('sigma') %>%
    ggplot(aes(x=MC_resp, y=sigma,
               fill=stat(ifelse(abs(y) < ROPE, '0', group)))) +
    stat_eye(position=position_dodge(width=1), n=10000) +
    geom_hline(yintercept=c(-ROPE, ROPE), linetype='dashed') +
    scale_fill_manual(values=c('gray80', PALETTE)) +
    ylab('Standard Deviation of Causal Judgment Contrasts:\nConfidence') +
    xlab('Discrete Causal Judgment') +
    theme_classic() + theme(legend.position='none')
ggsave('contrast-confidence-sd.png', width=6, height=4)

c3.sd <- em.sd %>% contrast('trt.vs.ctrl', interaction=TRUE) %>%
    describe_posterior(ci=0.95, rope_ci=0.95,
                       rope_range=c(-sdROPE, sdROPE)) %>%
    mutate(Parameter=paste0('SD, ', Parameter),
           P=pd_to_p(pd))
c3.sd

em.sd %>% contrast('trt.vs.ctrl', interaction=TRUE) %>%
    gather_emmeans_draws('sigma') %>%
    rename(MC_resp=MC_resp_trt.vs.ctrl,
           MC_conf=MC_conf_trt.vs.ctrl) %>%
    ggplot(aes(x=MC_resp, y=sigma,
               fill=stat(ifelse(abs(y) < ROPE, '0', group)))) +
    stat_eye(position=position_dodge(width=1), n=10000) +
    geom_hline(yintercept=c(-ROPE, ROPE), linetype='dashed') +
    scale_fill_manual(values=c('gray80', PALETTE[-1])) +
    ylab('Standard Deviation of Causal Judgment Contrasts:\nConfidence x Discrete Causal Judgment') +
    xlab('Discrete Causal Judgment') +
    theme_classic() + theme(legend.position='none')
ggsave('contrast-discreteXconfidence-sd.png', width=6, height=4)


## Export contrast stats
rbind(c2.mean, c2.sd) %>%
    separate(Parameter, into=c('Parameter', 'Confidence', 'DiscreteRating'),
             sep=', ') %>%
    rbind(rbind(c1.mean, c1.sd, c3.mean, c3.sd) %>%
          separate(Parameter, into=c('Parameter', 'DiscreteRating', 'Confidence'),
                   sep=', ')) %>% as.data.frame %>%
    mutate(Median=round(Median, 2), pd=round(pd, 2),
           ROPE_Percentage=round(ROPE_Percentage*100, 2),
           CI=sprintf('[%.2f, %.2f]', CI_low, CI_high),
           ROPE=sprintf('[%.2f, %.2f]', ROPE_low, ROPE_high),
           DiscreteRating=factor(str_replace_all(DiscreteRating,
                                                 c('partially caused'='PC',
                                                   'totally caused'='TC',
                                                   'did not cause'='DNC')),
                                 levels=c('DNC', 'PC', 'TC', 'DNC - PC', 'TC - PC')),
           Confidence=factor(str_replace_all(Confidence, c('0'='Low', '1'='High')),
                             levels=c('Low', 'High', 'High - Low'))) %>%
    select(Parameter, DiscreteRating, Confidence,
           Median, CI, pd, P, ROPE, ROPE_Percentage) %>%
    arrange(Parameter, Confidence, DiscreteRating) %>%
    write.csv('contrasts.csv', row.names=F)




## Plot model coefficient posteriors
mNormal %>% gather_draws(b_Intercept, b_MC_respdidnotcause,
                         b_MC_resptotallycaused,
                         b_MC_conf, `b_MC_respdidnotcause:MC_conf`,
                         `b_MC_resptotallycaused:MC_conf`) %>%
    ungroup() %>%
    mutate(.variable=
               factor(str_replace_all(.variable,
                                      c('b_'='', ':'=' : ',
                                        'MC_respdidnotcause'=
                                            'Discrete Causal Judgment\n[did not cause - partially caused]',
                                        'MC_resptotallycaused'=
                                            'Discrete Causal Judgment\n[totally caused - partially caused]',
                                        'MC_conf'='Confidence')),
                      levels=c('Discrete Causal Judgment\n[totally caused - partially caused] : Confidence',
                               'Discrete Causal Judgment\n[did not cause - partially caused] : Confidence',
                               'Confidence',
                               'Discrete Causal Judgment\n[totally caused - partially caused]',
                               'Discrete Causal Judgment\n[did not cause - partially caused]',
                               'Intercept'))) %>%
    ggplot(aes(y=.variable, x=.value,
               fill=stat(abs(x) > ROPE))) +
    xlab('Estimate') + ylab('') + stat_halfeye() +
    geom_vline(xintercept=c(-ROPE, ROPE), linetype='dashed') +
    scale_fill_manual(values=c('gray80', wes_palette("Darjeeling1", n=5)[5])) +
    theme_classic() + theme(legend.position='none') +
    plot_annotation(title='Coefficient Estimates for Experiment 1:',
                    subtitle='Mean Causal Judgment')
ggsave('coefficients_mu.png', width=8, height=6)

mNormal %>% gather_draws(b_sigma_Intercept, b_sigma_MC_respdidnotcause,
                         b_sigma_MC_resptotallycaused,
                         b_sigma_MC_conf, `b_sigma_MC_respdidnotcause:MC_conf`,
                         `b_sigma_MC_resptotallycaused:MC_conf`) %>%
    ungroup() %>%
    mutate(.variable=
               factor(str_replace_all(.variable,
                                      c('b_'='', ':'=' : ', 'sigma_'='',
                                        'MC_respdidnotcause'=
                                            'Discrete Causal Judgment\n[did not cause - partially caused]',
                                        'MC_resptotallycaused'=
                                            'Discrete Causal Judgment\n[totally caused - partially caused]',
                                        'MC_conf'='Confidence')),
                      levels=c('Discrete Causal Judgment\n[totally caused - partially caused] : Confidence',
                               'Discrete Causal Judgment\n[did not cause - partially caused] : Confidence',
                               'Confidence',
                               'Discrete Causal Judgment\n[totally caused - partially caused]',
                               'Discrete Causal Judgment\n[did not cause - partially caused]',
                               'Intercept'))) %>%
    ggplot(aes(y=.variable, x=.value,
               fill=stat(abs(x) > sdROPE))) +
    xlab('Estimate') + ylab('') + stat_halfeye() +
    geom_vline(xintercept=c(-sdROPE, sdROPE), linetype='dashed') +
    scale_fill_manual(values=c('gray80', wes_palette("Darjeeling1", n=5)[5])) +
    theme_classic() + theme(legend.position='none') +
    plot_annotation(title='Coefficient Estimates for Experiment 1:',
                            subtitle='Standard Deviation of Causal Judgment')
ggsave('coefficients_sigma.png', width=8, height=6)


group_effects <- mNormal %>% gather_draws(sd_id__Intercept, sd_id__sigma_Intercept,
                         `sd_vignette:structure:condition__Intercept`,
                         `sd_vignette:structure:condition__sigma_Intercept`,
                         cor_id__Intercept__sigma_Intercept,
                         `cor_vignette:structure:condition__Intercept__sigma_Intercept`) %>%
    ungroup() %>%
    mutate(group=ifelse(str_detect(.variable, 'id'), 'Participant', 'Vignette'),
           .variable=
               factor(str_replace_all(.variable,
                                      c('sd_id__Intercept'='σ(Intercept)',
                                        'sd_id__sigma_Intercept'='σ(Intercept_σ)',
                                        'sd_vignette:structure:condition__Intercept'=
                                            'σ(Intercept)',
                                        'sd_vignette:structure:condition__sigma_Intercept'=
                                            'σ(Intercept_σ)',
                                        'cor_id__Intercept__sigma_Intercept'=
                                            'correlation(Intercept, Intercept_σ)',
                                        'cor_vignette:structure:condition__Intercept__sigma_Intercept'=
                                            'correlation(Intercept, Intercept_σ)')),
                      levels=c('correlation(Intercept, Intercept_σ)',
                               'σ(Intercept_σ)', 'σ(Intercept)')))

group_effects %>%
    ggplot(aes(y=.variable, x=.value)) +
    xlab('Estimate') + ylab('') +
    stat_halfeye(normalize='xy', fill=wes_palette("Darjeeling1", n=5)[5]) +
    facet_grid(group ~ .) +
    theme_classic() + theme(legend.position='none')


ge1 <- group_effects %>%
    filter(!str_detect(.variable, 'correlation')) %>%
    ggplot(aes(y=.variable, x=.value)) +
    xlab('Estimate') + ylab('') + coord_cartesian(xlim=c(0, 0.6)) +
    stat_halfeye(normalize='xy', fill=wes_palette("Darjeeling1", n=5)[5]) +
    facet_grid(group ~ .) +
    theme_classic() + theme(legend.position='none')

ge2 <- group_effects %>%
    filter(str_detect(.variable, 'correlation')) %>%
    ggplot(aes(y=.variable, x=.value)) +
    xlab('Estimate') + ylab('') +
    stat_halfeye(fill=wes_palette("Darjeeling1", n=5)[5]) +
    facet_grid(group ~ .) +
    theme_classic() + theme(legend.position='none')

ge1 / ge2 + plot_annotation(title='Coefficient Estimates for Experiment 1:',
                            subtitle='Group-Level Effects',
                            tag_levels = 'A')
ggsave('coefficients_re.png', width=8, height=6)

group_effects %>%
    filter(!str_detect(.variable, 'correlation')) %>%
    pivot_wider(names_from=c(.variable, group), values_from=.value) %>%
    select(-.chain, -.iteration, -.draw) %>%
    describe_posterior(., ci=0.95, test=c())

group_effects %>%
    filter(str_detect(.variable, 'correlation')) %>%
    pivot_wider(names_from=c(.variable, group), values_from=.value) %>%
    select(-.chain, -.iteration, -.draw) %>%
    describe_posterior(., ci=0.95, rope_ci=0.95,
                       rope_range=c(-sd(unlist(.))*0.1, sd(unlist(.))*0.1)) %>%
    mutate(P=pd_to_p(pd))




judgments %>%
    mutate(vignette=factor(vignette, levels=c('B', 'T', 'C',
                                              'D', 'S', 'Bu',
                                              'I', 'W', 'Cof'))) %>%
    ggplot(aes(x=MC_conf, y=VAS_resp, color=MC_resp, fill=MC_resp)) +
    geom_smooth(method='lm', size=2, fullrange=TRUE, alpha=0.2) +
    facet_wrap(~ vignette,
                 labeller=labeller(vignette=c('B'='Battery', 'Bu'='Building',
                                              'C'='Computer', 'Cof'='Coffee',
                                              'D'='Dice', 'I'='Implosion',
                                              'S'='Sprinkler', 'T'='Train',
                                              'W'='Watch'))) +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1), expand=FALSE) +
    scale_x_continuous(labels=c('0', '.25', '5', '.75', '1')) +
    scale_y_continuous(labels=c('0', '.25', '5', '.75', '1')) +
    xlab('Confidence') + ylab('Mean Causal Causal Judgment') +
    theme_GC() + theme(panel.spacing=unit(1.5, "lines"),
                       legend.position='bottom')
ggsave('data_vignette.png', width=10, height=10)


judgments %>%
    mutate(condition=factor(ifelse(condition=='N', 'Normal',
                            ifelse(condition=='Ab', 'Abnormal',
                            ifelse(condition=='Ac', 'Action', 'Inaction'))),
                            levels=c('Normal', 'Abnormal', 'Inaction', 'Action')),
           structure=ifelse(structure=='JC', 'Conjunctive', 'Disjunctive')) %>%
    ggplot(aes(x=MC_conf, y=VAS_resp, color=MC_resp, fill=MC_resp)) +
    geom_smooth(method='lm', size=2, fullrange=TRUE, alpha=0.2) +
    facet_grid(structure ~ condition) +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1), expand=FALSE) +
    scale_x_continuous(labels=c('0', '.25', '5', '.75', '1')) +
    scale_y_continuous(labels=c('0', '.25', '5', '.75', '1')) +
    xlab('Confidence') + ylab('Mean Causal Causal Judgment') +
    theme_GC() + theme(panel.spacing=unit(1.5, "lines"),
                       legend.position='bottom')
ggsave('data_conditions.png', width=15, height=10)







#########################################################################################
##
## Use model comparison to test our effects (Table S1.1)
## 
mNormal.add <- brm(bf(VAS_resp ~ MC_resp + MC_conf + (1 |i| id) +
                      (1 |v| vignette:structure:condition),
                  sigma ~ MC_resp + MC_conf + (1 |i| id) +
                      (1 |v| vignette:structure:condition)),
               prior=c(set_prior('normal(0, 1.0)'),
                       set_prior('normal(0, 5.0)', dpar='sigma')),
               save_pars=save_pars(all=TRUE), data=judgments,
               file='mNormal2', cores=4, iter=5000, inits="0")

mNormal.str <- brm(bf(VAS_resp ~ MC_resp + (1 |i| id) +
                      (1 |v| vignette:structure:condition),
                  sigma ~ MC_resp + (1 |i| id) +
                      (1 |v| vignette:structure:condition)),
               prior=c(set_prior('normal(0, 1.0)'),
                       set_prior('normal(0, 5.0)', dpar='sigma')),
               save_pars=save_pars(all=TRUE), data=judgments,
               file='mNormal3', cores=4, iter=5000, inits="0")

mNormal.conf <- brm(bf(VAS_resp ~ MC_conf + (1 |i| id) +
                      (1 |v| vignette:structure:condition),
                  sigma ~ MC_conf + (1 |i| id) +
                      (1 |v| vignette:structure:condition)),
               prior=c(set_prior('normal(0, 1.0)'),
                       set_prior('normal(0, 5.0)', dpar='sigma')),
               save_pars=save_pars(all=TRUE), data=judgments,
               file='mNormal4', cores=4, iter=5000, inits="0")

mNormal.null <- brm(bf(VAS_resp ~ 1 + (1 |i| id) +
                           (1 |v| vignette:structure:condition),
                       sigma ~ 1 + (1 |i| id) +
                           (1 |v| vignette:structure:condition)),
                    save_pars=save_pars(all=TRUE), data=judgments,
                    file='mNormal5', cores=4, iter=5000, inits="0")

loo(mNormal, mNormal.add, mNormal.str, mNormal.conf, mNormal.null, moment_match=TRUE)
model_weights(mNormal, mNormal.add, mNormal.str, mNormal.conf, mNormal.null)




#########################################################################################
##
## Replicate analyses with confidence in *continuous* causal judgment
## 

## Estimate correlation between MC_conf and VAS_conf
mConf <- brm(mvbind(MC_conf, VAS_conf) ~ 1,
             save_pars=save_pars(all=TRUE), data=judgments,
             file='mConfCorr', cores=4, iter=5000)
summary(mConf)

mConf %>% spread_draws(rescor__MCconf__VASconf) %>%
    pull(rescor__MCconf__VASconf) %>%
    describe_posterior(ci=0.95, rope_ci=0.95) %>%
    mutate(P=pd_to_p(pd))

ggplot(judgments, aes(x=MC_conf, y=VAS_conf, color=MC_resp)) +
    geom_point() +
    geom_smooth(method='lm')

## replicate using VAS_conf instead of MC_conf
mNormal.VAS <- brm(bf(VAS_resp ~ MC_resp * VAS_conf + (1 |i| id) +
                      (1 |v| vignette:structure:condition),
                  sigma ~ MC_resp * VAS_conf + (1 |i| id) +
                      (1 |v| vignette:structure:condition)),
               prior=c(set_prior('normal(0, 1.0)'),
                       set_prior('normal(0, 5.0)', dpar='sigma')),
               save_pars=save_pars(all=TRUE), data=judgments,
               file='mNormalVAS', cores=4, iter=5000, inits="0")

mNormal.VAS.add <- brm(bf(VAS_resp ~ MC_resp + VAS_conf + (1 |i| id) +
                      (1 |v| vignette:structure:condition),
                  sigma ~ MC_resp + VAS_conf + (1 |i| id) +
                      (1 |v| vignette:structure:condition)),
               prior=c(set_prior('normal(0, 1.0)'),
                       set_prior('normal(0, 5.0)', dpar='sigma')),
               save_pars=save_pars(all=TRUE), data=judgments,
               file='mNormalVAS2', cores=4, iter=5000, inits="0")

mNormal.VAS.conf <- brm(bf(VAS_resp ~ VAS_conf + (1 |i| id) +
                      (1 |v| vignette:structure:condition),
                  sigma ~ VAS_conf + (1 |i| id) +
                      (1 |v| vignette:structure:condition)),
               prior=c(set_prior('normal(0, 1.0)'),
                       set_prior('normal(0, 5.0)', dpar='sigma')),
               save_pars=save_pars(all=TRUE), data=judgments,
               file='mNormalVAS4', cores=4, iter=5000, inits="0")

summary(mNormal.VAS)


loo(mNormal.VAS, mNormal.VAS.add, mNormal.str, mNormal.VAS.conf, mNormal.null, moment_match=TRUE)
model_weights(mNormal.VAS, mNormal.VAS.add, mNormal.str, mNormal.VAS.conf, mNormal.null)



## Marginal effects of confidence
emtrends(mNormal.VAS, ~MC_resp, var='VAS_conf')


## Plot marginal means & model predictions
figConfMarginal <- marginal_density(judgments, 'VAS_conf')
figRatingMarginal <- marginal_density(judgments, 'VAS_resp') + coord_flip()

## plot raw data
fig2A <- dplot(judgments, legend=FALSE)
fig2A

## plot model predictions
predictions <- add_predicted_draws(judgments, mNormal.VAS) %>% ungroup

fig2B <- predictions %>%
    ##filter(.draw %in% sample(min(.draw):max(.draw), 50)) %>%
    dplot(y='.prediction', ylab='Predicted Causal Judgment', bw.x=0.15, legend=FALSE)

figPredictionsMarginal <- marginal_density(predictions, '.prediction') +
    coord_flip()


fig2C <- mNormal.VAS %>%
    emmeans(~ MC_resp*VAS_conf, at=list(VAS_conf=seq(0, 1, 0.01))) %>%
    as.data.frame %>%
    ggplot(aes(x=VAS_conf, y=emmean,
               ymin=lower.HPD, ymax=upper.HPD,
               color=MC_resp, fill=MC_resp)) +
    geom_line(size=2) + geom_ribbon(alpha=0.3) +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
    xlab('Confidence') + ylab('Estimated Mean Causal Judgment') + theme_GC()


fig2D <- mNormal.VAS %>%
    emmeans(~ MC_resp*VAS_conf, at=list(VAS_conf=seq(0, 1, 0.01)),
            dpar='sigma', type='response') %>%
    as.data.frame %>%
    ggplot(aes(x=VAS_conf, y=response,
               ymin=lower.HPD, ymax=upper.HPD,
               color=MC_resp, fill=MC_resp)) +
    geom_line(size=2) + geom_ribbon(alpha=0.3) +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
    xlab('Confidence') + ylab('Estimated Standard Deviation\nof Causal Judgments') + theme_GC()



figConfMarginal + plot_spacer() + figConfMarginal + plot_spacer() +
    fig2A + figRatingMarginal + fig2B + figPredictionsMarginal +
        fig2C + plot_spacer() + fig2D + plot_spacer() +
        plot_layout(nrow=3, heights=c(0.3, 1, 1),
                    ncol=4, widths=c(1, 0.3, 1, 0.3),
                    guides='collect') +
        plot_annotation(tag_levels=list(c('', '', 'A', '', 'B', '', 'C', 'D'))) &
        theme(plot.tag=element_text(size=28),
              plot.tag.position=c(0.05, 1.1),
              legend.position='bottom')
ggsave('FigureS1-11.png', width=12.5, height=11)



#########################################################################################
##
## Test for norm effects
## 

ROPE.MC_resp <- .05
ROPE.MC_conf <- sd(judgments$MC_conf) * .1
ROPE.VAS_conf <- sd(judgments$VAS_conf) * .1

### Test for normality effect in continuous causal judgment / confidence
mNorm_VAS <- brm(bf(VAS_resp ~ structure*normality + (1 |v| vignette),
                    sigma ~ structure*normality + (1 |v| vignette)) +
                 bf(VAS_conf ~ structure*normality + (1 |v| vignette),
                    sigma ~ structure*normality + (1 |v| vignette)) +
                 set_rescor(TRUE),
                 prior=c(set_prior('normal(0, 1.0)', resp='VASresp'),
                         set_prior('normal(0, 5.0)', resp='VASresp', dpar='sigma'),
                         set_prior('normal(0, 1.0)', resp='VASconf'),
                         set_prior('normal(0, 5.0)', resp='VASconf', dpar='sigma')),
                 data=judgments, file='mNorm_VAS', iter=5000, backend='cmdstanr',
                 control=list(adapt_delta=.95), cores=4, inits="0", sample_prior="yes", save_pars=save_pars(all=TRUE))
summary(mNorm_VAS, priors=TRUE)

draws_VAS <- judgments %>%
    data_grid(structure, normality) %>%
    add_epred_draws(mNorm_VAS, re_formula=NA, dpar='sigma')

## contrast abnormal - normal
draws_VAS %>% compare_levels(.epred, by='normality') %>%
    pivot_wider(names_from=c(.category, structure), values_from=.epred) %>%
    select(VASresp_JC, VASresp_OD) %>%
    describe_posterior(ci=0.95, rope_ci=0.95, rope_range=c(-ROPE, ROPE)) %>%
    bind_rows(draws_VAS %>% compare_levels(.epred, by='normality') %>%
              pivot_wider(names_from=c(.category, structure), values_from=.epred) %>%
              select(VASconf_JC, VASconf_OD) %>%
              describe_posterior(ci=0.95, rope_ci=0.95, rope_range=c(-ROPE.VAS_conf, ROPE.VAS_conf))) %>%
    mutate(P=pd_to_p(pd))

## plot results
ggplot(draws_VAS, aes(x=structure, y=.epred, fill=normality)) +
    stat_halfeye(aes(side=ifelse(normality=='normal', 'left', 'right')),
                 point_interval=median_hdi, position=position_dodge(.25)) +
    theme_classic(base_size=14) + theme(axis.title.x=element_blank()) +
    facet_grid(.category~., scales='free_x',
               labeller=labeller(.category=c('VASresp'='Causal Judgment', 'VASconf'='Confidence'))) +
    ylim(0, 1) + ylab('Posterior Mean') +
    scale_x_discrete(labels=c('Conjunctive', 'Disjunctive')) + 
    scale_fill_manual(name='Normality', values=rev(wes_palette("Darjeeling1", n=2)),
                      labels=c('Normal', 'Abnormal'))
ggsave('FigureS1-13.png', width=7.5, height=5)





### Test for normality effect in discrete causal judgment / confidence
judgments <- judgments %>%
    mutate(MC_resp=factor(MC_resp, levels=c('did not cause', 'partially caused', 'totally caused'),
                          ordered=TRUE))
mNorm_MC <- brm(bf(MC_resp ~ structure*normality + (1 |v| vignette), family=cumulative) +
                bf(MC_conf ~ structure*normality + (1 |v| vignette),
                   sigma ~ structure*normality + (1 |v| vignette)),
                prior=c(set_prior('normal(0, 1.0)', resp='MCresp'),
                        set_prior('normal(0, 1.0)', resp='MCconf'),
                        set_prior('normal(0, 5.0)', resp='MCconf', dpar='sigma')),
                data=judgments, file='mNorm_MC3', iter=5000, backend='cmdstanr',
                control=list(adapt_delta=.95), inits="0",
                sample_prior="yes", save_pars=save_pars(all=TRUE), cores=4)
summary(mNorm_MC, prior=TRUE)

## significance testing
judgments %>%
    data_grid(structure, normality) %>%
    add_linpred_draws(mNorm_MC, re_formula=NA) %>%
    compare_levels(.linpred, by='normality') %>%
    pivot_wider(names_from=c(.category, structure), values_from=.linpred) %>%
    select(starts_with('MCresp')) %>%
    describe_posterior(ci=0.95, rope_ci=0.95, rope_range=c(-ROPE.MC_resp, ROPE.MC_resp)) %>%
    bind_rows(judgments %>%
              data_grid(structure, normality) %>%
              add_linpred_draws(mNorm_MC, re_formula=NA) %>%
              compare_levels(.linpred, by='normality') %>%
              pivot_wider(names_from=c(.category, structure), values_from=.linpred) %>%
              select(starts_with('MCconf')) %>%
              describe_posterior(ci=0.95, rope_ci=0.95, rope_range=c(-ROPE.MC_conf, ROPE.MC_conf))) %>%
    mutate(P=pd_to_p(pd))


## plot model results
draws_MC <- judgments %>%
    data_grid(structure, normality) %>%
    add_epred_draws(mNorm_MC, re_formula=NA) %>%
    mutate(.category=factor(.category,
                            levels=c('did not cause', 'partially caused', 'totally caused', 'MCconf'),
                            labels=c('P(Did Not Cause)', 'P(Partially Caused)', 'P(Totally Caused)', 'Confidence')),
           .resp=ifelse(.category=='Confidence', 'Confidence', 'Causal Judgment'))
p.conf <- draws_MC %>%
    filter(.category == 'Confidence') %>%
    ggplot(aes(x=structure, y=.epred, fill=normality)) +
    stat_halfeye(aes(side=ifelse(normality=='normal', 'left', 'right')),
                 point_interval=median_hdi, position=position_dodge(.25)) +
    theme_classic(base_size=18) + theme(axis.title.x=element_blank()) +
    ylim(0, 1) + ylab('Posterior Mean Confidence') +
    scale_x_discrete(name='Causal Structure', labels=c('Conjunctive', 'Disjunctive'),
                     expand=expansion(0, 0)) +
    scale_fill_manual(name='Normality', values=rev(wes_palette("Darjeeling1", n=2)),
                      labels=c('Normal', 'Abnormal'))
p.cause <- draws_MC %>%
    filter(.category != 'Confidence', .category != 'P(Totally Caused)') %>%
    pivot_wider(names_from=.category, values_from=.epred) %>%
    mutate(`P(Partially Caused)`=`P(Did Not Cause)`+`P(Partially Caused)`) %>%
    pivot_longer(`P(Did Not Cause)`:`P(Partially Caused)`, names_to='.category', values_to='.epred') %>%
    ggplot(aes(x=normality, y=.epred)) +
    stat_cdfinterval(aes(fill = stat(f)), thickness=1,
                     n=1000, point_interval=NULL) +
    stat_pointinterval(aes(group=.category)) +
    scale_fill_gradient2(low=PALETTE[2], mid=PALETTE[1], high=PALETTE[3], midpoint=0.5,
                         guide='legend', limits=c(0, 1), breaks=c(0, 0.5, 1),
                         name='Discrete\nCausal\nJudgment',
                         labels=c('Did Not Cause', 'Partially Caused', 'Totally Caused')) +
    facet_wrap(~ structure, labeller=labeller(structure=c('JC'='Conjunctive', 'OD'='Disjunctive'))) +
    scale_x_discrete(name='Normality', labels=c('Normal', 'Abnormal')) +
    ylab('Cumulative Probability\nof Response') + ylim(0, 1) +
    theme_classic(base_size=18)

p.cause / p.conf
ggsave('FigureS1-14.png', width=10, height=10, dpi=2000)




### Run a mediation analysis

## fit the null model (average total effect)
m <- brm(bf(VAS_resp ~ structure*normality + (1 |i| id) + (1 |v| vignette),
            sigma ~ structure*normality + (1 |i| id) + (1 |v| vignette)),
         prior=c(set_prior('normal(0, 1.0)'),
                 set_prior('normal(0, 5.0)', dpar='sigma')),
         save_pars=save_pars(all=TRUE), data=judgments %>% mutate(MC_resp=factor(MC_resp, ordered=FALSE)),
         cores=4, iter=5000, inits="0", file='mMed_null',
         backend='cmdstanr', control=list(adapt_delta=.99))

## fit the mediation models (average direct effect) for confidence in continuous/discrete judgment
m_med_MC <- brm(bf(VAS_resp ~ structure*normality + MC_conf:MC_resp + (1 |i| id) + (1 |v| vignette),
                   sigma ~ structure*normality + MC_conf:MC_resp + (1 |i| id) + (1 |v| vignette)),
         prior=c(set_prior('normal(0, 1.0)'),
                 set_prior('normal(0, 5.0)', dpar='sigma')),
         save_pars=save_pars(all=TRUE), data=judgments %>% mutate(MC_resp=factor(MC_resp, ordered=FALSE)),
         cores=4, iter=5000, inits="0", file='mMed_MC',
         backend='cmdstanr', control=list(adapt_delta=.99))
m_med_VAS <- brm(bf(VAS_resp ~ structure*normality + VAS_conf:MC_resp + (1 |i| id) + (1 |v| vignette),
                    sigma ~ structure*normality + VAS_conf:MC_resp + (1 |i| id) + (1 |v| vignette)),
         prior=c(set_prior('normal(0, 1.0)'),
                 set_prior('normal(0, 5.0)', dpar='sigma')),
         save_pars=save_pars(all=TRUE), data=judgments %>% mutate(MC_resp=factor(MC_resp, ordered=FALSE)),
         cores=4, iter=5000, inits="0", file='mMed_VAS',
         backend='cmdstanr', control=list(adapt_delta=.99))

## calculate prop. mediated (note: use ROPE = [-.1, .1] for prop. mediated)
m %>% emmeans(~normality|structure) %>% contrast(method='revpairwise') %>% gather_emmeans_draws(value='total') %>%
    full_join(m_med_MC %>% emmeans(~normality|structure) %>% contrast(method='revpairwise') %>% gather_emmeans_draws(value='mediation')) %>%
    mutate(prop_mediated=1-mediation/total) %>%
    pivot_wider(names_from=structure, values_from=total:prop_mediated) %>%
    select(-.chain, -.iteration, -.draw) %>%
    describe_posterior(ci=.95, rope_ci=.95, rope_range=c(-ROPE, ROPE)) %>%
    mutate(P=pd_to_p(pd))
m %>% emmeans(~normality|structure) %>% contrast(method='revpairwise') %>% gather_emmeans_draws(value='total') %>%
    full_join(m_med_VAS %>% emmeans(~normality|structure) %>% contrast(method='revpairwise') %>% gather_emmeans_draws(value='mediation')) %>%
    mutate(prop_mediated=1-mediation/total) %>%
    pivot_wider(names_from=structure, values_from=total:prop_mediated) %>%
    select(-.chain, -.iteration, -.draw) %>%
    describe_posterior(ci=.95, rope_ci=.95, rope_range=c(-ROPE, ROPE)) %>%
    mutate(P=pd_to_p(pd))
