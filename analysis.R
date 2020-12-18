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
judgments <- read.csv('data/processed_data.csv', header=TRUE)
judgments$MC_conf <- judgments$MC_conf / 100
judgments$VAS_resp <- judgments$VAS_resp / 100
judgments$VAS_conf <- judgments$VAS_conf / 100
judgments <- judgments %>%
    mutate(MC_resp=factor(MC_resp, levels=c('partially caused',
                                            'did not cause',
                                            'totally caused')))

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


## Fit heteroscedastic Gaussian model
mNormal <- brm(bf(VAS_resp ~ MC_resp * MC_conf + (1 |i| id) +
                      (1 |v| vignette:structure:condition),
                  sigma ~ MC_resp * MC_conf + (1 |i| id) +
                      (1 |v| vignette:structure:condition)),
               prior=c(set_prior('normal(0, 1.0)'),
                       set_prior('normal(0, 5.0)', dpar='sigma')),
               save_all_pars=TRUE, data=judgments,
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
fig2A
ggsave('data.png', width=12.5, height=10)

## plot model predictions
predictions <- add_predicted_draws(judgments, mNormal) %>% ungroup

fig2B <- predictions %>%
    filter(.draw %in% sample(min(.draw):max(.draw), 50)) %>%
    dplot(y='.prediction', ylab='Predicted Causal Rating', bw.x=0.15, legend=FALSE)

figPredictionsMarginal <- marginal_density(predictions, '.prediction') +
    coord_flip()

fig2C <- mNormal %>%
    emmeans(~ MC_resp*MC_conf, at=list(MC_conf=seq(0, 1, 0.01))) %>%
    as.data.frame %>%
    mplot(y='emmean', ymin='lower.HPD', ymax='upper.HPD',
          ylab='Estimated Mean Causal Rating')

fig2D <- mNormal %>%
    emmeans(~ MC_resp*MC_conf, at=list(MC_conf=seq(0, 1, 0.01)),
            dpar='sigma', type='response') %>%
    as.data.frame %>%
    mplot(y='response', ymin='lower.HPD', ymax='upper.HPD',
          ylab='Estimated Standard Deviation\nof Causal Ratings')


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

ggsave('Figure2.png', width=12.5, height=11)


for (c in 1:100) {
    predictions %>%
        filter(MC_conf <= c/100) %>%
        marginal_density('.prediction', void=FALSE, legend='bottom') +
        ggtitle('Causal Rating')
    ggsave(sprintf('plots/hist-%03d.png', c), width=10, height=5)
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
    mutate(Parameter=paste0('Mean, ', Parameter))
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
    ylab('Mean Causal Rating Contrasts: Discrete Rating') +
    theme_classic() + theme(legend.position='none')
ggsave('contrast-discrete.png', width=6, height=4)

c2.mean <- em %>% contrast('trt.vs.ctrl', simple='MC_conf') %>%
    describe_posterior(ci=0.95, rope_ci=0.95, rope_range=c(-ROPE, ROPE)) %>%
    mutate(Parameter=paste0('Mean, ', Parameter))
c2.mean

em %>% contrast('trt.vs.ctrl', simple='MC_conf') %>% gather_emmeans_draws %>%
    ggplot(aes(x=MC_resp, y=.value,
               fill=stat(ifelse(abs(y) < ROPE, '0', group)))) +
    stat_eye(position=position_dodge(width=1), n=10000) +
    geom_hline(yintercept=c(-ROPE, ROPE), linetype='dashed') +
    scale_fill_manual(values=c('gray80', PALETTE)) +
    ylab('Mean Causal Rating Contrasts: Confidence') +
    xlab('Discrete Rating') +
    theme_classic() + theme(legend.position='none')
ggsave('contrast-confidence.png', width=6, height=4)

c3.mean <- em %>% contrast('trt.vs.ctrl', interaction=TRUE) %>%
    describe_posterior(ci=0.95, rope_ci=0.95, rope_range=c(-ROPE, ROPE)) %>%
    mutate(Parameter=paste0('Mean, ', Parameter))
c3.mean



em %>% contrast('trt.vs.ctrl', interaction=TRUE) %>% gather_emmeans_draws %>%
    rename(MC_resp=MC_resp_trt.vs.ctrl,
           MC_conf=MC_conf_trt.vs.ctrl) %>%
    ggplot(aes(x=MC_resp, y=.value,
               fill=stat(ifelse(abs(y) < ROPE, '0', group)))) +
    stat_eye(position=position_dodge(width=1), n=10000) +
    geom_hline(yintercept=c(-ROPE, ROPE), linetype='dashed') +
    scale_fill_manual(values=c('gray80', PALETTE[-1])) +
    ylab('Mean Causal Rating Contrasts:\nConfidence x Discrete Rating') +
    xlab('Discrete Rating') +
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
    mutate(Parameter=paste0('SD, ', Parameter))
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
    ylab('Standard Deviation of Causal Rating Contrasts:\nDiscrete Rating') +
    theme_classic() + theme(legend.position='none')
ggsave('contrast-discrete-sd.png', width=6, height=4)



c2.sd <- em.sd %>% contrast('trt.vs.ctrl', simple='MC_conf') %>%
    describe_posterior(ci=0.95, rope_ci=0.95, rope_range=c(-sdROPE, sdROPE)) %>%
    mutate(Parameter=paste0('SD, ', Parameter))
c2.sd

em.sd %>% contrast('trt.vs.ctrl', simple='MC_conf') %>%
    gather_emmeans_draws('sigma') %>%
    ggplot(aes(x=MC_resp, y=sigma,
               fill=stat(ifelse(abs(y) < ROPE, '0', group)))) +
    stat_eye(position=position_dodge(width=1), n=10000) +
    geom_hline(yintercept=c(-ROPE, ROPE), linetype='dashed') +
    scale_fill_manual(values=c('gray80', PALETTE)) +
    ylab('Standard Deviation of Causal Rating Contrasts\nConfidence') +
    xlab('Discrete Rating') +
    theme_classic() + theme(legend.position='none')
ggsave('contrast-confidence-sd.png', width=6, height=4)

c3.sd <- em.sd %>% contrast('trt.vs.ctrl', interaction=TRUE) %>%
    describe_posterior(ci=0.95, rope_ci=0.95,
                       rope_range=c(-sdROPE, sdROPE)) %>%
    mutate(Parameter=paste0('SD, ', Parameter))
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
    ylab('Standard Deviation of Causal Rating Contrasts:\nConfidence x Discrete Rating') +
    xlab('Discrete Rating') +
    theme_classic() + theme(legend.position='none')
ggsave('contrast-discreteXconfidence-sd.png', width=6, height=4)



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
           Median, CI, pd, ROPE, ROPE_Percentage) %>%
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
                                            'Discrete Rating [did not cause - partially caused]',
                                        'MC_resptotallycaused'=
                                            'Discrete Rating [totally caused - partially caused]',
                                        'MC_conf'='Confidence')),
                      levels=c('Discrete Rating [totally caused - partially caused] : Confidence',
                               'Discrete Rating [did not cause - partially caused] : Confidence',
                               'Confidence',
                               'Discrete Rating [totally caused - partially caused]',
                               'Discrete Rating [did not cause - partially caused]',
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
                                            'Discrete Rating [did not cause - partially caused]',
                                        'MC_resptotallycaused'=
                                            'Discrete Rating [totally caused - partially caused]',
                                        'MC_conf'='Confidence')),
                      levels=c('Discrete Rating [totally caused - partially caused] : Confidence',
                               'Discrete Rating [did not cause - partially caused] : Confidence',
                               'Confidence',
                               'Discrete Rating [totally caused - partially caused]',
                               'Discrete Rating [did not cause - partially caused]',
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
                       rope_range=c(-sd(unlist(.))*0.1, sd(unlist(.))*0.1))




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
    xlab('Confidence') + ylab('Mean Causal Rating') +
    theme_GC() + theme(panel.spacing=unit(1.5, "lines"),
                       legend.position='bottom')
ggsave('data_vignette.png', width=10, height=10)
