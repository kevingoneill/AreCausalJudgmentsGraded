library(RColorBrewer)
library(wesanderson)
library(ggisoband)
library(viridis)
library(patchwork)

## A custom color palette
PALETTE <- c(wes_palette("Darjeeling1", n=2)[2],
             wes_palette("Darjeeling1", n=3)[3], wes_palette("Darjeeling1", n=1))

## A custom ggplot theme
theme_GC <- function(palette=PALETTE) {
    list(theme_classic(),
         theme(axis.title.x=element_text(size=20, margin=margin(t=10)),
               axis.title.y=element_text(size=20, margin=margin(r=10)),
               axis.text.x=element_text(size=16, margin=margin(t=5)),
               axis.text.y=element_text(size=16, margin=margin(r=5)),
               legend.text=element_text(size=15),
               legend.title=element_text(size=20)),
         scale_alpha_continuous(range = c(0, 1)),
         scale_fill_manual(values=palette, name='Discrete Causal Judgment'),
         scale_color_manual(values=palette, name='Discrete Causal Judgment'))
}

## A function to plot the marginal density of judgments$var,
##   grouping by judgments$group.
marginal_density <- function(judgments, var, group='MC_resp',
                             fill='white', void=TRUE, legend='none',
                             palette=PALETTE) {
    judgments <- judgments %>%
        mutate(!!var := pmax(0.0, pmin(1.0, pull(., var))))
    p <- ggplot(judgments) + aes_string(x=var) + xlab('') + ylab('')
    if (is.character(group))
        if (legend == 'none')
            p <- p + aes_string(color=group, fill=group) +
                scale_fill_manual(values=palette) +
                scale_color_manual(values=palette) +
                geom_density(bw=0.01, alpha=0.33, show.legend=FALSE)
        else
            p <- p + aes_string(color=group, fill=group) +
                scale_fill_manual(values=palette) +
                scale_color_manual(values=palette) +
                geom_density(bw=0.01, alpha=0.33)
    else
        p <- p + geom_density(bw=0.01, fill=fill)

    p <- p + xlim(0.0, 1.0) + theme_classic() +
        guides(fill=guide_legend(reverse=T),
               color=guide_legend(reverse=T))
    
    if (void)
        p <- p + theme_void()
    else
        p <- p + theme(axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank(),
                       axis.line.y=element_blank(),
                       axis.text.x=element_blank(),
                       plot.title=element_text(hjust = 0.5))
    return(p + theme(legend.position=legend))
}

## Plot the multivariate distribution along with marginal distributions
density_2d <- function(center, top, right) {
    ## Arrange plots
    p <- top + plot_spacer() + center + right +
        plot_layout(nrow=2, ncol=2, widths=c(1, 0.3), heights=c(0.3, 1),
                    guides='collect')
    return(p)
}

dplot <- function(j, palette=PALETTE, y='VAS_resp', ylab='Causal Judgment',
                  bw.x=NULL, bw.y=NULL, legend=TRUE) {
    j <- j %>% ungroup %>%
        mutate(!!y := pmax(0.0, pmin(1.0, pull(., !!y)))) %>%
        ggplot(aes(x=MC_conf, group=MC_resp, color=MC_resp, fill=MC_resp)) +
        aes_string(y=y) +
        geom_density_bands(aes(alpha=stat(ndensity)), h=c(bw.x, bw.y),
                           bins=9, show.legend=c(alpha=FALSE,
                                                 fill=legend)) +
        coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
        xlab('Confidence') + ylab(ylab) + theme_GC(palette)
}

mplot <- function(samples, min=0.0, max=1.0,
                  palette=PALETTE,
                  y='.value', ymin=paste0(y, '.lower'),
                  ymax=paste0(y, '.upper'), ylab='Causal Judgment') {
    samples %>%
        filter(MC_conf >= min & MC_conf <= max) %>%
        ggplot() +
        aes(x=MC_conf, color=MC_resp, fill=MC_resp) +
        aes_string(y=y, ymin=ymin, ymax=ymax) +
        geom_line(size=2) + geom_ribbon(alpha=0.3) +
        coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
        xlab('Confidence') + ylab(ylab) + theme_GC(palette)
}

