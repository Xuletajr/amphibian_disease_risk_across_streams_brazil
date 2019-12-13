#############################################################################
#####  Bd occupancy, prevalence and infection intensity in the Brazil's #####   
#####         Atlantic Rainforest---- Output analysis and plots         #####
#############################################################################

# Clear memory
rm(list=ls(all=TRUE))

# Check the filepath for the working directory
getwd()                                                              

# Load model output
load("output_bd.RData")

# Load the jagsUI library
library(jagsUI)

# Load ggplot2 package
library(ggplot2)

# Look at output
out_bd

# Look at traceplots
traceplot(out_bd)

######### ######### ######### ######### ######### 
######### Occupancy model output
######### ######### ######### ######### ######### 

### Number of occupied sites
out_bd$mean$mean.occ; out_bd$q2.5$mean.occ; out_bd$q97.5$mean.occ; out_bd$sd$mean.occ

hist(out_bd$sims.list$mean.occ, breaks=25, col="grey70", ylim=c(0,300), xlab="Number of occupied sites by Bd")
abline(v=out_bd$mean$mean.occ, lwd = 3, lty=2, col="red")

dat.occ <- data.frame (p = out_bd$sims.list$mean.occ)

#####
### Occupancy probability
round(mean(plogis(out_bd$sims.list$a0)),  3)
round(quantile(plogis(out_bd$sims.list$a0),prob=0.025),  3)
round(quantile(plogis(out_bd$sims.list$a0),prob=0.975), 3)

# Stream density effect
round(out_bd$mean$a2, 2); round(out_bd$q2.5$a2, 2); round(out_bd$q97.5$a2, 2)

# Proportion of stream density parameter posterior distribution including positive values
length(which(out_bd$sims.list$a2 > 0))/3000

# Amphibian richness effect
round(out_bd$mean$a3, 2); round(out_bd$q2.5$a3, 2); round(out_bd$q97.5$a3, 2)

# Proportion of richness parameter posterior distribution including negative values
length(which(out_bd$sims.list$a3 < 0))/3000

# Forest effect
round(out_bd$mean$a1, 2); round(out_bd$q2.5$a1, 2); round(out_bd$q97.5$a1, 2)

# Proportion of forest cover parameter posterior distribution including positive values
length(which(out_bd$sims.list$a1 > 0))/3000 # 

######### ######### ######### ######### ######### ######### #########
######### Predict effect of covariates on occupancy probability
######### ######### ######### ######### ######### ######### #########

# Getting forest cover values for prediction
(original.for.pred <- seq(30, 100, length.out = 50))
(for.pred <- (original.for.pred - mforest)/sdforest)

# Getting stream density values for prediction
sort(hydro)
(original.hydro.pred <- seq(329, 1121, length.out = 50))
(hydro.pred <- (original.hydro.pred  - mhydro)/sdhydro)

# Getting host species richness values for prediction
sort(ric)
(original.ric.pred <- seq(7.96, 18.66, length.out = 50))
(ric.pred <- (original.ric.pred  - mric)/sdric)

###
# Like Kery and Royle (2015)
nsamp <- length(out_bd$sims.list$a0)
pred.occ <- array(NA, dim=c(50, nsamp, 3))
dim(pred.occ)

for(i in 1:nsamp){
  pred.occ[,i,1] <- plogis(out_bd$sims.list$a0[i] + out_bd$sims.list$a1[i] * for.pred)
  pred.occ[,i,2] <- plogis(out_bd$sims.list$a0[i] + out_bd$sims.list$a2[i] * hydro.pred )
  pred.occ[,i,3] <- plogis(out_bd$sims.list$a0[i] + out_bd$sims.list$a3[i] * ric.pred )
}

# Get postrior mean, 95% CIs and plot
(occ.mean <- apply(pred.occ, c(1,3), mean))
(occ.ci <- apply(pred.occ, c(1,3), function(x) quantile(x,prob=c(0.025,0.975))))

# Putting all data in data frame
dat2.1 <- data.frame(occ.mean = occ.mean[,2],
                   cov.seq = original.hydro.pred,
                   LL = occ.ci[1,,2],
                   UL = occ.ci[2,,2],
                   covariates = factor(rep("Stream density", 50 )))
                   ,
                   covariates  = factor(rep(1:3, each = 50), levels = 1:3, 
                                        labels = c("Stream density", "Forest cover", 
                                                   "Amphibian richness")))
dat2.1

dat2.2 <- data.frame(occ.mean = occ.mean[,1],
                     cov.seq = original.for.pred,
                     LL =  occ.ci[1,,1],
                     UL =  occ.ci[2,,1],
                     covariates = factor(rep("Forest cover", 50 )))
                     
dat2.2                                                     
                                                                                                        

dat2.3 <- data.frame(occ.mean = occ.mean[,3],
                     cov.seq =  original.ric.pred,
                     LL = occ.ci[1,,3],
                     UL = occ.ci[2,,3],
                     covariates  = factor(rep("Amphibian richness", 50 ))) 
                                         
dat2.3

# Preparing plot for stream density. 
(plot_stream <- ggplot(dat2.1, aes(cov.seq, occ.mean)) +
    geom_ribbon(aes(ymin = LL, ymax = UL, fill = covariates), 
                alpha =.3) +
    geom_line(aes(colour = covariates), size = 1) + 
    scale_colour_manual("", values = "black") +
    scale_fill_manual("", values = "gray50")  +
    theme_bw () +
    theme(legend.position = "none") +
    labs(x = expression(atop("Length of stream (m)", paste("within a buffer"))), 
         y = expression(paste(italic("Bd"), " occurrence probability"))) +
    theme(axis.text.x =  element_text(size = 8, angle = 45, colour = "black", hjust = 1), 
          axis.text.y = element_text(size = 10, colour = "black"), 
          panel.grid = element_blank(), 
          strip.text.x = element_text(size = 12, color = "black"),
          strip.background = element_rect(fill="gray75")) )

# Preparing plot for forest cover
(plot_forest <- ggplot(dat2.2, aes(cov.seq, occ.mean)) +
    geom_ribbon(aes(ymin = LL, ymax = UL, fill = covariates), 
                alpha =.3) +
    geom_line(aes(colour = covariates), size = 1) + 
    scale_colour_manual("", values = "black") +
    scale_fill_manual("", values = "gray50")  +
    theme_bw () +
    theme(legend.position = "none") +
    xlab("Forest cover (%)") +
    ylab(NULL) +
    theme(axis.text.x = element_text(size = 8, angle = 45, colour = "black", hjust = 1), 
          axis.text.y = element_blank(), 
          panel.grid = element_blank(), 
          strip.text.x = element_text(size = 12, color = "black"),
          strip.background = element_rect(fill="gray75")) )

# Preparing plot for frog species richness
(plot_richness <- ggplot(dat2.3, aes(cov.seq, occ.mean)) +
    geom_ribbon(aes(ymin = LL, ymax = UL, fill = covariates), 
                alpha =.3) +
    geom_line(aes(colour = covariates), size = 1) + 
    scale_colour_manual("", values = "black") +
    scale_fill_manual("", values = "gray50")  +
    theme_bw () +
    theme(legend.position = "none") +
    xlab("Number of species") +
    ylab(NULL) +
    theme(axis.text.x = element_text(size = 8, angle = 45, colour = "black", hjust = 1), 
          axis.text.y = element_blank(), 
          panel.grid = element_blank(), 
          strip.text.x = element_text(size = 12, color = "black"),
          strip.background = element_rect(fill="gray75")) ) 

# Load ggplot2 package
library(ggpubr)

# Figure 2 Draft
tiff(
  "./output_figures/fig02.tiff",
  width     = 6.5,
  height    = 3,
  units     = "in",
  res       = 900,
  pointsize = 5.5
)

# Grouping the 3 plots in only one figure
ggarrange(plot_stream, plot_forest, plot_richness , 
          ncol = 3, nrow = 1, align = "hv", hjust = 100)

dev.off()

######### ######### ######### ######### ######### ######### #########
######### Infection intensity
######### ######### ######### ######### ######### ######### #########
# n1 parameter - forest cover
out_bd$mean$n1; out_bd$q2.5$n1; out_bd$q97.5$n1
length(which(out_bd$sims.list$n1 > 0))/3000

# n2 parameter - species richness
out_bd$mean$n2; out_bd$q2.5$n2; out_bd$q97.5$n2
length(which(out_bd$sims.list$n2 > 0))/3000

nsamp <- length(out_bd$sims.list$a0)
pred.int <- array(NA, dim=c(50, nsamp, 2))
dim(pred.int)

for(i in 1:nsamp){
  pred.int[,i,1] <- (out_bd$sims.list$n0[i] + out_bd$sims.list$n1[i] * for.pred)
  pred.int[,i,2] <- (out_bd$sims.list$n0[i] + out_bd$sims.list$n2[i] * hydro.pred )
}

# Get postrior mean, 95% CIs and plot
(int.mean <- apply(pred.int, c(1,3), mean))
(int.ci <- apply(pred.int, c(1,3), function(x) quantile(x,prob=c(0.025,0.975))))

# Putting all data in data frame
dat6 <- data.frame(int.mean = c(int.mean[,2], int.mean[,1]),
                   cov.seq = c(original.for.pred, 
                               original.ric.pred),
                   LL = c(int.ci[1,,2], int.ci[1,,1]),
                   UL = c(int.ci[2,,2], int.ci[2,,1]),
                   covariates  = factor(rep(1:2, each = 50), levels = 1:2, 
                                        labels = c("Forest cover", 
                                                   "Stream density")))


# Figure 6 - Figure was not included in the manuscript
tiff(
  "./output_figures/fig06.tiff",
  width     = 6.5,
  height    = 3,
  units     = "in",
  res       = 900,
  pointsize = 5.5
)

ggplot(dat6, aes(cov.seq, int.mean)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = covariates), alpha =.3) +
  geom_line(aes(colour = covariates), size = 1) + 
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5),
                     label = c(1, 
                               2.71, 
                               7.38, 
                               20.08, 
                               54.59, 
                               148.41))+
  scale_colour_manual("", values=c("black", "black")) +
  scale_fill_manual("", values=c("gray50", "gray50"))  +
  facet_wrap( ~ covariates, nrow = 1, scales="free_x", strip.position =) + 
  theme_bw () +
  theme(legend.position="none") +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.title = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x=element_text(size = 12, color = "black"),
        strip.background = element_rect( fill="gray75"),
        axis.title=element_text(size = 12))

dev.off()

######### ######### ######### ######### ######### ######### #########
######### Prevalence model output
######### ######### ######### ######### ######### ######### #########
# Naive Bd prevalence within a stream
bd.intensity1 <- read.table("./data/bd_load.txt", header=T)
bd.intensity1[bd.intensity1 > 0] <- 1

frog_swabbed <- bd.intensity1
frog_swabbed[frog_swabbed==0] <- 1
n.swabbed <- rowSums (frog_swabbed, na.rm = TRUE) 

prev <- rowSums(bd.intensity1, na.rm = TRUE)/n.swabbed

min(prev)
max(prev)
mean(prev)

# Prevalence of terrestrial-breeding species (AI-0) 
round(out_bd$mean$prev.AI0,  3); round(mean(plogis(out_bd$sims.list$b0 )),  3) 
round(out_bd$q2.5$prev.AI0,  3) 
round(out_bd$q97.5$prev.AI0, 3)

# Prevalence (or detection probability ) of arboreal species with aquatic larvae (AI-1) 
round(out_bd$mean$prev.AI1,  3)
round(out_bd$q2.5$prev.AI1,  3)
round(out_bd$q97.5$prev.AI1, 3)

# Prevalence (or detection probability ) of terrestrial species with aquatic larvae (AI-2) 
round(out_bd$mean$prev.AI2,  3)
round(out_bd$q2.5$prev.AI2,  3)
round(out_bd$q97.5$prev.AI2, 3)

# Infection intensity of terrestrial-breeding species (AI-0) 
out_bd$mean$int.AI0; out_bd$q2.5$int.AI0; out_bd$q97.5$int.AI0

# Infection intensity of arboreal species with aquatic larvae (AI-1) 
out_bd$mean$int.AI1; out_bd$q2.5$int.AI1; out_bd$q97.5$int.AI1

# Infection intensity of terrestrial species with aquatic larvae (AI-2)
out_bd$mean$int.AI2; out_bd$q2.5$int.AI2; out_bd$q97.5$int.AI2  

# The differences between parameters related to aquatic index at each MCMC iteration 
# following Ruiz-Guti?rrez et al. (2010). We computed the proportion of iterations 
# where one parameter was greater than the other. We report these values as
# the probability of parameter x > parameter z, written as Pr(x > z). 

# Difference in the prevalence parameter between AI
# Prevalence: AI-0 vs AI-1 - Pr(b0 > b1)
mean(out_bd$sims.list$b0 > out_bd$sims.list$b1)
length(which(out_bd$sims.list$b0 - out_bd$sims.list$b1 > 0))/3000

# Prevalence: AI-0 vs AI-2 - Pr(b0 > b2)
mean(out_bd$sims.list$b0 > out_bd$sims.list$b2)
length(which(out_bd$sims.list$b0 - out_bd$sims.list$b2 > 0))/3000

# Prevalence: AI-1 vs AI-2 - Pr(b2 > b1)
mean(out_bd$sims.list$b2 > out_bd$sims.list$b1)
length(which(out_bd$sims.list$b2 - out_bd$sims.list$b1 > 0))/3000

# Difference in the intensity parameter between AI
# Infection intensity: AI-0 vs AI-1 - Pr(b0 > b1)
mean(out_bd$sims.list$int.AI0 < out_bd$sims.list$int.AI1)
length(which(out_bd$sims.list$int.AI0 - out_bd$sims.list$int.AI1 > 0))/3000

# Infection intensity: AI-0 vs AI-2 - Pr(b0 > b2)
mean(out_bd$sims.list$int.AI0 < out_bd$sims.list$int.AI2)
length(which(out_bd$sims.list$int.AI0 - out_bd$sims.list$int.AI2 > 0))/3000

# Infection intensity: AI-2 vs AI-1 - Pr(b2 > b1)
mean(out_bd$sims.list$int.AI2 > out_bd$sims.list$int.AI1)
length(which(out_bd$sims.list$int.AI2 - out_bd$sims.list$int.AI1 > 0))/3000

### FIGURE 3 MANUSCRIPT
# Putting prevalence data in a data frame
dat3 <- data.frame (p = c(out_bd$sims.list$prev.AI0, 
                          out_bd$sims.list$prev.AI1, 
                          out_bd$sims.list$prev.AI2),
                    AI = c(rep("AI-0", length(out_bd$sims.list$prev.AI0)), 
                           rep("AI-1", length(out_bd$sims.list$prev.AI1)), 
                           rep("AI-2", length(out_bd$sims.list$prev.AI2))))

dat3

# Exporting figure 3
tiff(
  "./output_figures/fig03.tiff",
  width     = 4,
  height    = 5,
  units     = "in",
  res       = 900,
  pointsize = 4
) 

ggplot(dat3, aes(x = p, fill = AI, color = AI, group = AI)) + 
  geom_density( stat = "density", alpha = .5 ,size = .8, linetype = "dashed")+
  scale_fill_manual(values = c( "#E69F00" , "#009E73", "#56B4E9"))+
  scale_color_manual(values = c( "#E69F00" , "#009E73", "#56B4E9"))+
  theme_classic() +
  labs(y = "Density function", 
       x = expression(paste(italic("Bd"), " infection prevalence"))) +
  theme(axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=rel(1.5)),
        legend.key.size =  unit(0.25, "in"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14))

dev.off()

# Forest effect - prevalence
round(out_bd$mean$b3, 3); round(out_bd$q2.5$b3, 3); round(out_bd$q97.5$b3, 3)

# round(out_bd$q25$b3, 2); round(out_bd$q75$b3, 2); round(out_bd$sd$b3, 2)
length(which(out_bd$sims.list$b3 > 0))/3000 # 

######### ######### ######### ######### ######### ######### 
#########  Probability of detecting Bd in a population
#########        Two scenarios of forest cover
######### ######### ######### ######### ######### ######### 
# 1 - (1 - p)^n  # n is the number of individuals tested

prob.detec.ai0.f100 <- matrix(NA, nrow = length(out_bd$sims.list$beta), ncol = 35)
prob.detec.ai1.f100 <- matrix(NA, nrow = length(out_bd$sims.list$beta), ncol = 35)
prob.detec.ai2.f100 <- matrix(NA, nrow = length(out_bd$sims.list$beta), ncol = 35)

prob.detec.ai0.f30 <- matrix(NA, nrow = length(out_bd$sims.list$beta), ncol = 35)
prob.detec.ai1.f30 <- matrix(NA, nrow = length(out_bd$sims.list$beta), ncol = 35)
prob.detec.ai2.f30 <- matrix(NA, nrow = length(out_bd$sims.list$beta), ncol = 35)


# Terrestrial-breeding (AI-0) - Forest cover 30%
for (i in 1:length(out_bd$sims.list$beta)) {
  for (j in 1:35) { 
    prob.detec.ai0.f30[i,j]  <- 1 - (1 - plogis(out_bd$sims.list$beta[i] + out_bd$sims.list$b0[i] + 
                                                  out_bd$sims.list$b3[i] * -1.6782 ))^j
  }
}

# Terrestrial-breeding (AI-0 ) - Forest cover 100%
for (i in 1:length(out_bd$sims.list$beta)) {
  for (j in 1:35) { 
    prob.detec.ai0.f100[i,j]  <- 1 - (1 - plogis(out_bd$sims.list$beta[i] + out_bd$sims.list$b0[i] + 
                                                   out_bd$sims.list$b3[i] * 1.4138 ))^j
  }
}


# Aquatic-breeding with arboreal habit (AI-1) - Forest cover 30%
for (i in 1:length(out_bd$sims.list$beta)) {
  for (j in 1:35) { 
    prob.detec.ai1.f30[i,j]  <- 1 - (1 - plogis(out_bd$sims.list$beta[i] + out_bd$sims.list$b1[i] +
                                                  out_bd$sims.list$b3[i] * -1.6782))^j
  }
}

# Aquatic-breeding with arboreal habit (AI-1) - Forest cover 100%
for (i in 1:length(out_bd$sims.list$beta)) {
  for (j in 1:35) { 
    prob.detec.ai1.f100[i,j]  <- 1 - (1 - plogis(out_bd$sims.list$beta[i] + out_bd$sims.list$b1[i] +
                                                   out_bd$sims.list$b3[i] * 1.4138))^j
  }
}

# Aquatic-breeding with terrestrial habit (AI-2) - Forest cover 30%
for (i in 1:length(out_bd$sims.list$beta)) {
  for (j in 1:35) { 
    prob.detec.ai2.f30[i,j]  <- 1 - (1 - plogis(out_bd$sims.list$beta[i] + out_bd$sims.list$b2[i] + 
                                                  out_bd$sims.list$b3[i] * -1.6782))^j
  }
}

# Aquatic-breeding with terrestrial habit (AI-2) - Forest cover 100%
for (i in 1:length(out_bd$sims.list$beta)) {
  for (j in 1:35) { 
    prob.detec.ai2.f100[i,j]  <- 1 - (1 - plogis(out_bd$sims.list$beta[i] + out_bd$sims.list$b2[i] + 
                                                   out_bd$sims.list$b3[i] * 1.4138))^j
  }
}

# Data frame for 30% forest cover
dat5a <- data.frame ( y = c(apply(prob.detec.ai0.f30, 2 , mean), 
                            apply(prob.detec.ai1.f30, 2 , mean), 
                            apply(prob.detec.ai2.f30, 2 , mean)),
                      samples = rep(1:35, 3),
                      group = rep (c ("AI-0", "AI-1", "AI-2"), each=35))

# Data frame for 100% forest cover
dat5b <- data.frame ( y = c(apply(prob.detec.ai0.f100, 2 , mean), 
                            apply(prob.detec.ai1.f100, 2 , mean), 
                            apply(prob.detec.ai2.f100, 2 , mean)),
                      samples = rep(1:35, 3),
                      group = rep (c ("AI-0", "AI-1", "AI-2"), each=35))


ai.colour <-c( "#E69F00" , "#009E73", "#56B4E9")

# Figure S2 Draft (Appendix S2)
tiff(
  "fig5ea.tiff",
  width     = 5,
  height    = 4,
  units     = "in",
  res       = 600,
  pointsize = 5.5
)

ggplot(dat5a, aes(x=samples, y=y,  colour=group, group=group,  pch=group)) +
  geom_line(alpha = 0.9) +
  geom_point(size = 2.5,  alpha = 0.8) +
  scale_fill_manual(values = c( "#E69F00" , "#009E73", "#56B4E9"))+
  scale_color_manual(values = c( "#E69F00" , "#009E73", "#56B4E9"))+
  labs(x = "", y = "") +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0,1)) +
  scale_x_continuous(breaks = seq(0, 35, 5), limits = c(0,35)) +
  labs(shape="Aquatic index") +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  theme_bw( ) +
  #ggtitle("Forest cover - 30%")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=14, face="bold", hjust = 0.5))

dev.off()

tiff(
  "fig5eb.tiff",
  width     = 5,
  height    = 4,
  units     = "in",
  res       = 600,
  pointsize = 5.5
)


ggplot(dat5b, aes(x=samples, y=y, colour=group, group=group,  pch=group)) +
  geom_line(alpha = 0.9) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_fill_manual(values = c( "#E69F00" , "#009E73", "#56B4E9"))+
  scale_color_manual(values = c( "#E69F00" , "#009E73", "#56B4E9"))+
  labs(x = "", y = "") +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0,1)) +
  scale_x_continuous(breaks = seq(0, 35, 5), limits = c(0,35)) +
  labs(shape="Aquatic index") +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  theme_bw( ) +
  #ggtitle("Forest cover - 100%")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=14, face="bold", hjust = 0.5))


dev.off()

######   end----