#############################################################################
#####               Bd real data with covariates - Jose                  ####
#####    Bd occupancy and prevalence in the Brazil's Atlantic Forest     ####
#############################################################################

###############################################
###         Reading the data
###############################################

# Bd detection/non-detection matrix
bd.pres <- read.table("bd_pres.txt", header=T)   

# Bd intesity matrix
bd.intensity <- read.table("bd_load.txt", header=T)

# number of sites sampled
(nsites =  nrow(bd.pres))

# Number of amphibian individuals sampled by site (swab sampled)
frog_swabbed <- bd.pres
frog_swabbed[frog_swabbed==0] <- 1
n.swabbed <- rowSums (frog_swabbed, na.rm = TRUE) # different number of swas for sampling site

# Summary swabs statistics
sum(n.swabbed)  # Total number of swabs
mean(n.swabbed) # Mean for sampling site
sd(n.swabbed)   # standard deviation for sampling site

# Number of infected frogs
sum(bd.pres, na.rm=T)

# Remove 0s from intensity data and log transform
bd.intensity[bd.intensity == 0] <- NA
#bd.intensity <- log(bd.intensity)

# Total number of swabs collected
bd.pres2 <- ifelse(bd.pres == 0, 1, 1)
sum(bd.pres2, na.rm = TRUE)
################################################
###         Aquatic index
###############################################
# We designated a host aquatic index (AI) to each sampled anuran following 
# Becker et al. (2014). Partitioning the net effect of host diversity on an emerging amphibian pathogen
# Adapted from: Lips et al. (2003). Ecological traits predict amphibian population declines in Central America

# AI0 - Direct developers
bd.AI0 <- read.table("bd_AI0.txt", header=T) 
bd.AI0 <- as.matrix(bd.AI0)
bd.AI0[is.na(bd.AI0)] <- 0

# AI1 - Species breeding in aquatic habitats but occupying the arboreal stratum 
bd.AI1 <- read.table("bd_AI1.txt", header=T) 
bd.AI1 <- as.matrix(bd.AI1)
bd.AI1[is.na(bd.AI1)] <- 0

# AI2 - Species breeding in aquatic habitats and occupying the margins of streams and other bodies of water
bd.AI2 <- read.table("bd_AI2.txt", header=T)
bd.AI2 <- as.matrix(bd.AI2)
bd.AI2[is.na(bd.AI2)] <- 0

################################################
###         Julian date
###############################################
(date <- as.matrix(read.table("date_bd.txt", header = T)))
mdate <- mean(date, na.rm=TRUE)
sddate <- sd(date, na.rm=TRUE)
date1 <- (date-mdate) / sddate 
date2 <- date1 * date1
round(date1, digits=4)

################################################
###         Covariates
###############################################
# Read in the habitat data 
habitat <- read.table("habitat_covariates.txt", header=TRUE, na.strings=c("NA"))
head(habitat)

# Standardize the natural forest cover (buffer radius of 200 m)
forest <- as.vector (habitat$for200)
mforest <- mean(forest, na.rm=TRUE)
sdforest <- sd(forest, na.rm=TRUE)
forest1 <- as.vector( (forest-mforest) / sdforest )

# Standardize the stream density (buffer radius of 200 m)
hydro <- as.vector (habitat$hydro200_length)
mhydro <- mean(hydro, na.rm=TRUE)
sdhydro <- sd(hydro, na.rm=TRUE)
hydro1 <- as.vector( (hydro-mhydro) / sdhydro )

# Standardize the the standard deviation of slope (proxy of topography complexity)
slope <- as.vector (habitat$slope200_sd)
mslope <- mean(slope, na.rm=TRUE)
sdslope <- sd(slope, na.rm=TRUE)
slope1 <- as.vector( (slope-mslope) / sdslope )

# Standardize the estimated anuran richness
ric <- as.vector (habitat$ric.mean)
mric <- mean(ric, na.rm=TRUE)
sdric <- sd(ric, na.rm=TRUE)
ric1 <- as.vector( (ric-mric) / sdric )

# Standardize
edge <- as.vector (habitat$edge_forest)
medge <- mean(edge, na.rm=TRUE)
sdedge <- sd(edge, na.rm=TRUE)
edge1 <- as.vector( (edge - medge) / sdedge )

#----- Are any of the covariates correlated?
## Functions for pairs correlation plot
panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "skyblue1", ...)
}


panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor = 1, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use = "complete.obs")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  text(0.5, 0.5, paste("R^2", " = ", txt, sep = ""), cex = 2.5)
}

panel.lower <- function(x, y, ...){
  points(x, y, col = rgb(0, 0, 0, 0.5), pch = 19, cex=1.3)
  abline(lm(y~x), lwd = 1.5, col = "steelblue")
}


newdat <- data.frame(species_rich = ric1, slope = slope1, 
                     Stream_density = hydro1, forest = forest1, edge = edge1)

# Pairs correlation plot 
tiff(
  "fig02.tiff",
  width     = 5,
  height    = 5,
  units     = "in",
  res       = 900,
  pointsize = 5.5
)

pairs(newdat, 
      upper.panel = panel.cor, 
      diag.panel = panel.hist,
      lower.panel = panel.lower)

dev.off()

# Write the model code to a text file 
sink("Frog_Bd_BAF5.txt")
cat("
    model{
    
    ####################
    # Define priors
    ####################
    
    # Ocupancy covariate paramaters (psi)
    a0 ~ dnorm (0, 0.368)    # Intercept  
    a1 ~ dnorm (0, 0.368)    # Forest cover
    a2 ~ dnorm (0, 0.368)    # Hydrography density
    a3 ~ dnorm (0, 0.368)    # Anuran richness
    
    # Probability of infection covariate paramaters (p)
    beta ~ dnorm(0, 0.368) # Intercept for prevalence model
    b0 ~ dnorm (0,0.368)   # direct developers
    b1 ~ dnorm (0,0.368)   # species breeding in aquatic habitats but occupying the arboreal stratum
    b2 ~ dnorm (0,0.368)   # species breeding in aquatic habitats and occupying the margins of streams and other bodies of water
    b3 ~ dnorm (0,0.368)   # Forest cover
    b4 ~ dnorm (0,0.368)   # Date - linear
    b5 ~ dnorm (0,0.368)   # Date - quadratic
    
    # Mean infection intensity (x) of each site
    n0 ~ dnorm(0, 0.01)  # Intercept
    n1 ~ dnorm(0, 0.01)  # Forest cover 
    n2 ~ dnorm(0, 0.01)  # Stream density
    
    # Variation in infection intensity across sites
    tau <- 1/ (sd*sd)
    sd ~ dunif(0.001, 5)
    
    # Variation in infection intensity at each site
    tau.var <- 1/ (sd.var*sd.var)
    sd.var ~ dunif(0.001, 5)
    
    # Detection probability of Bd (informative priors by Miller et al. 2012)
    alpha_p ~ dunif(0.25, 1.32)
    beta_p ~ dunif(0.14, 0.51)
    
    # Measurement error of infection intensity
    tau.p <- 1/ (sd.p * sd.p)
    sd.p ~ dunif(1.47, 1.97)
    
    ##################################
    # Ecological model likelihood
    ##################################
    
    # Bd occupancy probability across sites  
    
    for(j in 1:nsites){  # For each site
    
    logit(psi[j]) <- a0 + a1 * forest1[j] + a2 * hydro1[j] + a3 * ric1[j]
    
    z[j] ~ dbern(psi[j])  # Bd occupancy probability across sites  
    
    # Bd Infection intensity across sites
    
    mu[j] <-  n0 + n1 * forest1[j] + n2 * ric1[j]
    
    mu_z[j] <- mu[j] * z[j]
    
    x[j] ~ dnorm(mu[j], tau) # Mean Bd infection intensity at site-level
    
    for (k in 1:n.swabbed[j]){   # Different number of frogs sampled (swabbed) at each site j
    
    # Probability of Bd infection at each site
    
    # Describes the relationship of prob. Bd infection and aquatic index, forest cover, and date    
    logit(p[j, k]) <- beta + b0 * bd.AI0[j, k] + b1 * bd.AI1[j, k] + b2 * bd.AI2[j, k] + 
    b3 * forest[j] + b4 * date1[j, k] + b5 * date2[j, k]
    
    mu_inf[j, k] <- p[j, k] * z[j] # There has to be Bd at a site to calculate prob. of infection
    
    detections[j, k] ~ dbern(mu_inf[j, k])  # Corrected Bd detections
    
    ##################################
    # Sampling model likelihood
    ##################################   
    
    observed[j,k] ~ dbern(p.eff[j,k])  # Observed Bd detections
    
    p.eff[j,k] <- path_p[j,k] * detections[j,k]
    
    logit(path_p[j,k]) <- alpha_p + beta_p * N[j,k]
    # At the individual scale- correct for imperfect pathogen detection by Miller et al. 2012
    # Bd detection probability depends on infection intensity N[j,k]
    
    # Variation in infection intensity at each site
    N[j,k] ~ dnorm(mu_x[j,k], tau.var)   # Corrected Bd infection intensity
    mu_x[j,k] <- x[j] * z[j]
    
    w[j,k] ~ dnorm(mu_N[j,k], tau.p)  # Observed Bd infection intensity
    mu_N[j,k] <- N[j,k] * observed[j,k]
    # Correct for measurement error in estimates
    # If we had the data, w would have a 3rd dimension with replicate qPCR or swabs, then this would be used to calculate the true value of Bd infection intensity per individual at each site (N[j,k])
    
    } #k
    
    } #j
    
    ##################################
    # Derived quantities
    ################################## 
    # Absence/presence
    mean.occ  <- sum(z[])                           # number of occupied sites by Bd
    occ.rate <- exp(a0) / (1 + exp(a0))                 # average occupancy probability
    prev.AI0 <- exp(beta + b0) / (1 + exp(beta + b0))   # prevalence of terrestrial breeding species
    prev.AI1 <- exp(beta + b1) / (1 + exp(beta + b1))   # prevalence of species breeding in aquatic habitats but occupying the arboreal
    prev.AI2 <- exp(beta + b2) / (1 + exp(beta + b2))   # prevalence of species breeding in aquatic habitats and occupying the margins of streams and other bodies of water
    
    # Infection intensity
    for(j in 1:nsites){
    Ncols0[j] <- mean(N[j,n.swabbed[j]] * bd.AI0[j,n.swabbed[j]])
    Ncols1[j] <- mean(N[j,n.swabbed[j]] * bd.AI1[j,n.swabbed[j]])
    Ncols2[j] <- mean(N[j,n.swabbed[j]] * bd.AI2[j,n.swabbed[j]])
    }
    int.AI0 <- mean(Ncols0[])
    int.AI1 <- mean(Ncols1[])
    int.AI2 <- mean(Ncols2[])
    # Finish writing the text file into a document called model
    }
    ", fill = TRUE)
sink()

#Create the necessary arguments to run the jags() command 
#Load all the data
sp.data = list (nsites = nsites,        # Number of sites
                n.swabbed = n.swabbed,  # Number of swabbed frogs per site 
                forest1 = forest1,      # forest cover
                hydro1 = hydro1,   # hydrology
                ric1 = ric1,      # species richness
                bd.AI0 = bd.AI0,  # AI = 0
                bd.AI1 = bd.AI1,  # AI = 1
                bd.AI2 = bd.AI2,  # AI = 2
                date1 = date1,    # date
                date2 = date2,    # date squared
                observed = bd.pres, # site x frog matrix with Bd presence/absence
                w = log(bd.intensity) # site x frog matrix with Bd infection intensity
)

#Specify the parameters to be monitored
sp.params = c("mean.occ", "occ.rate",
              "mean.x", "mean.N", "mean.w",
              "int.AI0", "int.AI1", "int.AI2",
              "a0", "a1", "a2", "a3",
              "beta", "b0", "b1", "b2", "b3", "b4", "b5",
              "prev.AI0", "prev.AI1", "prev.AI2",
              "sd.p", "sd", "sd.var",
              "n0", "n1", "n2"
)

# Specify the initial values
zst <-apply(bd.pres, 1, max, na.rm=T)
zst [zst == "-Inf"]<- 0
names(zst) <- NULL

xst <- apply((bd.intensity), 1, max, na.rm = TRUE)
xst[xst == "-Inf"] <- 1

dets <- bd.pres
dets[is.na(dets) == TRUE] <- 1

Nst <- as.matrix(bd.intensity)
Nst[is.na(Nst) == TRUE] <- 1

sp.inits = function (){ list(
  z = zst,
  x = log(xst),
  N = as.matrix(log(bd.intensity)),
  detections = as.matrix(bd.pres),
  
  alpha_p = runif(1, 0.25, 1.32),
  beta_p = runif(1, 0.14, 0.51),
  
  sd.p = runif(1, 1.47, 1.97),
  sd.var = runif(1, 0.001, 5),
  sd = runif(1, 0.001, 5),
  
  n0 = runif(1, 0, 5),
  n1 = runif(1, 0, 5),
  n2 = runif(1, 0, 5),
  
  a0 = runif(1, -3, 3),
  a1 = runif(1, -3, 3),
  a2 = runif(1, -3, 3),
  a3 = runif(1, -3, 3),
  
  beta = runif(1, -3, 3),
  b0  = runif(1, -3, 3),
  b1  = runif(1, -3, 3),
  b2  = runif(1, -3, 3),
  b3  = runif(1, -3, 3),
  b4  = runif(1, -3, 3)
)
}

# MCMC test settings
ni <- 55000; nt <- 50; nb <- 5000; nc <- 3; na = 50000

#Load the R2Jags library
library(jagsUI)

#Run the model and call the results out (around 6.7 minutes)
out5 <- jags(data = sp.data, inits = sp.inits, parameters.to.save = sp.params, 
             model = "Frog_Bd_BAF5.txt", n.chains = nc, n.thin = nt, n.iter = ni, 
             n.burnin = nb, parallel = TRUE, store.data = TRUE, n.adapt = na)

# Look at output
out1

# Look at traceplots of output
plot(out1)

traceplot(out1)
