# Assessing amphibian disease risk across tropical streams while accounting for imperfect pathogen detection

### José W. Ribeiro Jr, Tadeu Siqueira, Graziella V. DiRenzo, Carolina Lambertini, Mariana L. Lyra, Luís F. Toledo, Célio F. B. Haddad, C. Guilherme Becker

### Oecologia (submitted)

### Code DOI:

### Please contact the first author for questions about the code or data: José W. Ribeiro Jr (jwribeirojunior@gmail.com)
__________________________________________________________________________________________________________________________________________
## Abstract:
Ecologists studying emerging wildlife diseases need to confront the realism of imperfect pathogen detection across heterogeneous habitats to enhance conservation actions. For example, spatial risk assessments of amphibian disease caused by *Batrachochytrium dendrobatidis* (*Bd*) has largely ignored imperfect detection across sampling sites. Because changes in pathogenicity and host susceptibility could trigger recurrent population declines, it is imperative to understand how pathogen prevalence and occupancy vary across environmental gradients. We assessed how *Bd* occurrence, prevalence, and infection intensity in a diverse Neotropical landscape vary across streams in relation to abiotic and biotic predictors using a hierarchical Bayesian model that accounts for imperfect *Bd* detection caused by qPCR error. Our model indicated that the number of streams harboring *Bd*-infected frogs is higher than observed, with *Bd* likely being present at 55% more streams than it was detected. We found that terrestrial-breeders captured along streams had higher infection prevalence, but lower infection intensity, than aquatic-breeding species. We found a positive relationship between Bd occupancy and stream density, and a negative relationship between *Bd* occupancy and amphibian richness. Forest cover was a weak predictor of *Bd* prevalence and intensity. Lastly, we provide estimates for the minimum sampling effort needed to detect *Bd* in a given sampling site where *Bd* occurs, guiding cost-effective disease risk monitoring programs. Our study underscores that hierarchical Bayesian models that account for pathogen uncertainty in pathogen detection grant more precise estimations of occurrence, prevalence, and infection intensity, and evaluating the role of abiotic and biotic variables on pathogen spatial distributions.

## Data

### *Bd* presence
__*bd_pres.txt*__ -

### *Bd* load
__*bd_load.txt*__ -

### Aquatic index
__*bd_AI0.txt*__ -

__*bd_AI0.txt*__ - 

__*bd_AI0.txt*__ - 

### Survey date
__*date_bd.txt*__ - 

### Habitat covariates
habitat_covariates.txt - contains habitat covariate information for each site. 
1. __"stream"__ - contains sampling site id
2. __"long"__ - is the geographic longitude coordinate as UTM
3. __"lat"__ - is the geographic latitude coordinate as UTM
4. __"altitude"__ - altitude of each site
5. __"forest"__ - is the proportion of natural forest cover area within a buffer of 200 m radius
6. __"stream_length"__ - is the stream length network (m) within a buffer of 200 m radius
7. __“slope200_sd”__ - is the standard deviation slope within a buffer of 200 m radius, it was derived from the Digital Elevation Model raster image (30-m resolution) from Shuttle Radar Topography Mission (SRTM).
8. __“ric.mean”__ - estimates of amphibian species richness published by Ribeiro et al. (2018) − a study carried out in the same focal streams at the same time as the present study − as a fine-scale measure of local amphibian diversity. 
9. __“edge_forest”__ - is the forest edge (m) within a buffer of 200 m radius

## Codes
__*bd_occupancy_prevalence_Brazil_Atlantic_Forest.R*__ - R code to run hierarchical Bayesian model to estimate *Bd* occurrence, prevalence, and infection intensity. Contains code to import and reshape the data, summary statistics, and run the model file in JAGS.

__*output_frogBd_17_may_2019.R*__ - 
