# Assessing amphibian disease risk across tropical streams while accounting for imperfect pathogen detection

### José W. Ribeiro Jr, Tadeu Siqueira, Graziella V. DiRenzo, Carolina Lambertini, Mariana L. Lyra, Luís F. Toledo, Célio F. B. Haddad, C. Guilherme Becker

### Oecologia (submitted)

### Code DOI:

### Please contact the first author for questions about the code or data: José W. Ribeiro Jr (jwribeirojunior@gmail.com)
__________________________________________________________________________________________________________________________________________
## Abstract:
Ecologists studying emerging wildlife diseases need to confront the realism of imperfect pathogen detection across heterogeneous habitats to enhance conservation actions. For example, spatial risk assessments of amphibian disease caused by Batrachochytrium dendrobatidis (Bd) has largely ignored imperfect detection across sampling sites. Because changes in pathogenicity and host susceptibility could trigger recurrent population declines, it is imperative to understand how pathogen prevalence and occupancy vary across environmental gradients. We assessed how Bd occurrence, prevalence, and infection intensity in a diverse Neotropical landscape vary across streams in relation to abiotic and biotic predictors using a hierarchical Bayesian model that accounts for imperfect Bd detection caused by qPCR error. Our model indicated that the number of streams harboring Bd-infected frogs is higher than observed, with Bd likely being present at 55% more streams than it was detected. We found that terrestrial-breeders captured along streams had higher infection prevalence, but lower infection intensity, than aquatic-breeding species. We found a positive relationship between Bd occupancy and stream density, and a negative relationship between Bd occupancy and amphibian richness. Forest cover was a weak predictor of Bd prevalence and intensity. Lastly, we provide estimates for the minimum sampling effort needed to detect Bd in a given sampling site where Bd occurs, guiding cost-effective disease risk monitoring programs. Our study underscores that hierarchical Bayesian models that account for pathogen uncertainty in pathogen detection grant more precise estimations of occurrence, prevalence, and infection intensity, and evaluating the role of abiotic and biotic variables on pathogen spatial distributions.

## Data
### Habitat covariates
habitat_covariates.txt - contains habitat covariate information for each site. 
1. "stream" - contains sampling site id
2. "long" - is the geographic longitude coordinate as UTM
3. "lat" - is the geographic latitude coordinate as UTM
4. "altitude" - altitude of each site
5. "forest" - is the proportion of natural forest cover area within a buffer of 200 m radius
6. "stream_length" - is the stream length network (m) within a buffer of 200 m radius
7. “slope200_sd” - is the standard deviation slope within a buffer of 200 m radius, it was derived from the Digital Elevation Model raster image (30-m resolution) from Shuttle Radar Topography Mission (SRTM).
8. “ric.mean” - 
9. “edge_forest” - 
