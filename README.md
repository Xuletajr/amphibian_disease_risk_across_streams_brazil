# Assessing amphibian disease risk across tropical streams while accounting for imperfect pathogen detection

### José W. Ribeiro Jr, Tadeu Siqueira, Graziella V. DiRenzo, Carolina Lambertini, Mariana L. Lyra, Luís F. Toledo, Célio F. B. Haddad, C. Guilherme Becker

### Oecologia

### Code DOI: https://doi.org/10.1007/s00442-020-04646-4

### Please contact the first author for questions about the code or data: José W. Ribeiro Jr (jwribeirojunior@gmail.com)
__________________________________________________________________________________________________________________________________________
## Abstract:
Ecologists studying emerging wildlife diseases need to confront the realism of imperfect pathogen detection across heterogeneous habitats to aid in conservation decisions. For example, spatial risk assessments of amphibian disease caused
by Batrachochytrium dendrobatidis (Bd) has largely ignored imperfect pathogen detection across sampling sites. Because
changes in pathogenicity and host susceptibility could trigger recurrent population declines, it is imperative to understand how
pathogen prevalence and occupancy vary across environmental gradients. Here, we assessed how Bd occurrence, prevalence,
and infection intensity in a diverse Neotropical landscape vary across streams in relation to abiotic and biotic predictors
using a hierarchical Bayesian model that accounts for imperfect Bd detection caused by qPCR error. Our model indicated
that the number of streams harboring Bd-infected frogs is higher than observed, with Bd likely being present at ~ 43% more
streams than it was detected. We found that terrestrial-breeders captured along streams had higher Bd prevalence, but lower
infection intensity, than aquatic-breeding species. We found a positive relationship between Bd occupancy probability and
stream density, and a negative relationship between Bd occupancy probability and amphibian local richness. Forest cover
was a weak predictor of Bd occurrence and infection intensity. Finally, we provide estimates for the minimum number of
amphibian captures needed to determine the presence of Bd at a given site where Bd occurs, thus, providing guidence for costeffective disease risk monitoring programs.

## Code
__*bd_occupancy_prevalence_Brazil_Atlantic_Forest.R*__ - Cotains R code to import and reshape raw data, and run hierarchical Bayesian model (Miller et al. 2012, DiRenzo et al. 2018) to estimate *Bd* occurrence, prevalence, and infection intensity (via JAGS through R).

__*output_frogBd_BAF.R*__ - Contains R code to estimate posterior summary statistics and to create figures. 

## Data
The rows are the 49 sampled sites (streams). The number of individuals captured per stream ranged from 1 to 10 amphibians (columns). We followed the protocol described by Boyle et al. (2004) and extracted *Bd* DNA from each swab using PrepMan Ultra® (Applied Biosystems).  We quantified *Bd* infection loads (zoospore genomic equivalents: ZGE) in each swab using qPCR analysis with TaqMan assays (Applied Biosystems). Each plate contained *Bd* standards of 0.1, 1, 10, 100, and 1000 ZGE from strain CLFT 159. The samples that had *Bd* load equal to or greater than one were treated as *Bd*-positive (Kriger et al. 2006, 2007).

### *Bd* presence
__*bd_pres.txt*__ - contains the data of *Bd* occurrence in Brazilian Atlantic Forest streams. Presence = 1; Absence = 0; NA's = not sampled.

### *Bd* load
__*bd_load.txt*__ - contains the data of *Bd* infection intensity (load) in Brazilian Atlantic Forest streams. Load = > 0 (ZGE); Absence = 0; NA's = not sampled.

### Aquatic index
Aquatic index is a measure of amphibian water dependency and has been used to relate the probability a species will be infected by Bd in amphibian communities (Becker et al. 2014; Mesquita et al. 2017). We classified species into three distinct groups:

__*bd_AI0.txt*__ - terrestrial species with terrestrial eggs (AI-0). 

__*bd_AI1.txt*__ - arboreal species with an aquatic larval phase (AI-1).

__*bd_AI2.txt*__ - terrestrial species with an aquatic larval phase (AI-2).

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

## References
- Becker CG, Rodriguez D, Toledo LF et al (2014) Partitioning the net effect of host diversity on an emerging amphibian pathogen. Proc
R Soc B 281:20141796. https://doi.org/10.1098/rspb.2014.1796
- Boyle DG, Boyle DB, Olsen V et al (2004) Rapid quantitative detection of chytridiomycosis (Batrachochytrium dendrobatidis) in amphibian samples using real-time Taqman PCR assay. Dis Aquat Organ 60:141–148. https://doi.org/10.3354/dao060141
- DiRenzo GV, Grant EHC, Longo AV, et al (2018) Imperfect pathogen detection from non-invasive skin swabs biases disease inference. Methods Ecol Evol 9:380–389. doi: https://doi.org/10.1111/2041-210X.12868
- Kriger KM, Ashton KJ, Hines HB, Hero J (2007) On the biological relevance of a single Batrachochytrium dendrobatidis zoospore:
a reply to Smith (2007). Dis Aquat Organ 73:257–260. https://doi.org/10.3354/dao073257
- Kriger KM, Hines HB, Hyatt AD et al (2006) Techniques for detecting chytridiomycosis in wild frogs: comparing histology with realtime Taqman PCR. Dis Aquat Organ 71:141–148. https://doi.org/10.3354/dao071141
- Mesquita AFC, Lambertini C, Lyra M et al (2017) Low resistance to chytridiomycosis in direct-developing amphibians. Sci Rep
7:16605. https://doi.org/10.1038/s41598-017-16425-y
- Miller DAW, Talley BL, Lips KR, Campbell Grant EH (2012) Estimating patterns and drivers of infection prevalence and intensity when detection is imperfect and sampling error occurs. Methods Ecol Evol 3:850–859. doi: https://doi.org/10.1111/j.2041-210X.2012.00216.x
- Ribeiro JW, Siqueira T, Brejão GL, Zipkin EF (2018) Effects of agriculture and topography on tropical amphibian species and communities. Ecol Appl 28:1554–1564. doi: https://doi.org/10.1002/eap.1741
