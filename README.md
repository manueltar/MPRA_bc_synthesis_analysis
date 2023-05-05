# MPRA_bc_synthesis_Sample_alignment_and_counts MPRA analysis


$ bash MPRA_pipeline.sh \<path to leave output files\> \<path to the MPRA_programmed_pipeline\> \<path to the FINAL output files from sample_processing\> \<mem\> \<processors\> \<bsub queue\>


$ bash ~/Scripts/Wraper_scripts/298_MPRA_P_redo_pipeline_v3.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/GLOBAL_ANALYSIS_Feb_2022_FC01/ 16000 4 normal


$ bash MPRA_pipeline.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST_ANALYSIS/ /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/ /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST/ 16000 4 normal

## FILE DEPENDENCIES

untar the ALL_db.tsv.tar.gz in the Dependencies folder

## SOFTWARE DEPENDENCIES

-R-4.1.0, you have to set the path to Rscript at the top of the bash script: Rscript=/software/R-4.1.0/bin/Rscript
- Rpackages. Important! change the default lib.loc in every script as it direct to my folder in the Sanger farm. See per script but all told they are:

optparse
MPRAnalyze
reshape2
labeling
ggrepel
ggeasy
farver
cowplot
zoo
withr
viridisLite
tzdb
TxDb.Hsapiens.UCSC.hg19.knownGene
ttutils
tidyverse
sysfonts
svglite
Sushi
SummarizedExperiment
showtextdb
showtext
S4Vectors
R.utils
rtracklayer
rstudioapi
R.oo
R.methodsS3
plyr
org.Hs.eg.db
OrganismDbi
memoise
matrixStats
MatrixGenerics
markdown
liftOverlib.loc
lemon
labeling
jsonlite
IRanges
Homo.sapiens
gwascat
gtools
gridExtra
gridBase
grid
GO.db
glue
ggtext
ggridges
ggplot2
ggforce
ggeasy
GenomicRanges
GenomicFeatures
GenomeInfoDb
farver
extrafontdb
extrafont
dplyr
digest
digest
data.table
curl
crayon
cowplot
cli
broom
biomaRt
BiocGenerics
Biobase
backports
AnnotationDbi
 
