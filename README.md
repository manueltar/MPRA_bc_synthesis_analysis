# MPRA_bc_synthesis_Sample_alignment_and_counts MPRA analysis


$ bash MPRA_analysis.sh \<folder for output files\> \<folder for the deduplicated count files, see Sample processing\> \<Main folder for the analysis\> \<memory\> \<processors\> \<queue\>

$ bash MPRA_analysis.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline_analysis/Output_folder/ /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/Output_folder/ /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline_analysis/ 8000 2 normal


## FILE DEPENDENCIES

untar the ALL_db.tsv.tar.gz in the Dependencies folder

## SOFTWARE DEPENDENCIES

-R-4.1.0, you have to set the path to Rscript at the top of the bash script: Rscript=/software/R-4.1.0/bin/Rscript
- Rpackages. Important! change the default lib.loc in every script as it direct to my folder in the Sanger farm. See per script but all told they are:

"Biobase","BiocGenerics","BiocManager","cowplot","curl","digest","digest","extrafont","extrafontdb","farver","GENESIS","ggeasy","ggforce","ggrepel","ggridges","ggtext","glue","gtools","GWASTools","jsonlite","labeling","labeling","markdown","MatrixGenerics","matrixStats","memoise","MPRAnalyze","RColorBrewer","reshape2","R.methodsS3","R.oo","R.utils","showtext","showtextdb","SummarizedExperiment","svglite","sysfonts","ttutils","viridisLite","AnnotationDbi","backports","biomaRt","broom","cli","cowplot","crayon","data.table","dplyr","GenomeInfoDb","GenomicFeatures","GenomicRanges","ggnewscale","ggplot2","GO.db","grid","gridBase","gridExtra","gwascat","Homo.sapiens","IRanges","labeling","lemon"