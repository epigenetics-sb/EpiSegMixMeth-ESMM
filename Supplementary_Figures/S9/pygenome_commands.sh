# Figure S9 : ESMM vs ChromHMM tracks for AICDA locus in all cell types.
# ESMM 
pyGenomeTracks  --tracks episegmixmeth_relabelled.ini --region chr12:8650000-9000000 --dpi 300 --outFileName ESMM_AICDA_chr12-8650000-9000000.png 

# ChromHMM
pyGenomeTracks  --tracks chromhmm_relabelled.ini --region chr12:8650000-9000000 --dpi 300 --outFileName ChromHMM_AICDA_chr12-8650000-9000000.png 

