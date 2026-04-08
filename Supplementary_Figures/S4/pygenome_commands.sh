# Figure S4 : RNA seq tracks for AICDA locus for all cell types.
# RNA-seq +ve strand tracks for AICDA locus for all cell types.
pyGenomeTracks --tracks rna_seq_Unique_plusRaw.ini --region chr12:8527000-9650000 --dpi 300 --outFileName RNA_plus_AICDA_chr12-8650000-9000000.png --fontSize 20 --height 50 --trackLabelFraction 0 --plotWidth 50

# RNA-seq -ve strand tracks for AICDA locus for all cell types.
pyGenomeTracks  --tracks rna_seq_Unique_minusRaw.ini --region chr12:8527000-9650000 --dpi 300 --outFileName RNA_negative_AICDA_chr12-8650000-9000000.png --fontSize 20 --height 50 --trackLabelFraction 0 --plotWidth 50

# Segementation tracks for AICDA locus for all cell types.
pyGenomeTracks  --tracks episegmixmeth_relabelled.ini --region chr12:8527000-9650000 --dpi 300 --outFileName ESMM_AICDA_chr12-8650000-9000000.png --fontSize 20 --height 50 --trackLabelFraction 0 --plotWidth 50