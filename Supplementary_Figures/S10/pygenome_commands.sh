# Figure S10 : RNA seq tracks for AICDA locus for all cell types.
# NBC
pyGenomeTracks  --tracks NBC_long_region_HiC.ini --region chr2:55683901-61110893 --dpi 300 --outFileName NBC.png 

# MBC
pyGenomeTracks  --tracks MBC_long_region_HiC.ini --region chr2:55683901-61110893 --dpi 300 --outFileName MBC.png

# GCBC
pyGenomeTracks  --tracks GCBC_long_region_HiC.ini --region chr2:55683901-61110893 --dpi 300 --outFileName GCBC.png

# PC
pyGenomeTracks  --tracks PC_long_region_HiC.ini --region chr2:55683901-61110893 --dpi 300 --outFileName PC.png
