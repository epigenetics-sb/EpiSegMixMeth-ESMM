cd ./project/Analysis/epilogos


# Figure 1 
pyGenomeTracks  --tracks ../pygenome/singals/H3K4me3.signals.bigwig.ini --region chr12:8550000-9000000 -o ../pygenome/singals/H3K4me3.signals.bigwig.ini.png --dpi 300 --height 40 --trackLabelFraction 0 --plotWidth 40

pyGenomeTracks  --tracks ../pygenome/singals/H3K4me1.signals.bigwig.ini --region chr12:8550000-9000000 -o ../pygenome/singals/H3K4me1.signals.bigwig.ini.png --dpi 300 --height 40 --trackLabelFraction 0 --plotWidth 40

pyGenomeTracks  --tracks ../pygenome/singals/H3K27ac.signals.bigwig.ini --region chr12:8550000-9000000 -o ../pygenome/singals/H3K27ac.signals.bigwig.ini.png --dpi 300 --height 40 --trackLabelFraction 0 --plotWidth 40

pyGenomeTracks  --tracks ../pygenome/singals/H3K36me3.signals.bigwig.ini --region chr12:8550000-9000000 -o ../pygenome/singals/H3K36me3.signals.bigwig.ini.png --dpi 300 --height 40 --trackLabelFraction 0 --plotWidth 40

pyGenomeTracks  --tracks ../pygenome/singals/H3K27me3.signals.bigwig.ini --region chr12:8550000-9000000 -o ../pygenome/singals/H3K27me3.signals.bigwig.ini.png --dpi 300 --height 40 --trackLabelFraction 0 --plotWidth 40

pyGenomeTracks  --tracks ../pygenome/singals/H3K9me3.signals.bigwig.ini --region chr12:8550000-9000000 -o ../pygenome/singals/H3K9me3.signals.bigwig.ini.png --dpi 300 --height 40 --trackLabelFraction 0 --plotWidth 40

pyGenomeTracks  --tracks ../pygenome/singals/WGBS.signals.bigwig.ini --region chr12:8550000-9650000 -o ../pygenome/singals/WGBS.signals.bigwig.ini.png --dpi 300 --height 40 --trackLabelFraction 0 --plotWidth 40

pyGenomeTracks  --tracks ../pygenome/bedfiles/episegmixmeth.ini --region chr12:8550000-9650000 -o ../pygenome/bedfiles/episegmixmeth.ini.png --dpi 300 --height 40 --trackLabelFraction 0 --plotWidth 40

pyGenomeTracks  --tracks ../pygenome/bedfiles/episegmixmeth_relabelled.ini --region chr12:8550000-9650000 -o ../pygenome/bedfiles/episegmixmeth_relabelled.ini.png --dpi 300 --height 40 --trackLabelFraction 0 --plotWidth 40
