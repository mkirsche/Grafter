BINDIR=`dirname $(readlink -f "$0")`
OUTDIR=`pwd`
javac $BINDIR/src/*.java
java -cp $BINDIR/src PrepareIgv contig_file=/home/mkirsche/physalis/purged.fa agp_file=/home/mkirsche/eclipse-workspace/Grafter/test_ragtag.agp scaffold_file=/home/mkirsche/physalis/scaffolds_ragtag.fa ill_tdf_file=/home/mkirsche/physalis/M82_grafter_self.sr.tdf ont_tdf_file=/home/mkirsche/physalis/M82_grafter_self.ont.tdf

xvfb-run --auto-servernum /home/mkirsche/Downloads/IGV_Linux_2.8.2/igv.sh -b $BINDIR/igv/igv.bat

