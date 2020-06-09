BINDIR=`dirname $(readlink -f "$0")`
OUTDIR=`pwd`
readfilter=$1
minq=$2
newcontigsfile=$OUTDIR/'chm13_'$readfilter'_'$minq'_newscaffolds.fa'
outfile=$OUTDIR/'chm13_'$readfilter'_'$minq'_newscaffolds.fa'
contigsfn='/work-zfs/mschatz1/CHM13/drafts/20200602/simplified.nodes.fasta'
echo 'Compiling'
javac $BINDIR/src/*.java
echo 'Running scaffolding'
java -cp src IncludeContained aln_fn=/work-zfs/mschatz1/CHM13/minimap_ul_0602_$readfilter.paf fasta_fn=$contigsfn read_fn=/work-zfs/mschatz1/CHM13/raw/UL_ONT/rel5_$readfilter.fastq read_map_file=$OUTDIR/'chm13_'$readfilter'_'$minq'_usefulreadsmap.paf' contig_map_file=$OUTDIR/'chm13_'$readfilter'_'$minq'_usefulcontigsmap.paf' chm13_$readfilter_$minq out_file=$newcontigsfile minq=$minq
echo 'Stitching in new contigs'
java -cp $BINDIR/src StitchFasta $contigsfn $newcontigsfile $outfile
echo 'Done'
echo ''
