BINDIR=`dirname $(readlink -f "$0")`
OUTDIR=`pwd`
readfilter=$1
minq=$2
newcontigsfile=$OUTDIR/'chm13_'$readfilter'_'$minq'_newscaffolds.fa'
outfile=$OUTDIR/'chm13_'$readfilter'_'$minq'_allscaffolds.fa'
contigsfn='/scratch/groups/mschatz1/mkirsche/t2t-chm13.20200602.fasta'
echo 'Compiling'
javac $BINDIR/src/*.java
echo 'Running scaffolding'
java -cp src Main aln_fn='/scratch/groups/mschatz1/mkirsche/minimap_t2t-chm13.20200602_'$readfilter.paf fasta_fn=$contigsfn read_fn='/scratch/groups/mschatz1/mkirsche/'rel5_$readfilter.fastq read_map_file=$OUTDIR/'chm13_'$readfilter'_'$minq'_usefulreadsmap.paf' contig_map_file=$OUTDIR/'chm13_'$readfilter'_'$minq'_usefulcontigsmap.paf' chm13_$readfilter_$minq out_file=$newcontigsfile minq=$minq full_out_gfa_fn=$OUTDIR/'scaffolds_'$readfilter'_'$minq'.full.gfa' joins_out_gfa_fn=$OUTDIR/'scaffolds_'$readfilter'_'$minq'.joins.gfa' min_weight_supp=150 min_weight=10
echo 'Stitching in new contigs'
java -cp $BINDIR/src StitchFasta $contigsfn $newcontigsfile $outfile
echo 'Done'
echo ''
