BINDIR=`dirname $(readlink -f "$0")`
OUTDIR=`pwd`
readfilter='200k'
chr=$1
minq=$2
newcontigsfile=$OUTDIR/'chm13_'$chr'_'$readfilter'_'$minq'_newscaffolds.fa'
outfile=$OUTDIR/'chm13_'$chr'_'$readfilter'_'$minq'_allscaffolds.fa'
contigsfn='/scratch/groups/mschatz1/mkirsche/t2tglobus/chr'$chr'_pieces.fasta'
alnfn='/scratch/groups/mschatz1/mkirsche/t2tglobus/minimap/chr'$chr'_'$readfilter'.paf'
readfn='/scratch/groups/mschatz1/mkirsche/rel5_'$readfilter'.fastq'
echo 'Compiling'
javac $BINDIR/src/*.java
echo 'Running scaffolding'
java -cp src Main aln_fn=$alnfn fasta_fn=$contigsfn read_fn=$readfn read_map_file=$OUTDIR/'chm13_'$chr'_'$readfilter'_'$minq'_usefulreadsmap.paf' contig_map_file=$OUTDIR/'chm13_'$chr'_'$readfilter'_'$minq'_usefulcontigsmap.paf' out_file=$newcontigsfile minq=$minq full_out_gfa_fn=$OUTDIR/'scaffolds_'$chr'_'$readfilter'_'$minq'.full.gfa' joins_out_gfa_fn=$OUTDIR/'scaffolds_'$chr'_'$readfilter'_'$minq'.joins.gfa' min_weight_supp=150 min_weight=10
echo 'Stitching in new contigs'
java -cp $BINDIR/src StitchFasta $contigsfn $newcontigsfile $outfile
echo 'Done'
echo ''
