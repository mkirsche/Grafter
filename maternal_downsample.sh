WORKINGDIR=`pwd`
OUTDIR=$WORKINGDIR'/'$1
rm -r $OUTDIR
mkdir $OUTDIR
BINDIR=`dirname $(readlink -f "$0")`
len=100
prop=$2

readfn='/scratch/groups/mschatz1/mkirsche/ultralong/rel2_'$len'kplus.fastq'
contigsfn='/scratch/groups/mschatz1/mkirsche/ultralong/ccs/maternal_and_unknown.contigs.mmpoa.fa'
outfn=$OUTDIR'/maternal_and_unknown.contigs.mmpoa.scaffoldsgraph.fa'

echo 'Sampling reads'
filteredreadfn=$OUTDIR'/filteredreads.fastq'
/scratch/groups/mschatz1/mkirsche/github/seqtk/seqtk sample $readfn $prop > $filteredreadfn

echo 'Aligning with minimap2'
paffn=$OUTDIR'/alignments.paf'
$minimappath -t 32 -k 19 -w 19 $contigsfn $filteredreadfn  > $paffn

$BINDIR/run2.sh $contigsfn $OUTDIR'/filteredreads.fastq'  $paffn $outfn 1 $OUTDIR 2>&1 | tee $OUTDIR'/'maternal.log

java -cp ~/hashing/CCS/assembly_eval/ AssemblyStats $outfn | tee $OUTDIR'/'maternal_afterstatsgraph.txt;

$BINDIR/delta.sh $outfn $OUTDIR'/'after_maternalgraph;

