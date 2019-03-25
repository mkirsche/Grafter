WORKINGDIR=`pwd`
OUTDIR=$WORKINGDIR'/'$1
rm -r $OUTDIR
mkdir $OUTDIR
BINDIR=`dirname $(readlink -f "$0")`
len=100
prop=$2

echo 'Sampling reads'
readfn='/scratch/groups/mschatz1/mkirsche/ultralong/rel2_'$len'kplus.fastq'

/scratch/groups/mschatz1/mkirsche/github/seqtk/seqtk sample $readfn $prop > $OUTDIR'/filteredreads.fastq'

$BINDIR/run2.sh /scratch/groups/mschatz1/mkirsche/ultralong/ccs/maternal_and_unknown.contigs.mmpoa.fa $OUTDIR'/filteredreads.fastq'  '/scratch/groups/mschatz1/mkirsche/ultralong/ccs/rel2_'$len'kplus_ccs_mat.paf' $OUTDIR/maternal_and_unknown.contigs.mmpoa.scaffoldsgraph.fa 1 $OUTDIR 2>&1 | tee $OUTDIR'/'maternal.log

java -cp ~/hashing/CCS/assembly_eval/ AssemblyStats $OUTDIR'/'maternal_and_unknown.contigs.mmpoa.scaffoldsgraph.fa | tee $OUTDIR'/'maternal_afterstatsgraph.txt;

$BINDIR/delta.sh $OUTDIR'/'maternal_and_unknown.contigs.mmpoa.scaffoldsgraph.fa $OUTDIR'/'after_maternalgraph;

