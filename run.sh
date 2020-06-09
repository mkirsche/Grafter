contigsfn=$1
readsfn=$2
paffn=$3
outfile=$4
minq=$5
minimappath=/scratch/groups/mschatz1/mkirsche/github/minimap2/minimap2

if [ ! -f $paffn ]; then
    echo "PAF file not found - generating it"
    $minimappath -t 8 $contigsfn $readsfn  > $paffn
fi

BINDIR=`dirname $(readlink -f "$0")`

echo 'Compiling'
javac $BINDIR/scaffolding/*.java

usefulpaf=$paffn'_useful_'$minq'.paf'
echo 'Finding useful alignments'
java -cp $BINDIR scaffolding.FindUsefulScaffoldingAlignments aln_fn=$paffn out_file=$usefulpaf minq=$minq
echo 'Useful alignments output to '$usefulpaf

readmap=$readsfn'_usefulmap_'$minq'.paf'
contigmap=$contigsfn'_usefulmap_'$minq'.paf'
newcontigs=$contigsfn'_newcontigs_'$minq'.paf'
echo 'Scaffolding'
java -cp $BINDIR scaffolding.Scaffold aln_fn=$usefulpaf fasta_fn=$contigsfn read_fn=$readsfn read_map_file=$readmap contig_map_file=$contigmap out_file=$newcontigs
echo 'Scaffolds output to '$newcontigs

echo 'Integrating scaffolds into assembly'
java -cp $BINDIR scaffolding.StitchFasta $contigsfn $newcontigs $outfile
echo 'Final assembly output to '$outfile
