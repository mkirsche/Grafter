contigsfn=$1
readsfn=$2
paffn=$3

minimappath='~/github/minimap2/minimap2'

if [ ! -f $paffn ]; then
    echo "PAF file not found - generating it"
    $minimappath -t 8 $contigsfn $readsfn  > $paffn
fi

BINDIR=`dirname $(readlink -f "$0")`

echo 'Compiling'
javac $BINDIR/scaffolding/*.java

usefulpaf=$paffn'_useful.paf'
echo 'Finding useful alignments'
java -cp $BINDIR scaffolding.FindUsefulScaffoldingAlignments $paffn $usefulpaf

readmap=$readsfn'_usefulmap.paf'
contigmap=$contigsfn'_usefulmap.paf'
newcontigs=$contigsfn'_newcontigs.paf'
echo 'Scaffolding'
java -cp $BINDIR scaffolding.Scaffold $usefulpaf $contigsfn $readsfn $readmap $contigmap $newcontigs

outfile=$4
echo 'Integrating scaffolds into assembly'
java -cp $BINDIR scaffolding.Scaffold $contigsfn $newcontigs $outfile
