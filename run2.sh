contigsfn=$1
readsfn=$2
paffn=$3
outfile=$4

minimappath=/scratch/groups/mschatz1/mkirsche/github/minimap2/minimap2

if [ ! -f $paffn ]; then
    echo "PAF file not found - generating it"
    $minimappath -t 8 $contigsfn $readsfn  > $paffn
fi

BINDIR=`dirname $(readlink -f "$0")`

echo 'Compiling'
javac $BINDIR/scaffolding/*.java

readmap=$readsfn'_usefulmap2.paf'
contigmap=$contigsfn'_usefulmap2.paf'
newcontigs=$contigsfn'_newcontigs2.paf'
echo 'Scaffolding'
java -cp $BINDIR scaffolding.IncludeContained $paffn $contigsfn $readsfn $readmap $contigmap $newcontigs
echo 'Scaffolds output to '$newcontigs

echo 'Integrating scaffolds into assembly'
java -cp $BINDIR scaffolding.StitchFasta $contigsfn $newcontigs $outfile
echo 'Final assembly output to '$outfile
