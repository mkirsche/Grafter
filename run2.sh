contigsfn=$1
readsfn=$2
paffn=$3
outfile=$4
breakcontigs=$5

minimappath=/scratch/groups/mschatz1/mkirsche/github/minimap2/minimap2

if [ ! -f $paffn ]; then
    echo "PAF file not found - generating it"
    $minimappath -t 32 -k 19 -w 19 $contigsfn $readsfn  > $paffn
fi

BINDIR=`dirname $(readlink -f "$0")`

echo 'Compiling'
javac $BINDIR/scaffolding/*.java

readmap=$readsfn'_usefulmap2.paf'
contigmap=$contigsfn'_usefulmap2.paf'
newcontigs=$contigsfn'_newcontigs2.paf'

if [ "$breakcontigs" -eq "0" ]; then
    java -cp $BINDIR scaffolding.IncludeContained $paffn $contigsfn $readsfn $readmap $contigmap $newcontigs
    echo 'Scaffolds output to '$newcontigs
    echo 'Integrating scaffolds into assembly'
    java -cp $BINDIR scaffolding.StitchFasta $contigsfn $newcontigs $outfile
    echo 'Final assembly output to '$outfile
    exit;
fi

#brokencontigs=$contigsfn'.broken'
#echo 'Scaffolding'
#java -cp $BINDIR scaffolding.IncludeContained $paffn $contigsfn $readsfn $readmap $contigmap $newcontigs '--break --outputbroken='$brokencontigs
java -cp $BINDIR scaffolding.IncludeContained $paffn $contigsfn $readsfn $readmap $contigmap $newcontigs '--break --outputbroken='$brokencontigs
#echo 'Scaffolds with original mappings output to '$newcontigs

#brokenpaffn=$paffn'.broken.paf'
#readmapbroken=$readsfn'_usefulmap2.paf.broken'
#contigmapbroken=$contigsfn'_usefulmap2.paf.broken'
#newcontigsbroken=$contigsfn'_newcontigs2.paf.broken'
#if [ ! -f $brokenpaffn ]; then
#    echo "Broken PAF file not found - generating it"
#    $minimappath -t 32 -k 19 -w 19 $brokencontigs $readsfn  > $brokenpaffn
#fi
#java -cp $BINDIR scaffolding.IncludeContained $brokenpaffn $brokencontigs $readsfn $readmapbroken $contigmapbroken $newcontigsbroken
echo 'Scaffolds with updated mappings output to '$newcontigs

echo 'Integrating scaffolds into assembly'
java -cp $BINDIR scaffolding.StitchFasta $contigsfn $newcontigs $outfile
echo 'Final assembly output to '$outfile
