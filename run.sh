contigsfn=$1
readsfn=$2
paffn=$3

BINDIR=`dirname $(readlink -f "$0")`

javac $BINDIR/scaffolding/*.java

usefulpaf=$paffn'_useful.paf'
java -cp $BINDIR scaffolding.FindUsefulScaffoldingAlignments $paffn $usefulpaf

readmap=$readsfn'_usefulmap.paf'
contigmap=$contigsfn'_usefulmap.paf'
newcontigs=$contigsfn'_newcontigs.paf'
java -cp $BINDIR scaffolding.Scaffold $usefulpaf $contigsfn $readsfn $readmap $contigmap $newcontigs

outfile=$4
java -cp $BINDIR scaffolding.Scaffold $contigsfn $newcontigs $outfile
