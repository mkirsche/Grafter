minimappath=/scratch/groups/mschatz1/mkirsche/github/minimap2/minimap2
genomepath=/scratch/groups/mschatz1/mkirsche/svdb/genome.fa

assembly=$1
outprefix=$2

$minimappath -t 32 -x asm5 -c --cs $genomepath $assembly > $outprefix'.paf'
python paf2delta.py $outprefix'.paf'
