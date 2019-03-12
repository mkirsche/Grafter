minimappath=/scratch/groups/mschatz1/mkirsche/github/minimap2/minimap2
genomepath=/home-3/mkirsche@jhu.edu/csdata/genome.fa

assembly=$1
outprefix=$2

$minimappath -x asm5 -c --cs $genomepath $assembly > $outprefix'.paf'
python paf2delta.py $outprefix'.paf'
