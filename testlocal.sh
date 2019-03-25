./sync.sh
javac scaffolding/*.java

paffn=/home/mkirsche/eclipse-workspace/Ultralong/rel2_200kplus_ccs_mat.paf
fastafn=/home/mkirsche/eclipse-workspace/Ultralong/maternal_and_unknown.contigs.mmpoa.fa
readfn=/home/mkirsche/eclipse-workspace/Ultralong/rel2_200kplus.fastq
readmapfn=/home/mkirsche/eclipse-workspace/Ultralong/readmap_maternal.txt
contigmapfn=/home/mkirsche/eclipse-workspace/Ultralong/contigmap_maternal.txt
newcontigsfn=/home/mkirsche/eclipse-workspace/Ultralong/new_contigs.fa
outputfn=/home/mkirsche/eclipse-workspace/Ultralong/maternalscaffolds.fa
logfn=/home/mkirsche/eclipse-workspace/Ultralong/maternallog.txt

java -Xmx8192m scaffolding.IncludeContained $paffn $fastafn $readfn $readmapfn $contigmapfn $newcontigsfn 2>&1 | tee $logfn
java -Xmx4096m scaffolding.StitchFasta $fastafn $newcontigsfn $outputfn 2>&1 | tee -a $logfn
cat $logfn | grep -A 5 'Assembly Stats'


