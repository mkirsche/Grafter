./run2.sh /scratch/groups/mschatz1/mkirsche/ultralong/ccs/maternal_and_unknown.contigs.mmpoa.intra.chimera.broken.fa /scratch/groups/mschatz1/mkirsche/ultralong/rel2_200kplus.fastq /scratch/groups/mschatz1/mkirsche/ultralong/ccs/brokenmaternal.paf /scratch/groups/mschatz1/mkirsche/ultralong/ccs/maternal_and_unknown.contigs.mmpoa.broken.fa 0 2>&1 | tee maternal.log;

java -cp ~/hashing/CCS/assembly_eval/ AssemblyStats /scratch/groups/mschatz1/mkirsche/ultralong/ccs/maternal_and_unknown.contigs.mmpoa.broken.fa | tee maternal_afterstatsbroken.txt;

./delta.sh /scratch/groups/mschatz1/mkirsche/ultralong/ccs/maternal_and_unknown.contigs.mmpoa.broken.fa /scratch/groups/mschatz1/mkirsche/ultralong/ccs/after_maternalbroken;

