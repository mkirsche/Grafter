rm /scratch/groups/mschatz1/mkirsche/ultralong/ccs/paternal*_usefulmap2.paf
./run2.sh /scratch/groups/mschatz1/mkirsche/ultralong/ccs/paternal_and_unknown.contigs.mmpoa.fa /scratch/groups/mschatz1/mkirsche/ultralong/rel2_200kplus.fastq /scratch/groups/mschatz1/mkirsche/ultralong/ccs/rel2_200kplus_ccs_pat.paf /scratch/groups/mschatz1/mkirsche/ultralong/ccs/paternal_and_unknown.contigs.mmpoa.scaffolds2.fa | tee paternal.log;
./delta.sh /scratch/groups/mschatz1/mkirsche/ultralong/ccs/paternal_and_unknown.contigs.mmpoa.scaffolds2.fa /scratch/groups/mschatz1/mkirsche/ultralong/ccs/after_paternal2;
java -cp ~/hashing/CCS/assembly_eval/ AssemblyStats /scratch/groups/mschatz1/mkirsche/ultralong/ccs/paternal_and_unknown.contigs.mmpoa.scaffolds2.fa > paternal_afterstats2.txt;
