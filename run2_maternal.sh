rm /scratch/groups/mschatz1/mkirsche/ultralong/ccs/maternal*_usefulmap2.paf
./run2.sh /scratch/groups/mschatz1/mkirsche/ultralong/ccs/maternal_and_unknown.contigs.mmpoa.fa /scratch/groups/mschatz1/mkirsche/ultralong/rel2_200kplus.fastq /scratch/groups/mschatz1/mkirsche/ultralong/ccs/rel2_200kplus_ccs_mat.paf /scratch/groups/mschatz1/mkirsche/ultralong/ccs/maternal_and_unknown.contigs.mmpoa.scaffolds2.fa | tee maternal.log;

rm /scratch/groups/mschatz1/mkirsche/ultralong/ccs/maternal*_usefulmap2.paf
./run2.sh /scratch/groups/mschatz1/mkirsche/ultralong/ccs/maternal_and_unknown.contigs.mmpoa.scaffolds2.fa /scratch/groups/mschatz1/mkirsche/ultralong/rel2_200kplus.fastq /scratch/groups/mschatz1/mkirsche/ultralong/ccs/rel2_200kplus_ccs_mat.paf /scratch/groups/mschatz1/mkirsche/ultralong/ccs/maternal_and_unknown.contigs.mmpoa.scaffolds3.fa | tee maternal2.log;

rm /scratch/groups/mschatz1/mkirsche/ultralong/ccs/maternal*_usefulmap2.paf
./run2.sh /scratch/groups/mschatz1/mkirsche/ultralong/ccs/maternal_and_unknown.contigs.mmpoa.scaffolds3.fa /scratch/groups/mschatz1/mkirsche/ultralong/rel2_200kplus.fastq /scratch/groups/mschatz1/mkirsche/ultralong/ccs/rel2_200kplus_ccs_mat.paf /scratch/groups/mschatz1/mkirsche/ultralong/ccs/maternal_and_unknown.contigs.mmpoa.scaffolds4.fa | tee maternal2.log;

java -cp ~/hashing/CCS/assembly_eval/ AssemblyStats /scratch/groups/mschatz1/mkirsche/ultralong/ccs/maternal_and_unknown.contigs.mmpoa.scaffolds4.fa > maternal_afterstats2.txt;
./delta.sh /scratch/groups/mschatz1/mkirsche/ultralong/ccs/maternal_and_unknown.contigs.mmpoa.scaffolds4.fa /scratch/groups/mschatz1/mkirsche/ultralong/ccs/after_maternal2;

