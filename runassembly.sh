rm /scratch/groups/mschatz1/mkirsche/ultralong/ccs/maternal*_usefulmap2.paf
rm /scratch/groups/mschatz1/mkirsche/ultralong/*_usefulmap2.paf
./run2.sh /scratch/groups/mschatz1/mkirsche/ultralong/ccs/maternal_and_unknown.contigs.mmpoa.fa /scratch/groups/mschatz1/mkirsche/ultralong/rel2wt_50kplus.ctg.fa /scratch/groups/mschatz1/mkirsche/ultralong/ccs/MaternalToAssembly.paf /scratch/groups/mschatz1/mkirsche/ultralong/ccs/maternal_and_unknown.contigs.mmpoa.scaffoldsassemblygraph.fa | tee maternal.log;

java -cp ~/hashing/CCS/assembly_eval/ AssemblyStats /scratch/groups/mschatz1/mkirsche/ultralong/ccs/maternal_and_unknown.contigs.mmpoa.scaffoldsassemblygraph.fa | tee maternal_afterstatsassemblygraph.txt;

./delta.sh /scratch/groups/mschatz1/mkirsche/ultralong/ccs/maternal_and_unknown.contigs.mmpoa.scaffoldsassemblygraph.fa /scratch/groups/mschatz1/mkirsche/ultralong/ccs/after_maternalassemblygraph;

