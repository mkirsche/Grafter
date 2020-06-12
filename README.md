# Grafter
Grafter is used for scaffolding assemblies using alignments of ultralong reads which span multiple contigs.  It does so by finding reads which align to multiple contigs, constructing a graph of pairs of contigs which have one or more reads spanning them, and iteratively building up scaffolds based on the edge weights in this graph.

## Dependencies

A java runtime environment is required to compile and run Grafter.

## Compilation

`javac src/*.java`

## Running

``
java -cp src Main <args>
``

A usage menu can be viewed by running the program with no arguments, and is included below:

```

Grafter: A tool for scaffolding assemblies using alignments of ultralong reads

Usage: java -cp src Main [args]
  Example: java -cp src Main aln_fn=aln.paf fasta_fn=contigs.fasta read_fn=reads.fastq
    read_map_file=useful_reads.paf contig_map_file=useful_contigs.paf out_file=out.fasta

Required args:
  aln_fn          (String) - a file containing the alignments of ultralong reads to contigs
  fasta_fn        (String) - the contigs in FASTA format
  read_fn         (String) - the ultralong reads in FASTQ format
  outputbroken    (String) - where to output broken contigs
  read_map_file   (String) - where to output sequences of relevant reads
  contig_map_file (String) - Where to output sequences of relevant contigs
  out_file        (String) - the name of the file to output the scaffolded contigs to

Optional args:
  max_hanging      (int)    [1000]  - the maximum amount by which the end of a contig can exceed the alignment and still be joined
  minq             (int)    [40]    - the minimum quality of alignments needed to be kept
  min_read_supp    (int)    [1]     - number of reads which need to support an edge
  min_weight_supp  (float)  [15000] - total weight required for a pair of contigs to be joined
  min_weight       (float)  [1000]    - weight required for an overlap to count
  min_length       (int)    [3000]  - minimum length of alignments required on each read
  full_out_gfa_fn  (String) [none]  - where to write the scaffold graph in GFA format
  joins_out__gfa_fn (String) [none]  - where to write the scaffold graph in GFA format
  --break                          - allows original contigs to be broken
  --reuse_relevant_seqs            - reuse files with sequences of relevant reads and contigs

```

This produces a FASTA file with all scaffolds which consisted of 2 or more contigs, as well as any subcontigs which did not get rejoined if misassembly detection is enabled.  To obtain an updated assembly with contigs replaced by the scaffolds they are included in, run the following:

``
java -cp src StitchFasta [assembly_contigs.fa] [scaffolds.fa] [assembly_scaffolds.fa]
``


## Inputs

Grafter requires three inputs:

* A FASTA file containing the contigs from the original assembly
* A FASTQ file containing ultralong reads
* A PAF file containing alignments of the ultralong reads to the contigs
