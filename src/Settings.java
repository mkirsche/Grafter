
public class Settings {
	
	// The number of reads required to support the joining of two contigs
	static int MIN_READ_SUPPORT = 1;
	static double MIN_WEIGHT_SUPPORT = 15000;
	
	static double MAX_HANGING_PROP = 0.02;
	static int MAX_HANGING = 1000;
	static boolean VERBOSE = true;
	static boolean OUTPUT_BROKEN = false;
	static boolean ALLOW_BREAKS = false;
	
	static boolean PRINT_ORIENT = true;
	
	static double MIN_WEIGHT = 1000;
	
	static int MIN_ALIGNMENT_LENGTH = 3000;
	
	// The minimum alignment quality for an alignment to be kept
	static int MIN_QUALITY = 40;
	
	static int MAX_GAP = 10000;
	
	// File name of PAF file containing read-to-contig mappings
	static String pafFn = "";
	
	// File name of FASTA file containing contigs
	static String fastaFn = "";
	
	// File name of FASTQ file containing reads
	static String readFn = "";
	
	// Where to output sequences of relevant reads
	static String relevantReadSequenceFile = "";
	
	// Where to output sequences of relevant contigs
	static String relevantContigSequenceFile = "";
	
	// Where to output contigs that are formed
	static String outFn = "";
	
	// Where to output information about broken contigs because of possible misassemblies, if using that mode
	static String brokenOutputFile = "";
	
	// File name of GFA file containing an existing assembly graph
	static String graphFn = "";
	
	// Where to output the scaffold graph as a GFA file
	static String fullOutGfaFn = "";
	
	// Where to output the joins as a GFA file
	static String joinsOutGfaFn = "";
	
	// Where to output read intervals being used in the joins
	static String readMetadataFn = "";
	
	static boolean reuseRelevantSeqs = false;
	
	/*
	 * Parse command line arguments
	 */
	static void parseArgs(String[] args)
	{
		for(String arg : args)
		{
			int equalsIdx = arg.indexOf("=");
			if(equalsIdx == -1)
			{
				if(arg.toLowerCase().endsWith("break"))
				{
					Settings.ALLOW_BREAKS = true;
				}
				if(arg.toLowerCase().endsWith("reuse_relevant_seqs"))
				{
					Settings.reuseRelevantSeqs = true;
				}
			}
			else
			{
				String field = arg.substring(0, equalsIdx);
				String val = arg.substring(1 + equalsIdx);
				if(field.equalsIgnoreCase("aln_fn"))
				{
					Settings.pafFn = val;
				}
				if(field.equalsIgnoreCase("fasta_fn"))
				{
					Settings.fastaFn = val;
				}
				if(field.equalsIgnoreCase("read_fn"))
				{
					Settings.readFn = val;
				}
				if(field.equalsIgnoreCase("graph_fn"))
				{
					Settings.graphFn = val;
				}
				if(field.equalsIgnoreCase("full_out_gfa_fn"))
				{
					Settings.fullOutGfaFn = val;
				}
				if(field.equalsIgnoreCase("joins_out_gfa_fn"))
				{
					Settings.joinsOutGfaFn = val;
				}
				if(field.equalsIgnoreCase("read_metadata_fn"))
				{
					Settings.readMetadataFn = val;
				}
				if(field.equalsIgnoreCase("read_map_file"))
				{
					Settings.relevantReadSequenceFile = val;
				}
				if(field.equalsIgnoreCase("outputbroken"))
				{
					Settings.OUTPUT_BROKEN = true;
					Settings.brokenOutputFile = val;
				}
				if(field.equalsIgnoreCase("contig_map_file"))
				{
					Settings.relevantContigSequenceFile = val;
				}
				if(field.equalsIgnoreCase("out_file"))
				{
					Settings.outFn = val;
				}
				if(field.equalsIgnoreCase("max_hanging"))
				{
					Settings.MAX_HANGING = Integer.parseInt(val);
				}
				if(field.equalsIgnoreCase("minq"))
				{
					Settings.MIN_QUALITY = Integer.parseInt(val);
				}
				if(field.equalsIgnoreCase("min_weight_supp"))
				{
					Settings.MIN_WEIGHT_SUPPORT = Double.parseDouble(val);
				}
				if(field.equalsIgnoreCase("min_weight"))
				{
					Settings.MIN_WEIGHT = Double.parseDouble(val);
				}
				if(field.equalsIgnoreCase("min_length"))
				{
					Settings.MIN_ALIGNMENT_LENGTH = Integer.parseInt(val);
				}
				if(field.equalsIgnoreCase("min_read_supp"))
				{
					Settings.MIN_READ_SUPPORT = Integer.parseInt(val);
				}
			}
		}
		if(Settings.pafFn.length() == 0 || Settings.fastaFn.length() == 0)
		{
			usage();
			System.exit(1);
		}
		if(Settings.relevantReadSequenceFile.length() == 0 || Settings.relevantContigSequenceFile.length() == 0)
		{
			usage();
			System.exit(1);
		}
		if(Settings.readFn.length() == 0 || Settings.outFn.length() == 0)
		{
			usage();
			System.exit(1);
		}
	}
	
	/*
	 * Print a usage menu
	 */
	static void usage()
	{
		System.out.println();
		System.out.println("Grafter: A tool for scaffolding assemblies using alignments of ultralong reads");
		System.out.println();
		System.out.println("Usage: java -cp src Main [args]");
		System.out.println("  Example: java -cp src Main aln_fn=aln.paf fasta_fn=contigs.fasta read_fn=reads.fastq");
		System.out.println("    read_map_file=useful_reads.paf contig_map_file=useful_contigs.paf out_file=out.fasta");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  aln_fn          (String) - a file containing the alignments of ultralong reads to contigs");
		System.out.println("  fasta_fn        (String) - the contigs in FASTA format");
		System.out.println("  read_fn         (String) - the ultralong reads in FASTQ format");
		System.out.println("  outputbroken    (String) - where to output broken contigs");
		System.out.println("  read_map_file   (String) - where to output sequences of relevant reads");
		System.out.println("  contig_map_file (String) - Where to output sequences of relevant contigs");
		System.out.println("  out_file        (String) - the name of the file to output the scaffolded contigs to");
		System.out.println();
		System.out.println("Optional args:");
		System.out.println("  max_hanging      (int)    [1000]  - the maximum amount by which the end of a contig can exceed the alignment and still be joined");
		System.out.println("  minq             (int)    [40]    - the minimum quality of alignments needed to be kept");
		//System.out.println("  graph_fn         (String) [none] - a GFA file containing an assembly graph, causing only alignments which are validated by the graph to be kept");
		System.out.println("  min_read_supp    (int)    [1]     - number of reads which need to support an edge");
		System.out.println("  min_weight_supp  (float)  [15000] - total weight required for a pair of contigs to be joined");
		System.out.println("  min_weight       (float)  [1000]    - weight required for an overlap to count");
		System.out.println("  min_length       (int)    [3000]  - minimum length of alignments required on each read");
		System.out.println("  full_out_gfa_fn  (String) [none]  - where to write the scaffold graph in GFA format");
		System.out.println("  joins_out_gfa_fn (String) [none]  - where to write the scaffold graph in GFA format");
		System.out.println("  read_metadata_fn (String) [none]  - where to write the reads being used as a tsv");
		System.out.println("  --break                           - allows original contigs to be broken");
		System.out.println("  --reuse_relevant_seqs             - reuse files with sequences of relevant reads and contigs");
		System.out.println();
	}
}
