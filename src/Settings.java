
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
	
	static double MIN_WEIGHT = 10;
	
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
	static String outGfaFn = "";
	
	static boolean reuseRelevantSeqs = false;
}
