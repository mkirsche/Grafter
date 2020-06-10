import java.util.*;

import java.io.*;

public class IncludeContained {
	
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
		System.out.println("Ultralong Scaffolding - including contained nodes");
		System.out.println("Usage: java -cp src IncludeContained [args]");
		System.out.println("  Example: java -cp src IncludeContained aln_fn=aln.paf fasta_fn=contigs.fasta read_fn=reads.fastq");
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
		System.out.println("Optional args");
		System.out.println("  max_hanging (int)    [1000] - the maximum amount by which the end of a contig can exceed the alignment and still be joined");
		System.out.println("  minq        (int)    [40]   - the minimum quality of alignments needed to be kept");
		System.out.println("  graph_fn    (String) [none] - a GFA file containing an assembly graph, causing only alignments which are validated by the graph to be kept");
		System.out.println("  --break                     - allows original contigs to be broken");
		System.out.println("  --reuse_relevant_seqs       - reuse files with sequences of relevant reads and contigs");
		System.out.println();
	}
	
public static void main(String[] args) throws Exception
{
	// Read in command line parameters
	parseArgs(args);
	
	Scanner input = new Scanner(new FileInputStream(new File(Settings.pafFn)));
	PrintWriter out = new PrintWriter(new File(Settings.outFn));
	
	// Read in alignments and bucket by which read was aligned
	HashMap<String, ArrayList<SortablePafAlignment>> alignmentsPerRead = new HashMap<>();
	while(input.hasNext())
	{
		String line = input.nextLine();
		SortablePafAlignment cur = new SortablePafAlignment(line);
		
		double curThreshold = Math.min(.2 * cur.readLength, Settings.MIN_ALIGNMENT_LENGTH);
		
		// Filter out short alignments
		if(cur.readEnd - cur.readStart < curThreshold)
		{
			continue;
		}
		
		// Filter out low-quality alignments
		if(cur.mapq < Settings.MIN_QUALITY)
		{
			continue;
		}
		
		String readName = cur.readName;
		
		ReadUtils.addToMap(alignmentsPerRead, readName, cur);
	}
	
	// If performing misassembly correction, find a list of breakpoints based on novel adjacencies
	ArrayList<CorrectMisassemblies.NovelAdjacency> corrections = new ArrayList<CorrectMisassemblies.NovelAdjacency>();
	if(Settings.ALLOW_BREAKS)
	{
		corrections = CorrectMisassemblies.findMisassemblies(alignmentsPerRead);
		if(Settings.VERBOSE)
		{
			for(CorrectMisassemblies.NovelAdjacency na : corrections)
			{
				System.err.println(na);
			}
			System.err.println("Number of misassemblies: " + corrections.size());
		}
	}
	
	// Perform splitting as needed and remap reads to broken contigs
	HashSet<String> contigNames = new HashSet<String>();
	
	CorrectMisassemblies.ContigBreaker splitter = new CorrectMisassemblies.ContigBreaker(corrections, contigNames);
	
	System.err.println("Number of breaks: " + splitter.numBreaks);
	
	alignmentsPerRead = CorrectMisassemblies.remapAll(splitter, alignmentsPerRead);
	
	/*
	 * Get chains of unique mappings to reads and keep track of contigs/reads involved in them
	 */
	
	// Map from read to alignment chains it's involved in
	HashMap<String, ArrayList<ArrayList<SortablePafAlignment>>> chainsPerRead = new HashMap<>();
	
	// Set of read names involved in alignment chains
	HashSet<String> readNames = new HashSet<String>();
	
	// Iterate over reads, find chains of read
	for(String s : alignmentsPerRead.keySet())
	{
		if(alignmentsPerRead.get(s).size() == 1)
		{
			continue;
		}
		
		ArrayList<ArrayList<SortablePafAlignment>> chains = AlignmentGatherer.getUniqueMatches(alignmentsPerRead.get(s));
		
		if(chains.size() > 0)
		{
			for(ArrayList<SortablePafAlignment> l : chains)
			{
				for(SortablePafAlignment spa : l)
				{
					contigNames.add(spa.contigName);
				}
			}
			chainsPerRead.put(s, chains);
			readNames.add(s);
		}
	}
	
	/*
	 * Output broken assembly
	 */
	if(Settings.OUTPUT_BROKEN && corrections.size() > 0)
	{
		System.err.println("Outputting broken assembly");
		splitter.outputBrokenAssembly(Settings.fastaFn, Settings.brokenOutputFile);
	}
	
	/*
	 * Get sequences of relevant contigs/reads for merging
	 */
	
	// Relevant reads
	HashMap<String, String> readSequences = new HashMap<>(), contigSequences = new HashMap<>();
	if(!Settings.reuseRelevantSeqs || (readSequences = ReadUtils.readMap(Settings.relevantReadSequenceFile)).size() == 0)
	{
		System.err.println("Filtering reads");
		if(Settings.readFn.endsWith(".fa") || Settings.readFn.endsWith(".fasta"))
		{
			readSequences = ReadUtils.getFastaMap(Settings.readFn, readNames);
		}
		else
		{
			readSequences = ReadUtils.getFastqMap(Settings.readFn, readNames);
		}
		ReadUtils.writeMap(Settings.relevantReadSequenceFile, readSequences);
	}
	
	// Relevant contigs
	if(!Settings.reuseRelevantSeqs || (contigSequences = ReadUtils.readMap(Settings.relevantContigSequenceFile)).size() == 0)
	{
		System.err.println("Filtering contigs");
		contigSequences = ReadUtils.getFastaMap(Settings.fastaFn, contigNames);
		ReadUtils.writeMap(Settings.relevantContigSequenceFile, contigSequences);
	}
	
	if(Settings.VERBOSE)
	{
		System.err.println("Split contigs:\n" +splitter.subcontigMap.keySet());
	}
	
	// Adjust contig sequence map based on any splitting that happened
	ArrayList<String> keys = new ArrayList<String>();
	keys.addAll(contigSequences.keySet());
	
	for(String s : keys)
	{
		if(splitter.breakSequence(s, contigSequences.get(s)))
		{
			contigSequences.remove(s);
		}
	}
	for(String splitContigName : splitter.sequenceMap.keySet())
	{
		contigSequences.put(splitContigName, splitter.sequenceMap.get(splitContigName));
	}
	
	/*
	 * Compute k-mer frequencies across different reads which will be used to get better measures of overlap for graph-building
	 */
	System.err.println("Initializing frequency map for contig kmers");
	ContigKmerFrequencyMap freq = new ContigKmerFrequencyMap();
	
	// Add k-mers to index for overall counts and lengths of sequences
	System.err.println("Adding contig kmer frequencies");
	for(String s : contigSequences.keySet())
	{
		freq.addKmerCount(s, contigSequences.get(s));
	}
	
	// Index the k-mer counts of each sequence with a cumulative sum array for faster queries
	System.err.println("Indexing contig kmer frequencies");
	for(String s : contigSequences.keySet())
	{
		freq.addSumArray(s, contigSequences.get(s));
	}
	
	/*
	 * Add edges to the scaffold graph based on the chains of alignments
	 */
	ScaffoldGraph sg = new ScaffoldGraph();
	System.err.println("Joining contigs");
	int numMerged = 0;
	for(String readName : chainsPerRead.keySet())
	{
		ArrayList<ArrayList<SortablePafAlignment>> allChains = chainsPerRead.get(readName);
		for(ArrayList<SortablePafAlignment> chain : allChains)
		{
			addEdges(sg, chain, freq);
		}
	}
	
	/*
	 * Add a dummy edge between split contigs to give them the opportunity to be rejoined if they don't get joined with other things
	 */
	readSequences.put("undosplit", "A");
	for(String s : splitter.subcontigMap.keySet())
	{
		ArrayList<CorrectMisassemblies.ContigBreaker.Subcontig> subs = splitter.subcontigMap.get(s);
		int numSubcontigs = subs.size();
		for(int i = 0; i<numSubcontigs-1; i++)
		{
			sg.addEdge(subs.get(i).name, subs.get(i+1).name, "undosplit", 0, 0, 0, false, true, 1);
		}
	}
	
	// Output the contig overlap graph
	OutputScaffolds.outputGfa("gfaFn", sg, contigSequences);
	
	/*
	 * Run scaffolding on the graph
	 */
	ScaffoldGraph.Scaffolding results = sg.globalScaffolding();
	HashMap<String, ArrayDeque<String>> scaffoldContigs = results.scaffoldContigs;
	HashMap<String, ArrayDeque<ScaffoldGraph.Alignment>> scaffoldEdges = results.scaffoldEdges;
	HashSet<String> usedContigs = results.usedContigs;
	
	numMerged = results.numMerged;
	
	/*
	 * Output all scaffolds consisting of multiple contigs
	 */
	for(String s : scaffoldContigs.keySet())
	{
		String headerLine = OutputScaffolds.createHeaderLine(scaffoldContigs.get(s), splitter);
		if(Settings.VERBOSE)
		{
			System.err.println(headerLine);
		}
		out.println(headerLine);
		
		String seq = merge(scaffoldContigs.get(s), scaffoldEdges.get(s), readSequences, contigSequences);
		
		out.println(seq);
	}
	
	/*
	 * Output subcontigs which were not rejoined to anything
	 */
	for(String s : splitter.subcontigMap.keySet())
	{
		ArrayList<CorrectMisassemblies.ContigBreaker.Subcontig> cur = splitter.subcontigMap.get(s);
		for(CorrectMisassemblies.ContigBreaker.Subcontig sc : cur)
		{
			if(!usedContigs.contains(sc.name))
			{
				out.println(">" + sc.name + " " + sc.oldName);
				out.println(splitter.sequenceMap.get(sc.name));
			}
		}
	}
	System.err.println("Number of joins: " + numMerged);
	
	if(Settings.PRINT_ORIENT)
	{
		OutputScaffolds.printOrientations(scaffoldEdges, splitter);
	}
	
}

/*
 * Merges contigs together based on the alignments in a path of a scaffold graph
 */
static String merge(ArrayDeque<String> contigs, ArrayDeque<ScaffoldGraph.Alignment> als, HashMap<String, String> readMap, HashMap<String, String> relevantContigs)
{
	StringBuilder res = new StringBuilder();
	boolean first = true;
	for(ScaffoldGraph.Alignment spa : als)
	{
		if(first)
		{
			first = false;
			String curSeq = relevantContigs.get(contigs.peekFirst());
			if(Settings.VERBOSE)
			{
				System.err.println(contigs.peekFirst()+" "+curSeq.length());
			}
			if(spa.myContigPrefix) curSeq = ReadUtils.reverseComplement(curSeq);
			res.append(curSeq);
		}
		if(Settings.VERBOSE)
		{
			System.err.println(spa.to+" "+spa.myContigPrefix+" "+spa.theirContigPrefix+" "+spa.myReadEnd+" "+spa.theirReadStart+" "+spa.strand+" "+spa.read + " " + spa.weight);
		}
		int overlap = 0;
		if(spa.myReadEnd < spa.theirReadStart)
		{
			String readSeq = readMap.get(spa.read);
			if(!spa.theirContigPrefix)
			{
				readSeq = ReadUtils.reverseComplement(readSeq);
			}
			res.append(readSeq.substring(spa.myReadEnd, spa.theirReadStart));
		}
		else
		{
			overlap = spa.myReadEnd - spa.theirReadStart;
		}
		
		String curSeq = relevantContigs.get(spa.to);
				
		if(Settings.VERBOSE)
		{
			System.err.println(spa.to + " " + curSeq.length() + " " +overlap);
		}
		
		if(!spa.theirContigPrefix)
		{
			curSeq = ReadUtils.reverseComplement(curSeq);
		}
		res.append(curSeq.substring(overlap));
	}
	return res.toString();
}

/*
 * Add edges to a scaffold graph based on a chain of alignments to the same read
 */
static void addEdges(ScaffoldGraph sg, ArrayList<SortablePafAlignment> als, ContigKmerFrequencyMap freq)
{
	SortablePafAlignment last = null;
	boolean lastReversed = false;
	for(int i = 0; i<als.size(); i++)
	{
		SortablePafAlignment spa = als.get(i);
		
		boolean curReversed = false;
		
		boolean[] cse = AlignmentGatherer.contigStartEnd(spa);
		boolean[] rse = AlignmentGatherer.readStartEnd(spa);
		
		if(cse[0] && cse[1])
		{
			// Entire contig aligns - have to look at strand in alignment
			if(spa.strand == '-')
			{
				curReversed = true;
			}
		}
		else if(cse[0])
		{
			// Beginning of contig - we would expect suffix of read if same strand
			if(rse[0] && rse[1])
			{
				if(spa.strand == '-')
				{
					curReversed = true;
				}
			}
			else if(rse[0])
			{
				curReversed = true;
			}
		}
		else
		{
			// End of contig - we would expect prefix of read if same strand
			if(rse[0] && rse[1])
			{
				if(spa.strand == '-')
				{
					curReversed = true;
				}
			}
			else if(rse[1])
			{
				curReversed = true;
			}
		}
		
		if(last != null && !last.contigName.equals(spa.contigName))
		{
			int overlap = last.readEnd - spa.readStart;
			if(overlap <= spa.contigLength && overlap <= last.contigLength 
					&& overlap < .9 * Math.min(last.readEnd-last.readStart, spa.readEnd - spa.readStart))
			{
				double lastLength = last.contigEnd - last.contigStart;
				double curLength = spa.contigEnd - spa.contigStart;
				double weight = 2 * lastLength * curLength / (lastLength + curLength);
				double avgFreq1 = freq.getAverageFrequency(last.contigName, last.contigStart-1, last.contigEnd-1);
				double avgFreq2 = freq.getAverageFrequency(spa.contigName, spa.contigStart-1, spa.contigEnd-1);
				double penalty = CorrectMisassemblies.harmonicMean(avgFreq1, avgFreq2);
				System.out.println("repeat penalty: " + last.contigName+" "+spa.contigName+" "+penalty);
				weight /= penalty;
				if(weight >= Settings.MIN_WEIGHT)
				{
					System.out.println("repeat penalty: " + last.contigName+" "+spa.contigName+" "+penalty);
					sg.addEdge(last.contigName, spa.contigName, last.readName, last.readEnd, spa.readStart, spa.readLength, lastReversed, !curReversed, weight);
				}
			}
		}
		
		last = spa;
		lastReversed = curReversed;
	}
}


}
