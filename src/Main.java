import java.util.*;

import java.io.*;

public class Main {
	
public static void main(String[] args) throws Exception
{
	// Read in command line parameters
	Settings.parseArgs(args);
	
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
	if(Settings.fullOutGfaFn.length() > 0)
	{
		OutputScaffolds.outputGfa(Settings.fullOutGfaFn, sg, contigSequences);
	}
	
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

    /*
	 * Output GFA based on joins
	 */
	if(Settings.joinsOutGfaFn.length() > 0)
	{
		PrintWriter joinsOut = new PrintWriter(new File(Settings.joinsOutGfaFn));
		joinsOut.println("H\t1.0");
		for(String s : contigSequences.keySet())
		{
			joinsOut.println("S\t" + s + "\t*\tLN:" + contigSequences.get(s).length());
		}
		for(String contigKey : scaffoldEdges.keySet())
		{
			ArrayDeque<String> curContigNames = scaffoldContigs.get(contigKey);
			ArrayDeque<ScaffoldGraph.Alignment> curContigEdges = scaffoldEdges.get(contigKey);
			String from = curContigNames.pollFirst();
			for(ScaffoldGraph.Alignment aln : curContigEdges)
			{
				String to = curContigNames.pollFirst();
				char fromStrand = aln.myContigPrefix ? '-' : '+';
				char toStrand = aln.theirContigPrefix ? '-' : '+';
				out.printf("%s\t%s\t%s\t%s\t%s\t%s\n", "L", from, fromStrand, to, toStrand, "*");
				from = to;
			}
		}
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
				System.err.println("repeat penalty: " + last.contigName+" "+spa.contigName+" "+penalty);
				weight /= penalty;
				if(weight >= Settings.MIN_WEIGHT)
				{
					System.err.println("repeat penalty: " + last.contigName+" "+spa.contigName+" "+penalty);
					sg.addEdge(last.contigName, spa.contigName, last.readName, last.readEnd, spa.readStart, spa.readLength, lastReversed, !curReversed, weight);
				}
			}
		}
		
		last = spa;
		lastReversed = curReversed;
	}
}


}
