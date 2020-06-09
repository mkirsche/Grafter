import java.util.*;

import java.io.*;

public class IncludeContained {
	
	// The number of reads required to support the joining of two contigs
	static int minReadSupport = 1;
	static double minWeightSupport = 15000;
	
	static double maxHangingProportion = 0.02;
	static int maxHanging = 1000;
	static boolean fileMap = false;
	static boolean verbose = true;
	static boolean outputBroken = false;
	static boolean allowBreaks = false;
	
	static boolean printOrientations = true;
	
	static double minWeight = 10;
	
	static int minAlignmentLength = 3000;
	
	static int minQuality = 40;
	
	static String pafFn = "", fastaFn = "", readFn = "";
	static String readMapFile = "", contigMapFile = "", outFn = "";
	static String brokenOutputFile = "", graphFn = "";
	
	static void parseArgs(String[] args)
	{
		for(String arg : args)
		{
			int equalsIdx = arg.indexOf("=");
			if(equalsIdx == -1)
			{
				if(arg.toLowerCase().endsWith("break"))
				{
					allowBreaks = true;
				}
			}
			else
			{
				String field = arg.substring(0, equalsIdx);
				String val = arg.substring(1 + equalsIdx);
				if(field.equalsIgnoreCase("aln_fn"))
				{
					pafFn = val;
				}
				if(field.equalsIgnoreCase("fasta_fn"))
				{
					fastaFn = val;
				}
				if(field.equalsIgnoreCase("read_fn"))
				{
					readFn = val;
				}
				if(field.equalsIgnoreCase("graph_fn"))
				{
					graphFn = val;
				}
				if(field.equalsIgnoreCase("read_map_file"))
				{
					readMapFile = val;
				}
				if(field.equalsIgnoreCase("outputbroken"))
				{
					outputBroken = true;
					brokenOutputFile = val;
				}
				if(field.equalsIgnoreCase("contig_map_file"))
				{
					contigMapFile = val;
				}
				if(field.equalsIgnoreCase("out_file"))
				{
					outFn = val;
				}
				if(field.equalsIgnoreCase("max_hanging"))
				{
					maxHanging = Integer.parseInt(val);
				}
				if(field.equalsIgnoreCase("minq"))
				{
					minQuality = Integer.parseInt(val);
				}
			}
		}
		if(pafFn.length() == 0 || fastaFn.length() == 0 || readFn.length() == 0)
		{
			usage();
			System.exit(1);
		}
		if(readMapFile.length() == 0 || contigMapFile.length() == 0 || outFn.length() == 0)
		{
			usage();
			System.exit(1);
		}
	}
	
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
		System.out.println("  read_map_file   (String) - TODO");
		System.out.println("  contig_map_file (String) - TODO");
		System.out.println("  out_file        (String) - the name of the file to output the scaffolded contigs to");
		System.out.println();
		System.out.println("Optional args");
		System.out.println("  max_hanging (int) [1000]    - the maximum amount by which the end of a contig can exceed the alignment and still be joined");
		System.out.println("  minq        (int) [40]      - the minimum quality of alignments needed to be kept");
		System.out.println("  graph_fn    (String) [none] - a GFA file containing an assembly graph, causing only alignments which are validated by the graph to be kept");
		System.out.println("  --break                     - allows original contigs to be broken");
		System.out.println();
	}
	
public static void main(String[] args) throws IOException
{
	parseArgs(args);
	
	Scanner input = new Scanner(new FileInputStream(new File(pafFn)));
	PrintWriter out = new PrintWriter(new File(outFn));
	
	/*
	 * Read in alignments and bucket by which read was aligned
	 */
	HashMap<String, ArrayList<SortablePafAlignment>> alignmentsPerRead = new HashMap<>();
	while(input.hasNext())
	{
		String line = input.nextLine();
		SortablePafAlignment cur = new SortablePafAlignment(line);
		
		double curThreshold = Math.min(.2 * cur.readLength, minAlignmentLength);
		
		// Filter out short alignments
		if(cur.readEnd - cur.readStart < curThreshold)
		{
			continue;
		}
		
		// Filter out low-quality alignments
		if(cur.mapq < minQuality)
		{
			continue;
		}
		
		String readName = cur.readName;
		
		ReadUtils.addInit(alignmentsPerRead, readName, cur);
	}
	
	ArrayList<CorrectMisassemblies.NovelAdjacency> corrections = new ArrayList<CorrectMisassemblies.NovelAdjacency>();
	if(allowBreaks)
	{
		corrections = CorrectMisassemblies.findMisassemblies(alignmentsPerRead);
		if(verbose)
		{
			for(CorrectMisassemblies.NovelAdjacency na : corrections)
			{
				System.err.println(na);
			}
			System.err.println("Number of misassemblies: " + corrections.size());
		}
	}
	
	HashSet<String> contigNames = new HashSet<>();
	
	CorrectMisassemblies.ContigBreaker splitter = new CorrectMisassemblies.ContigBreaker(corrections, contigNames);
	
	System.err.println("Number of breaks: " + splitter.numBreaks);
	
	alignmentsPerRead = CorrectMisassemblies.remapAll(splitter, alignmentsPerRead);
	
	/*
	 * Get chains of unique mappings and keep track of contigs/reads involved in them
	 */
	HashMap<String, ArrayList<ArrayList<SortablePafAlignment>>> chainsPerRead = new HashMap<>();
	HashSet<String> readNames = new HashSet<String>();
	for(String s : alignmentsPerRead.keySet())
	{
		if(alignmentsPerRead.get(s).size() == 1)
		{
			continue;
		}
		
		ArrayList<ArrayList<SortablePafAlignment>> chains = getUniqueMatches(alignmentsPerRead.get(s));
		
		if(chains.size() > 0)
		{
			for(ArrayList<SortablePafAlignment> l : chains)
			{
				for(SortablePafAlignment spa : l)
				{
					contigNames.add(spa.contigName);
				}
			}
		}
		
		if(chains.size() > 0)
		{
			chainsPerRead.put(s, chains);
			readNames.add(s);
		}
	}
	
	/*
	 * Output broken assembly
	 */
	if(outputBroken && corrections.size() > 0)
	{
		System.err.println("Outputting broken assembly");
		splitter.outputBrokenAssembly(fastaFn, brokenOutputFile);
	}
	
	/*
	 * Get sequences of relevant contigs/reads for merging
	 */
	HashMap<String, String> readMap = new HashMap<>(), contigMap = new HashMap<>();
	if(!fileMap || (readMap = ReadUtils.readMap(readMapFile)).size() == 0)
	{
		System.err.println("Filtering reads");
		if(readFn.endsWith(".fa") || readFn.endsWith(".fasta"))
		{
			readMap = ReadUtils.getFastaMap(readFn, readNames);
		}
		else
		{
			readMap = ReadUtils.getFastqMap(readFn, readNames);
		}
		ReadUtils.writeMap(readMapFile, readMap);
	}
	if(!fileMap || (contigMap = ReadUtils.readMap(contigMapFile)).size() == 0)
	{
		System.err.println("Filtering contigs");
		contigMap = ReadUtils.getFastaMap(fastaFn, contigNames);
		ReadUtils.writeMap(contigMapFile, contigMap);
	}
	
	if(verbose)
	{
		System.err.println("Split contigs:\n" +splitter.subcontigMap.keySet());
	}
	ArrayList<String> keys = new ArrayList<String>();
	keys.addAll(contigMap.keySet());
	for(String s : keys)
	{
		if(splitter.breakSequence(s, contigMap.get(s)))
		{
			contigMap.remove(s);
		}
	}
	
	for(String splitContigName : splitter.sequenceMap.keySet())
	{
		contigMap.put(splitContigName, splitter.sequenceMap.get(splitContigName));
	}
	
	System.err.println("Initializing frequency map for contig kmers");
	FreqqyMap freq = new FreqqyMap();
	
	System.err.println("Adding contig kmer frequencies");
	for(String s : contigMap.keySet())
	{
		freq.addKmerCount(s, contigMap.get(s));
	}
	
	System.err.println("Indexing contig kmer frequencies");
	for(String s : contigMap.keySet())
	{
		freq.addSumArray(s, contigMap.get(s));
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
	
	contigMap.put("undosplit", "A");
	for(String s : splitter.subcontigMap.keySet())
	{
		ArrayList<CorrectMisassemblies.ContigBreaker.Subcontig> subs = splitter.subcontigMap.get(s);
		int numSubcontigs = subs.size();
		for(int i = 0; i<numSubcontigs-1; i++)
		{
			sg.addEdge(subs.get(i).name, subs.get(i+1).name, "undosplit", 0, 0, 0, false, true, 1);
		}
	}
	
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
		String headerLine = createHeaderLine(scaffoldContigs.get(s), splitter);
		if(verbose)
		{
			System.err.println(headerLine);
		}
		out.println(headerLine);
		
		String seq = merge(scaffoldContigs.get(s), scaffoldEdges.get(s), readMap, contigMap);
		
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
	
	if(printOrientations)
	{
		printOrientations(scaffoldEdges, splitter);
	}
	
}

static void printOrientations(HashMap<String, ArrayDeque<ScaffoldGraph.Alignment>> als, CorrectMisassemblies.ContigBreaker splitter) throws IOException
{
	PrintWriter out = new PrintWriter(new File("orientations.txt"));
	for(String s : als.keySet())
	{
		ArrayDeque<ScaffoldGraph.Alignment> cur = als.get(s);
		String contigName = cur.peekFirst().from;
		String oldName = splitter.sourceMap.containsKey(contigName) ? splitter.sourceMap.get(contigName) : contigName;
		out.print(oldName + " " + (cur.peekFirst().myContigPrefix ? '-' : '+'));
		for(ScaffoldGraph.Alignment spa : cur)
		{
			contigName = spa.to;
			oldName = splitter.sourceMap.containsKey(contigName) ? splitter.sourceMap.get(contigName) : contigName;
			out.print(" " + oldName + " " + (cur.peekFirst().theirContigPrefix ? '+' : '-'));
		}
		out.println();
	}
	out.close();
}

/*
 * Create a Fasta header line for a scaffold based on the names of contigs which make it up
 * format is >contigs1&contig2&... contig1 contig2 contig3 ...
 */
static String createHeaderLine(ArrayDeque<String> contigs, CorrectMisassemblies.ContigBreaker splitter)
{
	StringBuilder res = new StringBuilder("");
	HashSet<String> contigSet = new HashSet<String>();
	for(String s : contigs)
	{
		if(res.length() > 0) res.append("&");
		else res.append(">");
		res.append(s);
		contigSet.add(s);
	}
	
	for(String s : contigs)
	{
		if(splitter.sourceMap.containsKey(s))
		{
			res.append(" " + splitter.sourceMap.get(s));
		}
		else
		{
			res.append(" " + s);
		}
	}
	return res.toString();
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
			if(verbose)
			{
				System.err.println(contigs.peekFirst()+" "+curSeq.length());
			}
			if(spa.myContigPrefix) curSeq = ReadUtils.reverseComplement(curSeq);
			res.append(curSeq);
		}
		if(verbose)
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
				
		if(verbose)
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
static void addEdges(ScaffoldGraph sg, ArrayList<SortablePafAlignment> als, FreqqyMap freq)
{
	SortablePafAlignment last = null;
	boolean lastReversed = false;
	for(int i = 0; i<als.size(); i++)
	{
		SortablePafAlignment spa = als.get(i);
		
		boolean curReversed = false;
		
		boolean[] cse = contigStartEnd(spa);
		boolean[] rse = readStartEnd(spa);
		
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
				if(weight >= minWeight)
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

/*
 * General end checking for alignments
 */
static boolean[] startEnd(int startPos, int endPos, int length)
{
	double curMaxHanging = Math.min(maxHangingProportion*length, maxHanging);
	return new boolean[] {startPos < curMaxHanging, endPos + curMaxHanging >= length};
}

/*
 * Whether or not an alignment contains the start/end of the contig involved
 */
static boolean[] contigStartEnd(SortablePafAlignment pa)
{
	return startEnd(pa.contigStart, pa.contigEnd, pa.contigLength);
}

/*
 * Whether or not an alignment contains the start/end of the read involved
 */
static boolean[] readStartEnd(SortablePafAlignment pa)
{
	return startEnd(pa.readStart, pa.readEnd, pa.readLength);
}

/*
 * Compresses the alignments to a given read by combining alignments of the same contig into one
 * Also, filters out invalid alignments
 */
static int MAX_GAP = 10000;
static ArrayList<SortablePafAlignment> compress(ArrayList<SortablePafAlignment> alignments, boolean filterInvalid)
{
	int n = alignments.size();

	Comparator<SortablePafAlignment> byContigName = new Comparator<SortablePafAlignment>() {

		@Override
		public int compare(SortablePafAlignment a, SortablePafAlignment b) {
			if(a.contigName.equals(b.contigName))
			{
				return a.compareTo(b);
			}
			return a.contigName.compareTo(b.contigName);
		}
	};
	
	// Sort by contig name and break ties by read start position
	Collections.sort(alignments, byContigName);
	ArrayList<SortablePafAlignment> filtered = new ArrayList<>();
	for(int i = 0; i<n; i++)
	{
		// Find the end of the run of alignments of the current contig
		int j = i+1;
		while(j< n && alignments.get(i).contigName.equals(alignments.get(j).contigName))
		{
			j++;
		}
		
		// Now alignments[i:j) has all the alignments of this contig - combine or remove them
		boolean[] rse = new boolean[2], cse = new boolean[2];
		boolean gapFree = true;
		int lastReadEndPosition = alignments.get(i).readEnd;
		int lastContigEndPosition = alignments.get(i).contigEnd;
		for(int k = i; k<j; k++)
		{
			SortablePafAlignment cur = alignments.get(k);
			
			if(k > i && alignments.get(k-1).strand != alignments.get(k).strand)
			{
				j = k;
				break;
			}
			
			// Check for a gap between this alignment and the last one in either the read or contig
			if(cur.contigStart - lastContigEndPosition > MAX_GAP)
			{
				gapFree = false;
				break;
			}
			if(cur.readStart - lastReadEndPosition > MAX_GAP)
			{
				gapFree = false;
				break;
			}
			
			lastContigEndPosition = cur.contigEnd;
			lastReadEndPosition = cur.readEnd;
			boolean[] crse = readStartEnd(cur);
			boolean[] ccse = contigStartEnd(cur);
			for(int q = 0; q<2; q++)
			{
				rse[q] |= crse[q];
				cse[q] |= ccse[q];
			}
		}
		
		/*
		 * If the set of alignments had a large gap, ignore it 
		 */
		if(!gapFree)
		{
			i = j - 1;
			continue;
		}
		
		/*
		 * We have whether the alignment set covers the start/end of contig/read, so check that it's valid
		 */
		if(filterInvalid && !cse[0] && !cse[1])
		{
			/*
			 * Middle portion of contig aligns somewhere on read but neither end of it
			 * Only possible case is read contained in contig, making alignment useless
			 * Throw out the alignment and update i
			 */
			i = j - 1;
			continue;
		}
		else if(filterInvalid && !rse[0] && !rse[1])
		{
			/*
			 * Neither end of the read is involved, so contig must be contained in the read
			 * Filter out cases which don't reflect this
			 */
			if(!cse[0] || !cse[1])
			{
				i = j - 1;
				continue;
			}
		}
		else
		{
			/*
			 * Create consensus of all of the alignments by taking earliest start and latest end
			 */
			SortablePafAlignment total = alignments.get(i).copy();
			for(int k = i+1; k<j; k++)
			{
				SortablePafAlignment cur = alignments.get(k);
				total.contigStart = Math.min(total.contigStart, cur.contigStart);
				total.contigEnd = Math.max(total.contigEnd, cur.contigEnd);
				total.readStart = Math.min(total.readStart, cur.readStart);
				total.readEnd = Math.max(total.readEnd, cur.readEnd);
			}
			i = j - 1;
			filtered.add(total);
		}
	}
	
	Collections.sort(filtered);
	return filtered;
}

/*
 * Gets chains of unique matches to a read given the list of all of the alignments to it
 */
@SuppressWarnings("unchecked")
static ArrayList<ArrayList<SortablePafAlignment>> getUniqueMatches(ArrayList<SortablePafAlignment> alignments)
{
	/*
	 * Sort by start point
	 */
	Collections.sort(alignments);
	
	/*
	 * Compress all alignments of the same contig and remove invalid alignments
	 */
	alignments = compress(alignments, true);
	
	/*
	 * List of chains of alignments
	 */
	ArrayList<ArrayList<SortablePafAlignment>> res = new ArrayList<>();
	
	/*
	 * The list of alignments in the current chain
	 */
	ArrayList<SortablePafAlignment> cur = new ArrayList<>();
	for(int i = 0 ; i<alignments.size(); i++)
	{
		
		SortablePafAlignment a = alignments.get(i);
		
		/*
		 * Cases: 
		 *   1.) Contained in a previous alignment -> Ignore this alignment
		 *   2.) Overlaps last two alignments -> End chain here
		 *   3.) Valid continuation of chain
		 */
		if(cur.size() >= 1 && cur.get(cur.size()-1).readEnd >= a.readEnd)
		{
			// Contained in a previous alignment
			continue;
		}
		else if(cur.size() >= 2 && cur.get(cur.size() - 2).readEnd > a.readStart)
		{
			// Overlaps last two alignments
			cur.remove(cur.size()-1);
			if(cur.size() >= 2)
			{
				res.add((ArrayList<SortablePafAlignment>) cur.clone());
			}
			cur.clear();
		}
		else if(cur.size() >= 1 && cur.get(cur.size() - 1).readEnd + 100000 < a.readStart)
		{
			if(cur.size() >= 2)
			{
				res.add((ArrayList<SortablePafAlignment>) cur.clone());
			}
			cur.clear();
			cur.add(a);
		}
		else
		{
			// Valid continuation of chain
			cur.add(a);
		}
	}
	
	// Add leftover chain
	if(cur.size() >= 2)
	{
		res.add(cur);
	}
	
	return res;
}

/*
 * Alignment of a contig to an ultralong read - sortable by start position in the read 
 */
static class SortablePafAlignment implements Comparable<SortablePafAlignment>
{
	String readName, contigName;
	int readLength, readStart, readEnd;
	int contigLength, contigStart, contigEnd;
	int mapq;
	char strand;
	String line;
	// Call with a second parameter to denote that read and contig were flipped when calling minimap2
	SortablePafAlignment(String line, int backwards)
	{
		this.line = line;
		String[] ss = line.split("\t");
		contigName = ss[0];
		contigLength = Integer.parseInt(ss[1]);
		contigStart = Integer.parseInt(ss[2]);
		contigEnd = Integer.parseInt(ss[3]);
		strand = ss[4].charAt(0);
		readName = ss[5];
		readLength = Integer.parseInt(ss[6]);
		readStart = Integer.parseInt(ss[7]);
		readEnd = Integer.parseInt(ss[8]);
		mapq = Integer.parseInt(ss[11]);
	}
	SortablePafAlignment(String line)
	{
		this.line = line;
		String[] ss = line.split("\t");
		readName = ss[0];
		readLength = Integer.parseInt(ss[1]);
		readStart = Integer.parseInt(ss[2]);
		readEnd = Integer.parseInt(ss[3]);
		strand = ss[4].charAt(0);
		contigName = ss[5];
		contigLength = Integer.parseInt(ss[6]);
		contigStart = Integer.parseInt(ss[7]);
		contigEnd = Integer.parseInt(ss[8]);
		mapq = Integer.parseInt(ss[11]);
	}

	public int compareTo(SortablePafAlignment o) {
		if(readStart != o.readStart)
		{
			return readStart - o.readStart;
		}
		return readEnd - o.readEnd;
	}
	SortablePafAlignment copy()
	{
		SortablePafAlignment res = new SortablePafAlignment(line);
		if(!res.contigName.equals(contigName)) res.contigName = contigName;
		if(res.contigStart != contigStart) res.contigStart = contigStart;
		if(res.contigLength != contigLength) res.contigLength = contigLength;
		if(res.contigEnd != contigEnd) res.contigEnd = contigEnd;
		return res;
	}
	
}
}