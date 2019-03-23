package scaffolding;

import java.util.*;

import scaffolding.CorrectMisassemblies.ContigBreaker;

import java.io.*;

public class IncludeContained {
	
	// The number of reads required to support the joining of two contigs
	static int minReadSupport = 1;
	static double minWeightSupport = 50000;
	
	static double maxHangingProportion = 0.02;
	static int maxHanging = 100;
	static boolean fileMap = false;
	static boolean verbose = true;
	static boolean outputBroken = false;
	static boolean allowBreaks = true;
@SuppressWarnings("resource")
public static void main(String[] args) throws IOException
{
	/*
	 * File names for testing locally
	 */
	String pafFn = "rel2_200kplus_ccs_mat.paf";
	//String pafFn = "brokenmaternal.paf";
	String fastaFn = "maternal_and_unknown.contigs.mmpoa.fa";
	//String fastaFn = "maternal_and_unknown.contigs.mmpoa.intra.chimera.broken.fa";
	String readFn = "rel2_200kplus.fastq";
	String readMapFile = "readmap_maternal.txt";
	String contigMapFile = "contigmap_maternal.txt";
	String outFn = "new_contigs.fa";
	String brokenOutputFile = fastaFn;
	
	/*String pafFn = "AssemblyToMaternal.paf";
	String fastaFn = "maternal_and_unknown.contigs.mmpoa.fa";
	String readFn = "rel2wt_50kplus.ctg.fa";
	String readMapFile = "readmap_maternal.txt";
	String contigMapFile = "contigmap_maternal.txt";
	String outFn = "new_contigs.fa";*/
	
	/*
	 * Default files for testing on the server
	 */
	if(args.length > 0 && args[0].equals("--server"))
	{
		pafFn = "/scratch/groups/mschatz1/mkirsche/ultralong/ccs/rel2_200kplus_ccs.paf";
		fastaFn = "/scratch/groups/mschatz1/mkirsche/ultralong/ccs/paternal_and_unknown.contigs.mmpoa.fa";
		readFn = "/scratch/groups/mschatz1/mkirsche/ultralong/rel2_200kplus.fastq";
		readMapFile = "/scratch/groups/mschatz1/mkirsche/ultralong/ccs/readmap_paternal2.txt";
		contigMapFile = "/scratch/groups/mschatz1/mkirsche/ultralong/ccs/contigmap_paternal2.txt";
		outFn = "/scratch/groups/mschatz1/mkirsche/ultralong/ccs/paternal_newcontigs2.fa";
	}
	
	/*
	 * File names passed as arguments
	 */
	else if(args.length >= 6)
	{
		pafFn = args[0];
		fastaFn = args[1];
		readFn = args[2];
		readMapFile = args[3];
		contigMapFile = args[4];
		outFn = args[5];
		for(int i = 6; i<args.length; i++)
		{
			if(args[i].equals("--break"))
			{
				allowBreaks = true;
			}
			else if(args[i].startsWith("--outputbroken="))
			{
				outputBroken = true;
				brokenOutputFile = args[i].substring(args[i].indexOf('=')+1);
			}
		}
	}
	
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
		
		// Filter out short alignments
		if(cur.readEnd - cur.readStart < 5000)
		{
			continue;
		}
		
		String readName = cur.readName;
		
		addInit(alignmentsPerRead, readName, cur);
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
	if(!fileMap || (readMap = Scaffold.readMap(readMapFile)).size() == 0)
	{
		System.err.println("Filtering reads");
		if(readFn.endsWith(".fa") || readFn.endsWith(".fasta"))
		{
			readMap = Scaffold.getFastaMap(readFn, readNames);
		}
		else
		{
			readMap = Scaffold.getFastqMap(readFn, readNames);
		}
		Scaffold.writeMap(readMapFile, readMap);
	}
	if(!fileMap || (contigMap = Scaffold.readMap(contigMapFile)).size() == 0)
	{
		System.err.println("Filtering contigs");
		contigMap = Scaffold.getFastaMap(fastaFn, contigNames);
		Scaffold.writeMap(contigMapFile, contigMap);
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
			addEdges(sg, chain);
		}
	}
	
	// Map from first contig in a scaffold to info about the scaffold
	HashMap<String, ArrayDeque<String>> scaffoldContigs = new HashMap<>();
	HashMap<String, String> lastToFirst = new HashMap<String, String>();
	HashMap<String, ArrayDeque<ScaffoldGraph.Alignment>> scaffoldEdges = new HashMap<>();
	HashSet<String> usedContigs = new HashSet<String>();
	
	/*
	 * Create scaffolds by linking together contigs based on their best edges
	 */
	for(String s : sg.adj.keySet())
	{
		for(int strand = 0; strand < 2; strand++)
		{
			// Ignore if already at start or middle of a scaffold
			if(usedContigs.contains(s) && !lastToFirst.containsKey(s))
			{
				continue;
			}
			
			// Get the consensus edge of all edges going to the most highly supported contig
			ScaffoldGraph.Alignment best = consensus(s, sg.adj.get(s)[strand], usedContigs, scaffoldEdges, lastToFirst);
			if(best == null) continue;
			if(lastToFirst.containsKey(s) && scaffoldEdges.get(lastToFirst.get(s)).peekLast().theirContigPrefix == best.myContigPrefix)
			{
				continue;
			}
			String t = best.to;
			
			if(!usedContigs.contains(t))
			{
				// t is by itself in a contig
				
				if(usedContigs.contains(s))
				{
					// s is the last contig in some scaffold
					String firstContigInScaffold = lastToFirst.get(s);
					ArrayDeque<String> allContigsInScaffold = scaffoldContigs.get(firstContigInScaffold);
					ArrayDeque<ScaffoldGraph.Alignment> allEdgesInScaffold = scaffoldEdges.get(firstContigInScaffold);
					allContigsInScaffold.addLast(t);
					allEdgesInScaffold.add(best);
					
					usedContigs.add(t);
					lastToFirst.remove(s);
					lastToFirst.put(t, firstContigInScaffold);
				}
				else
				{
					// s hasn't been connected to anything yet
					scaffoldEdges.put(s, new ArrayDeque<ScaffoldGraph.Alignment>());
					scaffoldContigs.put(s, new ArrayDeque<String>());
					scaffoldContigs.get(s).addLast(s);
					scaffoldContigs.get(s).addLast(t);
					scaffoldEdges.get(s).addLast(best);
					
					usedContigs.add(s);
					usedContigs.add(t);
					lastToFirst.put(t, s);
				}
			}
			
			else
			{
				// In calculating best, already made sure it's the first or last in its scaffold
				// Move entire scaffold with t to the end of the scaffold with s
				boolean tFirst = scaffoldContigs.containsKey(t);
				String lastContigInScaffold = tFirst ? scaffoldContigs.get(t).peekLast() : lastToFirst.get(t);
				String tScaffoldKey = tFirst ? t : lastContigInScaffold;
				if(usedContigs.contains(s))
				{
					// s is the last contig in a scaffold, so append the scaffold with to after s
					String firstContigInScaffold = lastToFirst.get(s);
					lastToFirst.remove(s);
					lastToFirst.put(lastContigInScaffold, firstContigInScaffold);
					scaffoldEdges.get(firstContigInScaffold).addLast(best);
					ArrayDeque<ScaffoldGraph.Alignment> tScaffoldEdges = scaffoldEdges.get(tScaffoldKey);
					if(tFirst)
					{
						while(!tScaffoldEdges.isEmpty())
						{
							ScaffoldGraph.Alignment cur = tScaffoldEdges.pollFirst();
							scaffoldEdges.get(firstContigInScaffold).addLast(cur);
						}
					}
					else
					{
						String lastTo = tScaffoldKey;
						ArrayDeque<ScaffoldGraph.Alignment> toAdd = new ArrayDeque<ScaffoldGraph.Alignment>();
						while(!tScaffoldEdges.isEmpty())
						{
							ScaffoldGraph.Alignment cur = tScaffoldEdges.pollFirst();
							toAdd.addLast(cur.reverse(lastTo));
							lastTo = cur.to;
						}
						while(!toAdd.isEmpty())
						{
							scaffoldEdges.get(firstContigInScaffold).addLast(toAdd.pollLast());
						}
					}
					
					scaffoldEdges.remove(tScaffoldKey);
					
					ArrayDeque<String> tScaffoldContigs = scaffoldContigs.get(tScaffoldKey);
					while(!tScaffoldContigs.isEmpty())
					{
						scaffoldContigs.get(firstContigInScaffold).addLast(tFirst ? tScaffoldContigs.pollFirst() : tScaffoldContigs.pollLast());
					}
					scaffoldContigs.remove(tScaffoldKey);
					if(!tFirst)
					{
						lastToFirst.remove(t);
					}
				}
				else
				{
					// s is on its own, so add it to the beginning of the scaffold with t
					if(tFirst)
					{
						scaffoldEdges.put(s, scaffoldEdges.get(t));
						scaffoldEdges.get(s).addFirst(best);
						scaffoldEdges.remove(t);
						scaffoldContigs.put(s, scaffoldContigs.get(t));
						scaffoldContigs.get(s).addFirst(s);
						scaffoldContigs.remove(t);
						lastToFirst.put(lastContigInScaffold, s);
					}
					else
					{
						// t is at the end of its scaffold, so reverse scaffold and add it to new scaffold
						
						// Deal with edges
						scaffoldEdges.put(s, new ArrayDeque<ScaffoldGraph.Alignment>());
						scaffoldEdges.get(s).addFirst(best);
						String lastTo = tScaffoldKey;
						ArrayDeque<ScaffoldGraph.Alignment> toAdd = new ArrayDeque<ScaffoldGraph.Alignment>();
						while(!scaffoldEdges.get(tScaffoldKey).isEmpty())
						{
							ScaffoldGraph.Alignment cur = scaffoldEdges.get(tScaffoldKey).pollFirst();
							toAdd.addLast(cur.reverse(lastTo));
							lastTo = cur.to;
						}
						while(!toAdd.isEmpty())
						{
							scaffoldEdges.get(s).addLast(toAdd.pollLast());
						}
						scaffoldEdges.remove(tScaffoldKey);
						
						// Deal with list of contigs
						scaffoldContigs.put(s, new ArrayDeque<String>());
						scaffoldContigs.get(s).addFirst(s);
						while(!scaffoldContigs.get(tScaffoldKey).isEmpty())
						{
							scaffoldContigs.get(s).addLast(scaffoldContigs.get(tScaffoldKey).pollLast());
						}
						scaffoldContigs.remove(tScaffoldKey);
						
						// Deal with lastToFirst
						lastToFirst.remove(t);
						lastToFirst.put(tScaffoldKey, s);
					}
					usedContigs.add(s);					
				}
			}
			
			numMerged++;
		}
	}
	
	for(String s : scaffoldContigs.keySet())
	{
		if(verbose)
		{
			System.err.println(getHeaderLine(scaffoldContigs.get(s), splitter));
		}
		out.println(getHeaderLine(scaffoldContigs.get(s), splitter));
		
		String seq = merge(scaffoldContigs.get(s), scaffoldEdges.get(s), readMap, contigMap);
		
		out.println(seq);
	}
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
	
}

/*
 * Add a key, value pair to a map from string to list, but initialize list if it's not already there
 */
static <T> void addInit(HashMap<String, ArrayList<T>> map, String key, T val)
{
	if(!map.containsKey(key)) map.put(key, new ArrayList<T>());
	map.get(key).add(val);
}

/*
 * Create a Fasta header line for a scaffold based on the names of contigs which make it up
 * format is >contigs1&contig2&... contig1 contig2 contig3 ...
 */
static String getHeaderLine(ArrayDeque<String> contigs, ContigBreaker splitter)
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
			if(spa.myContigPrefix) curSeq = Scaffold.reverseComplement(curSeq);
			res.append(curSeq);
		}
		if(verbose)
		{
			System.err.println(spa.to+" "+spa.myContigPrefix+" "+spa.theirContigPrefix+" "+spa.myReadEnd+" "+spa.theirReadStart+" "+spa.strand+" "+" "+spa.read);
		}
		int overlap = 0;
		if(spa.myReadEnd < spa.theirReadStart)
		{
			String readSeq = readMap.get(spa.read);
			if(!spa.theirContigPrefix)
			{
				readSeq = Scaffold.reverseComplement(readSeq);
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
			curSeq = Scaffold.reverseComplement(curSeq);
		}
		res.append(curSeq.substring(overlap));
	}
	return res.toString();
}

/*
 * Gets the best alignment of another contig to follow a given contig
 */
static ScaffoldGraph.Alignment consensus(String from, ArrayList<ScaffoldGraph.Alignment> als, HashSet<String> usedContigs, HashMap<String, ArrayDeque<ScaffoldGraph.Alignment>> scaffoldEdges, HashMap<String, String> lastToFirst)
{
	HashMap<String, ArrayList<ScaffoldGraph.Alignment>> prefEdges = new HashMap<>();
	HashMap<String, ArrayList<ScaffoldGraph.Alignment>> suffEdges = new HashMap<>();
	for(ScaffoldGraph.Alignment a : als)
	{
		if(a.myContigPrefix)
		{
			addInit(prefEdges, a.to, a);
		}
		else
		{
			addInit(suffEdges, a.to, a);
		}
	}
	ArrayList<ScaffoldGraph.Alignment> best = null;
	double bestTotalWeight = 0;
	for(String to : prefEdges.keySet())
	{
		// Make sure the destination is on one of the ends of its scaffold
		if(usedContigs.contains(to) && !scaffoldEdges.containsKey(to) && !lastToFirst.containsKey(to)) continue;
		
		// Make sure that the destination isn't the beginning of the scaffold the edge is coming from
		if(scaffoldEdges.containsKey(to) && scaffoldEdges.get(to).peekLast().to.equals(from)) continue;
		
		// Make sure that the destination isn't the end of the scaffold the edge is coming from
		if(lastToFirst.containsKey(to) && lastToFirst.get(to).equals(from)) continue;
		
		ArrayList<ScaffoldGraph.Alignment> al = prefEdges.get(to);
		ArrayList<ScaffoldGraph.Alignment> valid = new ArrayList<>();
		double totalWeight = 0;
		for(ScaffoldGraph.Alignment a : al)
		{
			// Make sure this edge doesn't use the same side of the destination as an existing edge
			if(scaffoldEdges.containsKey(to) && scaffoldEdges.get(to).peekFirst().myContigPrefix == a.theirContigPrefix)
			{
				continue;
			}
			if(lastToFirst.containsKey(to) && scaffoldEdges.get(lastToFirst.get(to)).peekLast().theirContigPrefix == a.theirContigPrefix)
			{
				continue;
			}
			
			totalWeight += a.weight;
			valid.add(a);
		}
		
		if(best == null || totalWeight > bestTotalWeight)
		{
			best = valid;
			bestTotalWeight = totalWeight;
		}
	}
	for(String to : suffEdges.keySet())
	{
		// Make sure the destination is on one of the ends of its scaffold
		if(usedContigs.contains(to) && !scaffoldEdges.containsKey(to) && !lastToFirst.containsKey(to)) continue;
				
		// Make sure that the destination isn't the beginning of the scaffold the edge is coming from
		if(scaffoldEdges.containsKey(to) && scaffoldEdges.get(to).peekLast().to.equals(from)) continue;
				
		// Make sure that the destination isn't the end of the scaffold the edge is coming from
		if(lastToFirst.containsKey(to) && lastToFirst.get(to).equals(from)) continue;
		
		ArrayList<ScaffoldGraph.Alignment> valid = new ArrayList<>();
		ArrayList<ScaffoldGraph.Alignment> al = suffEdges.get(to);
		double totalWeight = 0;
		for(ScaffoldGraph.Alignment a : al)
		{
			// Make sure this edge doesn't use the same side of the destination as an existing edge
			if(scaffoldEdges.containsKey(to) && scaffoldEdges.get(to).peekFirst().myContigPrefix == a.theirContigPrefix)
			{
				continue;
			}
			if(lastToFirst.containsKey(to) && scaffoldEdges.get(lastToFirst.get(to)).peekLast().theirContigPrefix == a.theirContigPrefix)
			{
				continue;
			}
			totalWeight += a.weight;
			valid.add(a);
		}
		if(best == null || totalWeight > bestTotalWeight)
		{
			best = valid;
			bestTotalWeight = totalWeight;
		}
	}
	
	if(best == null || best.size() == 0) return null;
	
	if(best.size() < minReadSupport || bestTotalWeight < minWeightSupport) return null;
	
	
	ScaffoldGraph.Alignment res = best.get(0);
	
	return res;
}

/*
 * Add edges to a scaffold graph based on a chain of alignments to the same read
 */
static void addEdges(ScaffoldGraph sg, ArrayList<SortablePafAlignment> als)
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
			if(overlap <= spa.contigLength && overlap <= last.contigLength)
			{
				double lastLength = last.contigEnd - last.contigStart;
				double curLength = spa.contigEnd - spa.contigStart;
				double weight = 2 * lastLength * curLength / (lastLength + curLength);
				sg.addEdge(last.contigName, spa.contigName, last.readName, last.readEnd, spa.readStart, spa.readLength, lastReversed, !curReversed, weight);
			}
		}
		
		last = spa;
		lastReversed = curReversed;
	}
}

/*
 * Whether or not an alignment contains the start/end of the cintig involved
 */
static boolean[] contigStartEnd(FindUsefulScaffoldingAlignments.PafAlignment pa)
{
	int length = pa.contigEnd - pa.contigStart;
	double curMaxHanging = Math.min(maxHangingProportion*length, maxHanging);
	return new boolean[] {pa.contigStart < curMaxHanging, pa.contigEnd + curMaxHanging >= pa.contigLength};
}

/*
 * Whether or not an alignment contains the start/end of the read involved
 */
static boolean[] readStartEnd(FindUsefulScaffoldingAlignments.PafAlignment pa)
{
	int length = pa.readEnd - pa.readStart;
	return new boolean[] {pa.readStart < maxHanging * length, pa.readEnd + maxHanging * length >= pa.readLength};
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
		boolean sameStrand = true;
		int lastReadEndPosition = alignments.get(i).readEnd;
		int lastContigEndPosition = alignments.get(i).contigEnd;
		for(int k = i; k<j; k++)
		{
			SortablePafAlignment cur = alignments.get(k);
			
			if(k > i && alignments.get(k-1).strand != alignments.get(k).strand)
			{
				sameStrand = false;
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
			 * Create of all of the alignments by taking earliest start and latest end
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
static class SortablePafAlignment extends FindUsefulScaffoldingAlignments.PafAlignment implements Comparable<SortablePafAlignment>
{
	SortablePafAlignment(String line) {
		super(line);
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
