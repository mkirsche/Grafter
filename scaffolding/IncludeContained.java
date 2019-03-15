package scaffolding;

import java.util.*;

import java.io.*;

public class IncludeContained {
	static int maxHanging = 100;
@SuppressWarnings("resource")
public static void main(String[] args) throws IOException
{
	String pafFn = "rel2_200kplus_ccs_mat.paf";
	String fastaFn = "maternal_and_unknown.contigs.mmpoa.fa";
	String readFn = "rel2_200kplus.fastq";
	String readMapFile = "readmap_maternal.txt";
	String contigMapFile = "contigmap_maternal.txt";
	String outFn = "new_contigs.fa";
	
	if(args.length > 0 && args[0].equals("--server"))
	{
		pafFn = "/scratch/groups/mschatz1/mkirsche/ultralong/ccs/rel2_200kplus_ccs.paf";
		fastaFn = "/scratch/groups/mschatz1/mkirsche/ultralong/ccs/paternal_and_unknown.contigs.mmpoa.fa";
		readFn = "/scratch/groups/mschatz1/mkirsche/ultralong/rel2_200kplus.fastq";
		readMapFile = "/scratch/groups/mschatz1/mkirsche/ultralong/ccs/readmap_paternal2.txt";
		contigMapFile = "/scratch/groups/mschatz1/mkirsche/ultralong/ccs/contigmap_paternal2.txt";
		outFn = "/scratch/groups/mschatz1/mkirsche/ultralong/ccs/paternal_newcontigs2.fa";
	}
	
	else if(args.length >= 6)
	{
		pafFn = args[0];
		fastaFn = args[1];
		readFn = args[2];
		readMapFile = args[3];
		contigMapFile = args[4];
		outFn = args[5];
	}
	
	Scanner input = new Scanner(new FileInputStream(new File(pafFn)));
	HashMap<String, ArrayList<SortablePafAlignment>> alignmentsPerRead = new HashMap<>();
	PrintWriter out = new PrintWriter(new File(outFn));
	
	/*
	 * Read in alignments and bucket by which read was aligned
	 */
	while(input.hasNext())
	{
		String line = input.nextLine();
		SortablePafAlignment cur = new SortablePafAlignment(line);
		
		// Filter out short alignments
		if(cur.readEnd - cur.readStart < 10000)
		{
			continue;
		}
		String readName = cur.readName;
		if(!alignmentsPerRead.containsKey(readName))
		{
			alignmentsPerRead.put(readName, new ArrayList<SortablePafAlignment>());
		}
		alignmentsPerRead.get(readName).add(cur);
	}
	
	/*
	 * Get chains of unique mappings and list of relevant contigs
	 */
	HashMap<String, ArrayList<ArrayList<SortablePafAlignment>>> uniqueMap = new HashMap<>();
	HashSet<String> contigNames = new HashSet<String>();
	HashSet<String> readNames = new HashSet<String>();
	HashSet<String> alreadyJoined = new HashSet<String>();
	for(String s : alignmentsPerRead.keySet())
	{
		if(alignmentsPerRead.get(s).size() == 1)
		{
			continue;
		}
		
		ArrayList<ArrayList<SortablePafAlignment>> uniques = getUniqueMatches(alignmentsPerRead.get(s), alreadyJoined);
		
		if(uniques.size() > 0)
		{
			for(ArrayList<SortablePafAlignment> l : uniques)
			{
				for(SortablePafAlignment spa : l)
				{
					contigNames.add(spa.contigName);
				}
			}
		}
		
		if(uniques.size() > 0)
		{
			uniqueMap.put(s, uniques);
			readNames.add(s);
		}
	}
	
	
	/*
	 * Get sequences of relevant contigs/reads for merging
	 */
	HashMap<String, String> readMap, contigMap;
	if((readMap = Scaffold.readMap(readMapFile)).size() == 0)
	{
		System.err.println("Filtering reads");
		readMap = Scaffold.getFastqMap(readFn, readNames);
		Scaffold.writeMap(readMapFile, readMap);
	}
	if((contigMap = Scaffold.readMap(contigMapFile)).size() == 0)
	{
		System.err.println("Filtering contigs");
		contigMap = Scaffold.getFastaMap(fastaFn, contigNames);
		Scaffold.writeMap(contigMapFile, contigMap);
	}
	
	/*
	 * Perform string merging
	 */
	ScaffoldGraph sg = new ScaffoldGraph();
	System.err.println("Joining contigs");
	int numMerged = 0;
	HashMap<String, String> mergedContigs = new HashMap<String, String>();
	HashMap<String, HashSet<String>> components = new HashMap<>();
	for(String readName : uniqueMap.keySet())
	{
		//String readSeq = readMap.get(readName);
		ArrayList<ArrayList<SortablePafAlignment>> allAlignments = uniqueMap.get(readName);
		for(ArrayList<SortablePafAlignment> als : allAlignments)
		{
			addEdges(sg, als);
			//String seq = merge(als, readSeq, contigMap);
			String scaffoldName = "";
			for(int i = 0; i<als.size(); i++)
			{
				scaffoldName += als.get(i).contigName;
				if(i < als.size() - 1) scaffoldName += "&";
			}
			//numMerged += als.size() - 1;
			//mergedContigs.put(scaffoldName, seq);
			if(!components.containsKey(scaffoldName))
			{
				components.put(scaffoldName, new HashSet<String>());
			}
			for(int i = 0; i<als.size(); i++)
			{
				components.get(scaffoldName).add(als.get(i).contigName);
			}
		}
	}
	
	// Map from first contig in a scaffold to info about the scaffold
	HashMap<String, ArrayDeque<String>> scaffoldContigs = new HashMap<>();
	HashMap<String, String> lastToFirst = new HashMap<String, String>();
	HashMap<String, ArrayDeque<ScaffoldGraph.Alignment>> scaffoldEdges = new HashMap<>();
	HashSet<String> usedContigs = new HashSet<String>();
	
	for(String s : sg.adj.keySet())
	{
		// Ignore if already at start or middle of a scaffold
		if(usedContigs.contains(s) && !lastToFirst.containsKey(s))
		{
			continue;
		}
		ScaffoldGraph.Alignment best = consensus(s, sg.adj.get(s), usedContigs, scaffoldEdges);
		if(best == null) continue;
		String to = best.to;
		
		// TODO update 4 variables above
		if(!usedContigs.contains(to))
		{
			if(usedContigs.contains(s))
			{
				// s is the last contig in some scaffold
				String firstContigInScaffold = lastToFirst.get(s);
				ArrayDeque<String> allContigsInScaffold = scaffoldContigs.get(firstContigInScaffold);
				ArrayDeque<ScaffoldGraph.Alignment> allEdgesInScaffold = scaffoldEdges.get(firstContigInScaffold);
				allContigsInScaffold.addLast(to);
				allEdgesInScaffold.add(best);
				
				usedContigs.add(to);
				lastToFirst.remove(s);
				lastToFirst.put(to, firstContigInScaffold);
			}
			else
			{
				// s hasn't been connected to anything yet
				scaffoldEdges.put(s, new ArrayDeque<ScaffoldGraph.Alignment>());
				scaffoldContigs.put(s, new ArrayDeque<String>());
				scaffoldContigs.get(s).addLast(s);
				scaffoldContigs.get(s).addLast(to);
				scaffoldEdges.get(s).addLast(best);
				
				usedContigs.add(s);
				usedContigs.add(to);
				lastToFirst.put(to, s);
			}
		}
		
		else
		{
			// In calculating best, already made sure it's the first in its scaffold
			// Move entire scaffold with to to the end of the scaffold with s
			String lastContigInScaffold = scaffoldContigs.get(to).peekLast();
			if(usedContigs.contains(s))
			{
				String firstContigInScaffold = lastToFirst.get(s);
				lastToFirst.remove(s);
				lastToFirst.put(lastContigInScaffold, firstContigInScaffold);
				scaffoldEdges.get(firstContigInScaffold).addLast(best);
				while(!scaffoldEdges.get(to).isEmpty())
				{
					scaffoldEdges.get(firstContigInScaffold).addLast(scaffoldEdges.get(to).pollFirst());
				}
				scaffoldEdges.remove(to);
				
				while(!scaffoldContigs.get(to).isEmpty())
				{
					scaffoldContigs.get(firstContigInScaffold).addLast(scaffoldContigs.get(to).pollFirst());
				}
				scaffoldContigs.remove(to);
			}
			else
			{
				// s is on its own
				scaffoldEdges.put(s, scaffoldEdges.get(to));
				scaffoldEdges.get(s).addFirst(best);
				scaffoldEdges.remove(to);
				
				scaffoldContigs.put(s, scaffoldContigs.get(to));
				scaffoldContigs.get(s).addFirst(s);
				scaffoldContigs.remove(to);
				
				usedContigs.add(s);
				
				lastToFirst.put(lastContigInScaffold, s);
			}
		}
		
		numMerged++;
	}
	
	/*
	 * Output all merged contigs
	 */
/*	for(String s : mergedContigs.keySet())
	{
		out.print(">" + s);
		for(String sub : components.get(s)) out.print(" " + sub);
		out.println();
		out.println(mergedContigs.get(s));
	}
	*/
	for(String s : scaffoldContigs.keySet())
	{
		out.println(getHeaderLine(scaffoldContigs.get(s)));
		
		String seq = merge(scaffoldContigs.get(s), scaffoldEdges.get(s), readMap, contigMap);
		
		out.println(seq);
	}
	System.err.println("Number of joins: " + numMerged);
	
}

static String getHeaderLine(ArrayDeque<String> contigs)
{
	StringBuilder res = new StringBuilder("");
	for(String s : contigs)
	{
		if(res.length() > 0) res.append("&");
		else res.append(">");
		res.append(s);
	}
	
	for(String s : contigs) res.append(" " + s);
	
	return res.toString();
}

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
			if(spa.myContigPrefix) curSeq = Scaffold.reverseComplement(curSeq);
			res.append(curSeq);
		}
		int overlap = 0;
		if(spa.myReadEnd < spa.theirReadStart)
		{
			String readSeq = readMap.get(spa.read);
			res.append(readSeq.substring(spa.myReadEnd, spa.theirReadStart));
		}
		else
		{
			overlap = spa.myReadEnd - spa.theirReadStart;
		}
		
		String curSeq = relevantContigs.get(spa.to);
		
		if(!spa.theirContigPrefix)
		{
			curSeq = Scaffold.reverseComplement(curSeq);
		}
		res.append(curSeq.substring(overlap));
	}
	return res.toString();
}

/*
 * Gets the best alignment of a contig to follow a given contig
 */
static ScaffoldGraph.Alignment consensus(String from, ArrayList<ScaffoldGraph.Alignment> als, HashSet<String> usedContigs, HashMap<String, ArrayDeque<ScaffoldGraph.Alignment>> scaffoldEdges)
{
	HashMap<String, ArrayList<ScaffoldGraph.Alignment>> prefEdges = new HashMap<>();
	HashMap<String, ArrayList<ScaffoldGraph.Alignment>> suffEdges = new HashMap<>();
	for(ScaffoldGraph.Alignment a : als)
	{
		if(a.myContigPrefix)
		{
			if(!prefEdges.containsKey(a.to)) prefEdges.put(a.to, new ArrayList<ScaffoldGraph.Alignment>());
			prefEdges.get(a.to).add(a);
		}
		else
		{
			if(!suffEdges.containsKey(a.to)) suffEdges.put(a.to, new ArrayList<ScaffoldGraph.Alignment>());
			suffEdges.get(a.to).add(a);
		}
	}
	ArrayList<ScaffoldGraph.Alignment> best = null;
	for(String to : prefEdges.keySet())
	{
		if(usedContigs.contains(to) && !scaffoldEdges.containsKey(to)) continue;
		if(scaffoldEdges.containsKey(to) && scaffoldEdges.get(to).peekLast().to.equals(from)) continue;
		ArrayList<ScaffoldGraph.Alignment> al = prefEdges.get(to);
		ArrayList<ScaffoldGraph.Alignment> valid = new ArrayList<>();
		for(ScaffoldGraph.Alignment a : al)
		{
			if(scaffoldEdges.containsKey(to) && scaffoldEdges.get(to).peekFirst().myContigPrefix == a.theirContigPrefix)
			{
				continue;
			}
			valid.add(a);
		}
		
		if(best == null || valid.size() > best.size())
		{
			best = valid;
		}
	}
	for(String to : suffEdges.keySet())
	{
		if(usedContigs.contains(to) && !scaffoldEdges.containsKey(to)) continue;
		if(scaffoldEdges.containsKey(to) && scaffoldEdges.get(to).peekLast().to.equals(from)) continue;
		ArrayList<ScaffoldGraph.Alignment> valid = new ArrayList<>();
		ArrayList<ScaffoldGraph.Alignment> al = suffEdges.get(to);
		for(ScaffoldGraph.Alignment a : al)
		{
			if(scaffoldEdges.containsKey(to) && scaffoldEdges.get(to).peekFirst().myContigPrefix == a.theirContigPrefix)
			{
				continue;
			}
			valid.add(a);
		}
		if(best == null || valid.size() > best.size())
		{
			best = valid;
		}
	}
	
	if(best == null || best.size() == 0) return null;
	
	ScaffoldGraph.Alignment res = best.get(0);
	
	return res;
}

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
		
		if(last != null)
		{
			System.out.println("Edge: " + last.contigName+" "+spa.contigName);
			int overlap = last.readEnd - spa.readStart;
			if(overlap <= spa.contigLength && overlap <= last.contigLength)
				sg.addEdge(last.contigName, spa.contigName, last.readName, last.readEnd, spa.readStart, lastReversed, !curReversed);
		}
		
		last = spa;
		lastReversed = curReversed;
	}
}

static boolean[] contigStartEnd(FindUsefulScaffoldingAlignments.PafAlignment pa)
{
	return new boolean[] {pa.contigStart < maxHanging, pa.contigEnd + maxHanging >= pa.contigLength};
}

static boolean[] readStartEnd(FindUsefulScaffoldingAlignments.PafAlignment pa)
{
	return new boolean[] {pa.readStart < maxHanging, pa.readEnd + maxHanging >= pa.readLength};
}

static ArrayList<SortablePafAlignment> compress(ArrayList<SortablePafAlignment> alignments)
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
	Collections.sort(alignments, byContigName);
	ArrayList<SortablePafAlignment> filtered = new ArrayList<>();
	for(int i = 0; i<n; i++)
	{
		int j = i+1;
		while(j< n && alignments.get(i).contigName.equals(alignments.get(j).contigName))
		{
			j++;
		}
		// Now alignments[i:j) has all the alignments of this contig - combine or remove them
		boolean[] rse = new boolean[2], cse = new boolean[2];
		boolean okay = true;
		int lrs = alignments.get(i).readStart, lcs = alignments.get(i).contigStart;
		for(int k = i; k<j; k++)
		{
			SortablePafAlignment cur = alignments.get(k);
			int contigDist = cur.contigStart - lcs;
			int readDist = cur.readStart - lrs;
			if(readDist > 3 * contigDist || contigDist > 3 * readDist) okay = false;
			lcs = cur.contigStart;
			lrs = cur.readStart;
			boolean[] crse = readStartEnd(cur);
			boolean[] ccse = contigStartEnd(cur);
			for(int q = 0; q<2; q++)
			{
				rse[q] |= crse[q];
				cse[q] |= ccse[q];
			}
		}
		if(!okay)
		{
			i = j - 1;
			continue;
		}
		// Have whether the alignment set covers the start/end of contig/read
		if(!cse[0] && !cse[1])
		{
			/*
			 * Middle portion of contig aligns somewhere on read but neither end of it
			 * Only possible case is read contained in contig, making alignment useless
			 * Throw out the alignment and update i
			 */
			i = j-1;
			continue;
		}
		else if(!rse[0] && !rse[1])
		{
			/*
			 * Neither end of the contig is involved, so it must be contained in the read
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
			SortablePafAlignment total = alignments.get(i).copy();
			for(int k = i+1; k<j; k++)
			{
				SortablePafAlignment cur = alignments.get(k);
				total.contigStart = Math.min(total.contigStart, cur.contigStart);
				total.contigEnd = Math.max(total.contigEnd, cur.contigEnd);
				total.readStart = Math.min(total.readStart, cur.readStart);
				total.readEnd = Math.max(total.readEnd, cur.readEnd);
			}
			i = j-1;
			filtered.add(total);
		}
	}
	
	Collections.sort(filtered);
	return filtered;
}

@SuppressWarnings("unchecked")
static ArrayList<ArrayList<SortablePafAlignment>> getUniqueMatches(ArrayList<SortablePafAlignment> alignments, HashSet<String> alreadyJoined)
{
	Collections.sort(alignments);
	alignments = compress(alignments);
	ArrayList<ArrayList<SortablePafAlignment>> res = new ArrayList<>();
	
	ArrayList<SortablePafAlignment> cur = new ArrayList<>();
	HashSet<String> using = new HashSet<String>();
	boolean found = false;
	for(int i = 0 ; i<alignments.size(); i++)
	{
		
		SortablePafAlignment a = alignments.get(i);
		
		// Types: 0 is addition to alignment chain, 1 is invalid (overlapping last two alignments),
		// 2 is contained/not contig end so ignore this alignment
		int type = -1;
		boolean[] contigStartEnd = contigStartEnd(a);
		boolean[] readStartEnd = readStartEnd(a);
		/*if(alreadyJoined.contains(a.contigName) || using.contains(a.contigName))
		{
			vals.add(0);
			type = 2;
		}*/
		/*if(!contigStartEnd[0] && !contigStartEnd[1])
		{
			vals.add(1);
			type = 2;
		}
		else if((!contigStartEnd[1] || !contigStartEnd[0]) && (!readStartEnd[0] && !readStartEnd[1]))
		{
			vals.add(2);
			// Middle of read aligning to one end of contig -> invalid
			type = 2;
			cur.clear();
			using.clear();
		}*/
		if(cur.size() >= 2 && cur.get(cur.size() - 2).readEnd > a.readStart)
		{
			type = 1;
		}
		else if(cur.size() >= 1 && cur.get(cur.size()-1).readEnd >= a.readEnd)
		{
			type = 2;
		}
		else
		{
			type = 0;
		}
		
		if(type == 0)
		{
			using.add(a.contigName);
			cur.add(a);
		}
		else if(type == 1)
		{
			// Ambiguous whether two contigs ago goes with this one or the last one
			cur.remove(cur.size()-1);
			if(cur.size() >= 2)
			{
				res.add((ArrayList<SortablePafAlignment>) cur.clone());
				for(SortablePafAlignment spa : cur)
				{
					alreadyJoined.add(spa.contigName);
				}
			}
			cur.clear();
			using.clear();
		}
		else if(type == 2)
		{
			continue;
		}
	}
	
//	if(found)
//	{
//		System.out.println("cs: " + cur.size());
//		System.out.println(types+" "+vals);
//	}
	
	if(cur.size() >= 2)
	{
		res.add(cur);
		for(SortablePafAlignment spa : cur)
		{
			alreadyJoined.add(spa.contigName);
		}
	}
	
	if(found) System.out.println(res.size());
	
	return res;
}
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
		return new SortablePafAlignment(line);
	}
	
}
}
