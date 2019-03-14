package scaffolding;

import java.util.*;

import java.io.*;

public class IncludeContained {
	static int maxHanging = 100;
@SuppressWarnings("resource")
public static void main(String[] args) throws IOException
{
	String pafFn = "rel2_200kplus_ccs_mat.paf";
	String fastaFn = "";
	String readFn = "";
	String readMapFile = "readmap_paternal.txt";
	String contigMapFile = "contigmap_paternal.txt";
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
			//System.err.print("Chain sizes for " + s+":");
			for(ArrayList<SortablePafAlignment> l : uniques)
			{
				for(SortablePafAlignment spa : l)
				{
					contigNames.add(spa.contigName);
				}
				//System.err.print(" "+l.size());
			}
			//System.err.println();
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
	System.err.println("Joining contigs");
	int numMerged = 0;
	HashMap<String, String> mergedContigs = new HashMap<String, String>();
	for(String readName : uniqueMap.keySet())
	{
		String readSeq = readMap.get(readName);
		ArrayList<ArrayList<SortablePafAlignment>> allAlignments = uniqueMap.get(readName);
		for(ArrayList<SortablePafAlignment> als : allAlignments)
		{
			String seq = merge(als, readSeq, contigMap);
			String scaffoldName = "";
			for(int i = 0; i<als.size(); i++)
			{
				scaffoldName += als.get(i).contigName;
				if(i < als.size() - 1) scaffoldName += "&";
			}
			numMerged += als.size() - 1;
			mergedContigs.put(scaffoldName, seq);
		}
	}
	
	/*
	 * Output all merged contigs
	 */
	for(String s : mergedContigs.keySet())
	{
		out.println(">" + s);
		out.println(mergedContigs.get(s));
	}
	System.err.println("Number of joins: " + numMerged);
	
}

static String merge(ArrayList<SortablePafAlignment> als, String readSeq, HashMap<String, String> relevantContigs)
{
	StringBuilder res = new StringBuilder();
	SortablePafAlignment last = null;
	for(int i = 0; i<als.size(); i++)
	{
		SortablePafAlignment spa = als.get(i);
		int overlap = 0;
		if(last != null)
		{
			if(last.readEnd < spa.readStart)
			{
				res.append(readSeq.substring(last.readEnd, spa.readStart));
			}
			else
			{
				overlap = last.readEnd - spa.readStart;
			}
		}
		
		String curSeq = relevantContigs.get(spa.contigName);
		
		boolean[] cse = contigStartEnd(spa);
		boolean[] rse =readStartEnd(spa);
		
		if(cse[0] && cse[1])
		{
			// Entire contig aligns - have to look at strand in alignment
			if(spa.strand == '-')
			{
				curSeq = Scaffold.reverseComplement(curSeq);
			}
		}
		else if(cse[0])
		{
			// Beginning of contig - we would expect suffix of read if same strand
			if(rse[0] && rse[1])
			{
				if(spa.strand == '-')
				{
					curSeq = Scaffold.reverseComplement(curSeq);
				}
			}
			else if(rse[0])
			{
				curSeq = Scaffold.reverseComplement(curSeq);
			}
		}
		else
		{
			// End of contig - we would expect prefix of read if same strand
			if(rse[0] && rse[1])
			{
				if(spa.strand == '-')
				{
					curSeq = Scaffold.reverseComplement(curSeq);
				}
			}
			else if(rse[1])
			{
				curSeq = Scaffold.reverseComplement(curSeq);
			}
		}
		
		res.append(curSeq.substring(overlap));
	}
	return res.toString();
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
		for(int k = i; k<j; k++)
		{
			SortablePafAlignment cur = alignments.get(k);
			boolean[] crse = readStartEnd(cur);
			boolean[] ccse = contigStartEnd(cur);
			for(int q = 0; q<2; q++)
			{
				rse[q] |= crse[q];
				cse[q] |= ccse[q];
			}
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
			filtered.add(total);
		}
	}
	
	Collections.sort(filtered);
	return filtered;
	/*
	Collections.sort(alignments);
	ArrayList<SortablePafAlignment> res = new ArrayList<SortablePafAlignment>();
	SortablePafAlignment last = alignments.get(0);
	for(int i = 1; i<n; i++)
	{
		SortablePafAlignment cur = alignments.get(i);
		if(cur.contigName.equals(last.contigName))
		{
			last.contigEnd = cur.contigEnd;
			last.readEnd = cur.readEnd;
		}
		else
		{
			res.add(last);
			last = cur;
		}
	}
	res.add(last);
	return res;
	*/
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
	ArrayList<Integer> types = new ArrayList<Integer>(), vals = new ArrayList<Integer>();
	for(int i = 0 ; i<alignments.size(); i++)
	{
		
		SortablePafAlignment a = alignments.get(i);
		//System.out.println(a.contigName);
		if(a.contigName.contains("tig00084613"))
		{
			found = true;
			System.out.println("go");
			for(SortablePafAlignment aa : alignments) System.out.println(aa.contigName+" "+aa.readStart+" "+aa.readEnd + "("+aa.contigStart+", " + aa.contigEnd+")");
		}
		
		// Types: 0 is addition to alignment chain, 1 is invalid (overlapping last two alignments),
		// 2 is contained/not contig end so ignore this alignment
		int type = -1;
		boolean[] contigStartEnd = contigStartEnd(a);
		boolean[] readStartEnd = readStartEnd(a);
		if(alreadyJoined.contains(a.contigName) || using.contains(a.contigName))
		{
			vals.add(0);
			type = 2;
		}
		else if(!contigStartEnd[0] && !contigStartEnd[1])
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
		}
		else if(cur.size() >= 2 && cur.get(cur.size() - 2).readEnd > a.readStart)
		{
			vals.add(3);
			type = 1;
		}
		else if(cur.size() >= 1 && cur.get(cur.size()-1).readEnd >= a.readEnd)
		{
			vals.add(4);
			type = 2;
		}
		else
		{
			vals.add(5);
			type = 0;
		}
		
		types.add(type);
		
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
	
	if(found)
	{
		System.out.println("cs: " + cur.size());
		System.out.println(types+" "+vals);
	}
	
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
