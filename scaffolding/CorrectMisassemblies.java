/*
 * Code for correcting misassemblies
 */

package scaffolding;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

import scaffolding.IncludeContained.SortablePafAlignment;

/*
 * Module for correcting misassemblies as part of the ultralong read scaffolder
 */
public class CorrectMisassemblies {
	
	// How far away misassemblies must be from the end of a contig to be worth breaking
	static int buffer = 25000;
	
	// The minimum ratio of evidence for vs. against a misassembly to believe it
	static double evidenceRatio = 2; // Value is 1.5 for full data
	
	// The maximum total evidence which can go against a misassembly and still have it be considered
	static double maxEvidence = 450000 / evidenceRatio;
	
	// The minimum weight of an alignment for it to be considered a piece of evidence towards a misassembly
	static double minSingleAlignmentWeight = 20000;
	
	// The number of reads which must support a misassembly
	static int minInversionSupport = 1;
	static int minChimeraSupport = 2;
	static int minSplitSupport = 4; // Higher because only requires one endpoint
	
/*
 * Finds inversions based on alignments of contigs to ultralong reads
 * An inversion is defined as the alignments from a contig changing strand
 */
static ArrayList<NovelAdjacency> findInversions(ArrayList<IncludeContained.SortablePafAlignment> alignments)
{
	ArrayList<NovelAdjacency> res = new ArrayList<>();
	
	int n = alignments.size();

	// Group by contig name, and sort each group by contig start position
	Comparator<SortablePafAlignment> byContigName = new Comparator<SortablePafAlignment>() {

		@Override
		public int compare(SortablePafAlignment a, SortablePafAlignment b) {
			if(a.contigName.equals(b.contigName))
			{
				return a.contigStart - b.contigStart;
			}
			return a.contigName.compareTo(b.contigName);
		}
	};
	Collections.sort(alignments, byContigName);
	
	for(int i = 0; i<n; i++)
	{
		// Find the end of the run of alignments of the current contig
		int j = i+1;
		while(j < n && alignments.get(i).contigName.equals(alignments.get(j).contigName))
		{
			j++;
		}
		
		// Now alignments[i:j) has all the alignments of this contig - look for changes in strand
		for(int k = i+1; k<j; k++)
		{
			SortablePafAlignment last = alignments.get(k-1), cur = alignments.get(k);
			
			// Make sure the alignments are close together on the contig
			if(last.strand != cur.strand && cur.contigStart < last.contigEnd + buffer && cur.readStart < last.readEnd + buffer)
			{
				// Weight the misassembly by the harmonic mean of alignment lengths - discard if weight is too small
				double weight = harmonicMean(last.contigEnd - last.contigStart, cur.contigEnd - cur.contigStart);
				if(weight >= minSingleAlignmentWeight)
				{
					boolean lastPrefix = last.strand == '-';
					boolean curPrefix = cur.strand == '+';
					res.add(new NovelAdjacency(last.contigName, cur.contigName, lastPrefix ? last.contigStart : last.contigEnd, 
						curPrefix ? cur.contigStart : cur.contigEnd, last.contigLength, cur.contigLength, last.readName, weight, 0));
				}
			}
		}
	}
	
	return res;
}
/*
 * Takes all alignments to a read and looks for evidence of chimeric contigs
 * This is where the middle of one contig should align to a different contig rather than the rest of its given contig
 */
static ArrayList<NovelAdjacency> findChimeras(ArrayList<IncludeContained.SortablePafAlignment> alignments)
{
	ArrayList<NovelAdjacency> res = new ArrayList<NovelAdjacency>();
	
	// Combine alignments of the same contig to a single read, but do not filter out invalid ones
	ArrayList<IncludeContained.SortablePafAlignment> compressed = IncludeContained.compress(alignments, false);
	if(compressed.size() < 2) return res;
	IncludeContained.SortablePafAlignment last = compressed.get(0);
	for(int i = 1; i < compressed.size(); i++)
	{
		IncludeContained.SortablePafAlignment cur = compressed.get(i);
		
		boolean lastPrefix = last.strand == '-';
		boolean curPrefix = cur.strand == '+';
		
		// For each contig involved, see if there is an attempt to join a position in the middle of it
		boolean lastNonEnd = false, curNonEnd = false;
		
		if(lastPrefix && last.contigStart > buffer)
		{
			lastNonEnd = true;
		}
		else if(!lastPrefix && last.contigEnd + buffer < last.contigLength)
		{
			lastNonEnd = true;
		}
		
		if(curPrefix && cur.contigStart > buffer)
		{
			curNonEnd = true;
		}
		
		else if(!curPrefix && cur.contigEnd + buffer < cur.contigLength)
		{
			curNonEnd = true;
		}
		
		// Make sure that there is a chimera and that the alignments don't overlap
		if((curNonEnd || lastNonEnd) && cur.readStart <= last.readEnd + 1000)
		{
			res.add(new NovelAdjacency(last.contigName, cur.contigName, lastPrefix ? last.contigStart : last.contigEnd, 
					curPrefix ? cur.contigStart : cur.contigEnd, last.contigLength, cur.contigLength, last.readName, 
					harmonicMean(last.contigEnd - last.contigStart, cur.contigEnd - cur.contigStart), 1));
		}
		
		last = cur;
	}
	return res;
}
/*
 * Find split alignments - this is where a lot of reads have alignments starting or ending at the same location of a contig
 * In this case, it may make sense to break the contig at that position
 * Note that strict thresholds are used here because the alignments are noisy and easily interrupted by repeats
 */
static ArrayList<NovelAdjacency> findSplitAlignments(HashMap<String, ArrayList<IncludeContained.SortablePafAlignment>> byContig)
{
	int maxEndpointDist = 100;
	ArrayList<NovelAdjacency> res = new ArrayList<>();
	
	for(String contigName : byContig.keySet())
	{
		TreeMap<Integer, Double> endpointWeights = new TreeMap<>();
		TreeMap<Integer, Integer> endpointFrequency = new TreeMap<>();
		TreeMap<Integer, String> readSupport = new TreeMap<Integer, String>();
		ArrayList<IncludeContained.SortablePafAlignment> als = byContig.get(contigName);
		for(IncludeContained.SortablePafAlignment spa : als)
		{
			double curWeight = spa.contigEnd - spa.contigStart;
			int[] ends = new int[] {spa.contigStart, spa.contigEnd};
			int contigLength = spa.contigLength;
			for(int endpoint : ends)
			{
				if(endpoint < buffer*2 || endpoint + buffer*2 > contigLength)
				{
					continue;
				}
				Integer floor = endpointWeights.floorKey(endpoint);
				Integer ceiling = endpointWeights.ceilingKey(endpoint);
				int floorDist = floor == null ? (int)1e9 : endpoint - floor;
				int ceilingDist = ceiling == null ? (int)1e9 : ceiling - endpoint;
				if(floorDist <= ceilingDist && floorDist < maxEndpointDist)
				{
					// Consider this the same as the previous endpoint
					endpointWeights.put(floor, endpointWeights.get(floor) + curWeight);
					endpointFrequency.put(floor, 1 + endpointFrequency.get(floor));
				}
				else if(ceilingDist < maxEndpointDist)
				{
					// Consider this the same as the next endpoint
					endpointWeights.put(ceiling, endpointWeights.get(ceiling) + curWeight);
					endpointFrequency.put(ceiling, 1 + endpointFrequency.get(ceiling));
				}
				else
				{
					endpointWeights.put(endpoint, curWeight);
					endpointFrequency.put(endpoint, 1);
					readSupport.put(endpoint, spa.readName);
				}
			}
		}
		for(int endpoint : endpointWeights.keySet())
		{
			if(endpointFrequency.get(endpoint) < 5)
			{
				continue;
			}
			if(endpointWeights.get(endpoint) < 100000)
			{
				continue;
			}
			NovelAdjacency toAdd =(new NovelAdjacency(contigName, contigName, endpoint, endpoint, 
					als.get(0).contigLength, als.get(0).contigLength, readSupport.get(endpoint), 
					endpointWeights.get(endpoint), 2));
			toAdd.support = endpointFrequency.get(endpoint);
			res.add(toAdd);
		}
	}
	
	return res;
}
static boolean check(NovelAdjacency na, HashMap<String, ArrayList<IncludeContained.SortablePafAlignment>> byContig)
{
	double evidence = 0.0;
	ArrayList<IncludeContained.SortablePafAlignment> contig1aln = byContig.get(na.contig1);
	for(IncludeContained.SortablePafAlignment spa : contig1aln)
	{
		if(spa.contigStart < na.pos1 - 10000 && spa.contigEnd > na.pos1 + 10000)
		{
			evidence += harmonicMean(na.pos1 - spa.contigStart, spa.contigEnd - na.pos1);
		}
	}
	ArrayList<IncludeContained.SortablePafAlignment> contig2aln = byContig.get(na.contig2);
	for(IncludeContained.SortablePafAlignment spa : contig2aln)
	{
		if(spa.contigStart < na.pos2 && spa.contigEnd > na.pos2)
		{
			evidence += harmonicMean(na.pos2 - spa.contigStart, spa.contigEnd - na.pos2);
		}
	}
	System.out.println(na.contig1+" "+na.contig2+" "+na.weight+" "+na.pos1+" "+na.pos2+" "+evidence);
	return evidence < maxEvidence && (evidence * evidenceRatio < na.weight || (na.contig1.equals(na.contig2) && evidence * evidenceRatio < na.weight));
}

/*
 * Given a list of novel adjacencies, combine those between the same contigs which are at very similar positions
 * Also, filter out those which have a lot of alignments spanning their supposed split points 
 */
static ArrayList<NovelAdjacency> compressAndFilter(ArrayList<NovelAdjacency> nas, boolean filter, HashMap<String, ArrayList<IncludeContained.SortablePafAlignment>> byContig)
{
	ArrayList<NovelAdjacency> res = new ArrayList<NovelAdjacency>();
	for(int i = 0; i<nas.size(); i++)
	{
		NovelAdjacency cur = nas.get(i);
		int j = i+1;
		int totSupport = nas.get(i).support;
		while(j < nas.size())
		{
			NovelAdjacency next = nas.get(j);
			if(!next.contig1.equals(cur.contig1)) break;
			if(!next.contig2.equals(cur.contig2)) break;
			if(Math.abs(next.pos1 - cur.pos1) > 10000) break;
			if(Math.abs(next.pos2 - cur.pos2) > 10000) break;
			cur.weight += next.weight;
			totSupport += next.support;
			cur.pos1 = (int)(cur.pos1 * (j-i) + next.pos1) / (j - i + 1);
			cur.pos2 = (int)(cur.pos2 * (j-i) + next.pos2) / (j - i + 1);
			j++;
		}
		boolean hasSupport = (cur.type == 0 && totSupport >= minInversionSupport) ||
				(cur.type == 1 && totSupport >= minChimeraSupport) || cur.type == 2;
		if(!filter || (cur.weight > 20000 && j >= i+3 && check(cur, byContig)))
		{
			res.add(cur);
		}
		else if(hasSupport && check(cur, byContig))
		{
			res.add(cur);
		}
		i = j-1;
	}
	return res;
}
static ArrayList<NovelAdjacency> findMisassemblies(HashMap<String, ArrayList<IncludeContained.SortablePafAlignment>> alignmentsPerRead)
{
	HashMap<String, ArrayList<IncludeContained.SortablePafAlignment>> byContig = reindex(alignmentsPerRead);
	ArrayList<NovelAdjacency> corrections = new ArrayList<CorrectMisassemblies.NovelAdjacency>();
	for(String s : alignmentsPerRead.keySet())
	{
		ArrayList<CorrectMisassemblies.NovelAdjacency> tmp = 
				CorrectMisassemblies.findChimeras(alignmentsPerRead.get(s));
		
		ArrayList<NovelAdjacency> inv = findInversions(alignmentsPerRead.get(s));
		tmp.addAll(inv);
		
		if(tmp.size() > 0)
		{
			for(CorrectMisassemblies.NovelAdjacency na : tmp)
			{
				corrections.add(na);
			}
		}
	}
	ArrayList<NovelAdjacency> splitAlignments = findSplitAlignments(byContig);
	corrections.addAll(splitAlignments);
	Collections.sort(corrections);
	
	return compressAndFilter(corrections, true, byContig);
}
static double harmonicMean(double x, double y)
{
	return 2.0 * x * y / (x+y); 
}
/*
 * Take a set of alignments grouped by read and instead group them by contig
 */
static HashMap<String, ArrayList<IncludeContained.SortablePafAlignment>> reindex(HashMap<String, ArrayList<IncludeContained.SortablePafAlignment>> alignments)
{
	HashMap<String, ArrayList<IncludeContained.SortablePafAlignment>> res = new HashMap<>();
	
	for(String s : alignments.keySet())
	{
		ArrayList<SortablePafAlignment> curList = alignments.get(s);
		for(SortablePafAlignment spa : curList)
		{
			String newKey = spa.contigName;
			ReadUtils.addInit(res, newKey, spa);
		}
	}
	return res;
}
/*
 * A novel adjacency (or misassembly) is a point or pair of points where an existing assembly should be broken
 */
static class NovelAdjacency implements Comparable<NovelAdjacency>
{
	// The contigs involved in the novel adjacency
	String contig1, contig2;
	
	// The position in each contig where the misassembly is present
	int pos1, pos2;
	
	// A read whose alignments support the misassembly
	String read;
	
	// The weight of the misassembly giving some measure of the amount of evidence supportingits presence
	double weight;
	
	// 0 is inversion, 1 is chimera, 2 is split
	int type;
	
	// Whether or not the prefix of each contig is involved in the misassembly
	int length1, length2;
	int support;
	NovelAdjacency(String c1, String c2, int p1, int p2, int l1, int l2, String rr, double ww, int tt)
	{
		support = 1;
		contig1 = c1;
		contig2 = c2;
		pos1 = p1;
		pos2 = p2;
		length1 = l1;
		length2 = l2;
		read = rr;
		weight = ww;
		type = tt;
		if(contig1.compareTo(contig2) > 0 || (contig1.equals(contig2) && pos1 > pos2))
		{
			String tmp = contig1;
			contig1 = contig2;
			contig2 = tmp;
			int tmppos = pos1;
			pos1 = pos2;
			pos2 = tmppos;
			int tmplen = length1;
			length1 = length2;
			length2 = tmplen;
		}
	}
	public String toString()
	{
		return "Novel adjacency: " + contig1 + " " + pos1 + " " + length1 + " " 
				+ contig2 + " " + pos2 + " " + length2 + " " + read + " "
				+ weight + " " + type;
	}
	@Override
	public int compareTo(NovelAdjacency o) {
		if(!contig1.equals(o.contig1)) return contig1.compareTo(o.contig1);
		if(!contig2.equals(o.contig2)) return contig2.compareTo(o.contig2);
		return pos1 - o.pos1;
	}
}
static HashMap<String, ArrayList<IncludeContained.SortablePafAlignment>> remapAll(ContigBreaker splitter, HashMap<String, ArrayList<IncludeContained.SortablePafAlignment>> alignmentsPerRead)
{
	HashMap<String, ArrayList<IncludeContained.SortablePafAlignment>> res = new HashMap<>();
	for(String readName : alignmentsPerRead.keySet())
	{
		for(IncludeContained.SortablePafAlignment spa : alignmentsPerRead.get(readName))
		{
			IncludeContained.SortablePafAlignment cur = splitter.remap(spa);
			if(cur == null)
			{
				continue;
			}
			ReadUtils.addInit(res, readName, cur);
		}
	}
	return res;
}
/*
 * Handling information for breaking contigs in light of misassemblies
 */
static class ContigBreaker
{
	HashMap<String, ArrayList<Subcontig>> subcontigMap;
	HashMap<String, ArrayList<Integer>> breakPositionMap;
	HashMap<String, Integer> lengthMap;
	HashMap<String, String> sequenceMap;
	HashMap<String, String> sourceMap;
	int numBreaks;
	ContigBreaker(ArrayList<NovelAdjacency> nas, HashSet<String> contigNames)
	{
		subcontigMap = new HashMap<>();
		breakPositionMap = new HashMap<>();
		lengthMap = new HashMap<>();
		sequenceMap = new HashMap<>();
		sourceMap = new HashMap<>();
		for(NovelAdjacency na : nas)
		{
			if(na.pos1 > buffer && na.pos1 + buffer < na.length1)
			{
				ReadUtils.addInit(breakPositionMap, na.contig1, na.pos1);
				lengthMap.put(na.contig1, na.length1);
			}
			if(na.contig2.equals(na.contig1) && na.pos1 == na.pos2)
			{
				continue;
			}
			if(na.pos2 > buffer && na.pos2 + buffer < na.length2)
			{
				ReadUtils.addInit(breakPositionMap, na.contig2, na.pos2);
				lengthMap.put(na.contig2, na.length2);
			}
		}
		numBreaks = 0;
		for(String contigName : breakPositionMap.keySet())
		{
			contigNames.add(contigName);
			ArrayList<Subcontig> scs = destroy(contigName, breakPositionMap.get(contigName));
			numBreaks += scs.size() - 1;
			for(Subcontig sc : scs)
			{
				sourceMap.put(sc.name, contigName);
			}
			subcontigMap.put(contigName, scs);
		}
	}
	boolean breakSequence(String contigName, String sequence)
	{
		if(subcontigMap.containsKey(contigName))
		{
			ArrayList<Subcontig> scs = subcontigMap.get(contigName);
			for(Subcontig sc : scs)
			{
				sequenceMap.put(sc.name, sequence.substring(sc.startPos, sc.endPos));
			}
			return true;
		}
		return false;
	}
	IncludeContained.SortablePafAlignment remap(IncludeContained.SortablePafAlignment old)
	{
		IncludeContained.SortablePafAlignment res = old.copy();
		if(!subcontigMap.containsKey(old.contigName))
		{
			return res;
		}
		ArrayList<Subcontig> scs = subcontigMap.get(old.contigName);
		for(int i = 0; i<scs.size(); i++)
		{
			Subcontig sc = scs.get(i);
			boolean afterStart = old.contigStart > sc.startPos - buffer;
			boolean beforeEnd = old.contigEnd < sc.endPos + buffer;
			if(afterStart && beforeEnd)
			{
				res.contigName = sc.name;
				res.contigLength = sc.endPos - sc.startPos;
				res.contigStart = Math.min(res.contigLength, Math.max(1, old.contigStart - sc.startPos));
				res.contigEnd = Math.min(sc.endPos, old.contigEnd) - sc.startPos;
				return res;
			}
		}
		return null;
	}
	static ArrayList<Integer> filterBreaks(ArrayList<Integer> breaks, int length)
	{
		ArrayList<Integer> res = new ArrayList<>();
		for(int i = 0; i < breaks.size(); i++)
		{
			int cur = breaks.get(i);
			if(cur > length)
			{
				cur = length;
			}
			if(i == 0 || cur > breaks.get(i-1) + buffer)
			{
				res.add(cur);
			}
		}
		return res;
	}
	ArrayList<Subcontig> destroy(String contigName, ArrayList<Integer> breakPositions)
	{
		ArrayList<Subcontig> res = new ArrayList<>();
		Collections.sort(breakPositions);
		breakPositions = filterBreaks(breakPositions, lengthMap.get(contigName));
		res.add(new Subcontig(0, breakPositions.get(0), contigName + "_" + 1, contigName));
		
		for(int i = 0; i<breakPositions.size(); i++)
		{
			int startPos = breakPositions.get(i);
			int endPos = (i == breakPositions.size() - 1) ? lengthMap.get(contigName) : breakPositions.get(i+1);
			res.add(new Subcontig(startPos, endPos, contigName + "_" + (i+2), contigName));
		}
		return res;
	}
	@SuppressWarnings("resource")
	void outputBrokenAssembly(String fn, String ofn) throws IOException
	{
		BufferedReader br = new BufferedReader(new FileReader(new File(fn)));
		PrintWriter out = new PrintWriter(new File(ofn));
		
		String contigName = br.readLine().split(" ")[0].substring(1);
		StringBuilder seq = new StringBuilder("");
		
		while(true)
		{
			try {
				String line = br.readLine();
				if(line.charAt(0) == '>')
				{
					// process last read
					print(contigName, seq.toString(), out);
					seq = new StringBuilder("");
					// new read name
					contigName = line.split(" ")[0].substring(1);
				}
				else
				{
					seq.append(line);
				}
			} catch(Exception e) {
				break;
			}
		}
		print(contigName, seq.toString(), out);
		out.close();
	}
	void print(String contigName, String seq, PrintWriter out)
	{
		if(breakSequence(contigName, seq))
		{
			ArrayList<Subcontig> scs = subcontigMap.get(contigName);
			for(Subcontig sc : scs)
			{
				out.println(">"+sc.name+"\n"+sequenceMap.get(sc.name));
			}
		}
		else
		{
			out.println(">"+contigName+"\n"+seq);
		}
		
	}
	@SuppressWarnings("resource")
	static HashMap<String, String> getFastaMap(String fn, HashSet<String> names) throws IOException
	{
		HashMap<String, String> res = new HashMap<String, String>();
		BufferedReader br = new BufferedReader(new FileReader(new File(fn)));
		String readName = br.readLine().split(" ")[0].substring(1);
		StringBuilder seq = new StringBuilder("");
		boolean useful = names.contains(readName);
		while(true)
		{
			try {
				String line = br.readLine();
				if(line.charAt(0) == '>')
				{
					// process last read
					if(useful)
					{
						res.put(readName, seq.toString());
						seq = new StringBuilder("");
					}
					// new read name
					readName = line.split(" ")[0].substring(1);
					useful = names.contains(readName);
				}
				else if(useful)
				{
					seq.append(line);
				}
			} catch(Exception e) {
				break;
			}
		}
		if(useful)
		{
			res.put(readName, seq.toString());
		}
		return res;
	}
	static class Subcontig
	{
		int startPos, endPos;
		String name;
		String oldName;
		Subcontig(int ss, int ee, String nn, String oo)
		{
			startPos = ss;
			endPos = ee;
			name = nn;
			oldName = oo;
		}
	}
}
}