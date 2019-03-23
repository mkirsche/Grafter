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

public class CorrectMisassemblies {
	
static ArrayList<NovelAdjacency> findInversions(ArrayList<IncludeContained.SortablePafAlignment> alignments)
{
	ArrayList<NovelAdjacency> res = new ArrayList<>();
	
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
	for(int i = 0; i<n; i++)
	{
		// Find the end of the run of alignments of the current contig
		int j = i+1;
		while(j< n && alignments.get(i).contigName.equals(alignments.get(j).contigName))
		{
			j++;
		}
		
		// Now alignments[i:j) has all the alignments of this contig - look for changes in strand
		for(int k = i+1; k<j; k++)
		{
			SortablePafAlignment last = alignments.get(k-1), cur = alignments.get(k);
			boolean lastPrefix = last.strand == '-';
			boolean curPrefix = cur.strand == '+';
			if(last.strand != cur.strand)
			{
				double weight = harmonicMean(last.contigEnd - last.contigStart, cur.contigEnd - cur.contigStart);
				if(weight >= 20000)
				{
					res.add(new NovelAdjacency(last.contigName, cur.contigName, lastPrefix ? last.contigStart : last.contigEnd, 
						curPrefix ? cur.contigStart : cur.contigEnd, last.contigLength, cur.contigLength, last.readName, 
								weight, lastPrefix, curPrefix));
				}
			}
		}
	}
	
	return res;
}
	static int buffer = 50000;	
/*
 * Takes all alignments to a read and looks for evidence of chimeric contigs
 */
static ArrayList<NovelAdjacency> findChimeras(ArrayList<IncludeContained.SortablePafAlignment> alignments)
{
	ArrayList<NovelAdjacency> res = new ArrayList<NovelAdjacency>();
	ArrayList<NovelAdjacency> inv = findInversions(alignments);
	res.addAll(inv);
	ArrayList<IncludeContained.SortablePafAlignment> compressed = IncludeContained.compress(alignments, false);
	if(compressed.size() < 2) return res;
	IncludeContained.SortablePafAlignment last = compressed.get(0);
	for(int i = 1; i < compressed.size(); i++)
	{
		IncludeContained.SortablePafAlignment cur = compressed.get(i);
		
		boolean lastPrefix = last.strand == '-';
		boolean curPrefix = cur.strand == '+';
		
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
		
		if(curNonEnd || lastNonEnd && cur.readStart <= last.readEnd + 1000)
		{
			res.add(new NovelAdjacency(last.contigName, cur.contigName, lastPrefix ? last.contigStart : last.contigEnd, 
					curPrefix ? cur.contigStart : cur.contigEnd, last.contigLength, cur.contigLength, last.readName, 
					harmonicMean(last.contigEnd - last.contigStart, cur.contigEnd - cur.contigStart),
					lastPrefix, curPrefix));
		}
		
		last = cur;
	}
	return res;
}
static ArrayList<NovelAdjacency> findMisassemblies(HashMap<String, ArrayList<IncludeContained.SortablePafAlignment>> alignmentsPerRead)
{
	HashMap<String, ArrayList<IncludeContained.SortablePafAlignment>> byContig = reindex(alignmentsPerRead);
	ArrayList<CorrectMisassemblies.NovelAdjacency> corrections = new ArrayList<CorrectMisassemblies.NovelAdjacency>();
	for(String s : alignmentsPerRead.keySet())
	{
		ArrayList<CorrectMisassemblies.NovelAdjacency> tmp = 
				CorrectMisassemblies.findChimeras(alignmentsPerRead.get(s));
		if(tmp.size() > 0)
		{
			for(CorrectMisassemblies.NovelAdjacency na : tmp)
			{
				corrections.add(na);
			}
		}
	}
	Collections.sort(corrections);
	
	return compressAndFilter(corrections, true, byContig);
}
static boolean check(NovelAdjacency na, HashMap<String, ArrayList<IncludeContained.SortablePafAlignment>> byContig)
{
	double evidence = 0.0;
	ArrayList<IncludeContained.SortablePafAlignment> contig1aln = byContig.get(na.contig1);
	for(IncludeContained.SortablePafAlignment spa : contig1aln)
	{
		if(spa.contigStart < na.pos1 && spa.contigEnd > na.pos1)
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
	//System.out.println(na.contig1+" "+na.contig2+" "+na.weight+" "+evidence);
	return evidence *  1.5 < na.weight;
}
static ArrayList<NovelAdjacency> compressAndFilter(ArrayList<NovelAdjacency> nas, boolean filter, HashMap<String, ArrayList<IncludeContained.SortablePafAlignment>> byContig)
{
	ArrayList<NovelAdjacency> res = new ArrayList<NovelAdjacency>();
	for(int i = 0; i<nas.size(); i++)
	{
		NovelAdjacency cur = nas.get(i);
		int j = i+1;
		while(j < nas.size())
		{
			NovelAdjacency next = nas.get(j);
			if(!next.contig1.equals(cur.contig1)) break;
			if(!next.contig2.equals(cur.contig2)) break;
			if(Math.abs(next.pos1 - cur.pos1) > 10000) break;
			if(Math.abs(next.pos2 - cur.pos2) > 10000) break;
			cur.weight += next.weight;
			cur.pos1 = (int)(cur.pos1 * (j-i) + next.pos1) / (j - i + 1);
			cur.pos2 = (int)(cur.pos2 * (j-i) + next.pos2) / (j - i + 1);
			j++;
		}
		if(!filter || (cur.weight > 20000 && j >= i+3 && check(cur, byContig)))
		{
			res.add(cur);
		}
		i = j-1;
	}
	return res;
}
static double harmonicMean(double x, double y)
{
	return 2.0 * x * y / (x+y); 
}
static HashMap<String, ArrayList<IncludeContained.SortablePafAlignment>> reindex(HashMap<String, ArrayList<IncludeContained.SortablePafAlignment>> alignments)
{
	HashMap<String, ArrayList<IncludeContained.SortablePafAlignment>> res = new HashMap<>();
	
	for(String s : alignments.keySet())
	{
		ArrayList<SortablePafAlignment> curList = alignments.get(s);
		for(SortablePafAlignment spa : curList)
		{
			String newKey = spa.contigName;
			IncludeContained.addInit(res, newKey, spa);
		}
	}
	return res;
}
static class NovelAdjacency implements Comparable<NovelAdjacency>
{
	String contig1, contig2;
	int pos1, pos2;
	int strand;
	String read;
	double weight;
	boolean prefix1, prefix2;
	int length1, length2;
	NovelAdjacency(String c1, String c2, int p1, int p2, int l1, int l2, String rr, double ww, boolean pr1, boolean pr2)
	{
		contig1 = c1;
		contig2 = c2;
		pos1 = p1;
		pos2 = p2;
		length1 = l1;
		length2 = l2;
		read = rr;
		weight = ww;
		prefix1 = pr1;
		prefix2 = pr2;
		if(contig1.compareTo(contig2) > 0 || (contig1.equals(contig2) && pos1 > pos2))
		{
			String tmp = contig1;
			contig1 = contig2;
			contig2 = tmp;
			int tmppos = pos1;
			pos1 = pos2;
			pos2 = tmppos;
			boolean tmppr = prefix1;
			prefix1 = prefix2;
			prefix2 = tmppr;
			int tmplen = length1;
			length1 = length2;
			length2 = tmplen;
		}
	}
	public String toString()
	{
		return "Novel adjacency: " + contig1 + " " + pos1 + " " + length1 + " " 
				+ contig2 + " " + pos2 + " " + length2 + " " + read + " "
				+ weight + " " + prefix1 + " " + prefix2;
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
			IncludeContained.addInit(res, readName, cur);
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
				IncludeContained.addInit(breakPositionMap, na.contig1, na.pos1);
				lengthMap.put(na.contig1, na.length1);
			}
			if(na.pos2 > buffer && na.pos2 + buffer < na.length2)
			{
				IncludeContained.addInit(breakPositionMap, na.contig2, na.pos2);
				lengthMap.put(na.contig2, na.length2);
			}
		}
		for(String contigName : breakPositionMap.keySet())
		{
			contigNames.add(contigName);
			ArrayList<Subcontig> scs = destroy(contigName, breakPositionMap.get(contigName));
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
				res.contigStart = Math.max(1, old.contigStart - sc.startPos);
				res.contigEnd = Math.min(sc.endPos, old.contigEnd) - sc.startPos;
				res.contigLength = sc.endPos - sc.startPos;
				return res;
			}
		}
		return null;
	}
	ArrayList<Subcontig> destroy(String contigName, ArrayList<Integer> breakPositions)
	{
		ArrayList<Subcontig> res = new ArrayList<>();
		Collections.sort(breakPositions);
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