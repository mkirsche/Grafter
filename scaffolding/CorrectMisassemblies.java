/*
 * Code for correcting misassemblies - not yet working
 * TODO Make stricter threshold for calling a misassembly
 */

package scaffolding;

import java.util.*;

import scaffolding.IncludeContained.SortablePafAlignment;

public class CorrectMisassemblies {
static class BrokenContig
{
	String oldName;
	ArrayList<Integer> splits;
	BrokenContig(String oldName, ArrayList<Integer> splits)
	{
		this.oldName = oldName;
		this.splits = splits;
	}
}
static ArrayList<BrokenContig> correctAllContigs(HashMap<String, ArrayList<IncludeContained.SortablePafAlignment>> byRead)
{
	HashMap<String, ArrayList<IncludeContained.SortablePafAlignment>> byContig = reindex(byRead);
	
	ArrayList<BrokenContig> res = new ArrayList<BrokenContig>();
	
	for(String s : byContig.keySet())
	{
		BrokenContig cur = correctContig(s, byContig.get(s));
		if(cur != null)
		{
			res.add(cur);
		}
	}
	
	return res;
}
static BrokenContig correctContig(String name, ArrayList<IncludeContained.SortablePafAlignment> alignments)
{
	BrokenContig res = null;
	if(alignments.size() < 25)
	{
		return res;
	}
	
	// Sort by read start position
	Collections.sort(alignments, new Comparator<SortablePafAlignment>() {

		@Override
		public int compare(SortablePafAlignment a, SortablePafAlignment b) {
			if(a.contigStart == b.contigStart)
			{
				return a.compareTo(b);
			}
			return Integer.compare(a.contigStart, b.contigStart);
		}
	});
	
	ArrayList<Integer> splits = new ArrayList<Integer>();
	int lastEnd = -1;
	HashMap<String, Integer> firstSeenEachRead = new HashMap<String, Integer>();
	for(SortablePafAlignment spa : alignments)
	{
		if(lastEnd == -1 || spa.contigStart < lastEnd + 100000)
		{
			lastEnd = Math.max(lastEnd, spa.contigEnd);
		}
		else
		{
			if(firstSeenEachRead.containsKey(spa.contigName) && firstSeenEachRead.get(spa.contigName) < splits.size())
			{
				while(splits.size() > firstSeenEachRead.get(spa.contigName))
				{
					splits.remove(splits.size()-1);
				}
				for(String s : firstSeenEachRead.keySet())
				{
					if(firstSeenEachRead.get(s) > splits.size())
					{
						firstSeenEachRead.put(s, splits.size());
					}
				}
			}
			else
			{
				splits.add(lastEnd);
				splits.add(spa.contigStart);
			}
			lastEnd = Math.max(lastEnd, spa.contigEnd);
		}
		if(!firstSeenEachRead.containsKey(spa.readName))
		{
			firstSeenEachRead.put(spa.readName, splits.size());
		}
	}
	if(splits.size() > 0)
	{
		System.out.println(name+" "+splits);
	}
	
	if(splits.size() > 0)
	{
		res = new BrokenContig(name, splits);
	}
	
	return res;
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
}