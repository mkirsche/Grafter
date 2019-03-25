package scaffolding;

import java.util.*;

public class ScaffoldGraphBuilder {
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
				ReadUtils.addInit(prefEdges, a.to, a);
			}
			else
			{
				ReadUtils.addInit(suffEdges, a.to, a);
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
		
		if(best.size() < IncludeContained.minReadSupport || bestTotalWeight < IncludeContained.minWeightSupport) return null;
		
		
		ScaffoldGraph.Alignment res = best.get(0);
		
		return res;
	}
}
