import java.util.*;

public class ScaffoldGraphBuilder {
	
	static boolean stillValid(ScaffoldGraph.Alignment sga, ScaffoldGraph.Scaffolding curScaffolds, HashMap<String, String> lastToFirst)
	{
		HashSet<String> usedContigs = curScaffolds.usedContigs;
		HashMap<String, ArrayDeque<ScaffoldGraph.Alignment>> scaffoldEdges = curScaffolds.scaffoldEdges;
		String to = sga.to;
		String from = sga.from;
		
		// Make sure from is end of its scaffold
		if(usedContigs.contains(from) && !lastToFirst.containsKey(from))
		{
			return false;
		}
		
		// Make sure the destination is on one of the ends of its scaffold
		if(usedContigs.contains(to) && !scaffoldEdges.containsKey(to) && !lastToFirst.containsKey(to))
		{
			return false;
		}
					
		// Make sure that the destination isn't the beginning of the scaffold the edge is coming from
		if(scaffoldEdges.containsKey(to) && scaffoldEdges.get(to).peekLast().to.equals(from))
		{
			return false;
		}
					
		// Make sure that the destination isn't the end of the scaffold the edge is coming from
		if(lastToFirst.containsKey(to) && lastToFirst.get(to).equals(from))
		{
			return false;
		}
		
		// Make sure this edge doesn't use the same side of the destination as an existing edge
		if(scaffoldEdges.containsKey(to) && scaffoldEdges.get(to).peekFirst().myContigPrefix == sga.theirContigPrefix)
		{
			return false;
		}
		if(lastToFirst.containsKey(to) && scaffoldEdges.get(lastToFirst.get(to)).peekLast().theirContigPrefix == sga.theirContigPrefix)
		{
			return false;
		}
		return true;
	}
	
	static PriorityQueue<ScaffoldGraph.Alignment> getAllSortedEdges(ScaffoldGraph sg)
	{
		PriorityQueue<ScaffoldGraph.Alignment> res = new PriorityQueue<>();
		for(String s : sg.adj.keySet())
		{
			if(Settings.VERBOSE)
			{
				System.err.println("Searching for edges from " + s);
			}
			for(int strand = 0; strand < 2; strand++)
			{
				ArrayList<ScaffoldGraph.Alignment> als = sg.adj.get(s)[strand];
				HashMap<String, ArrayList<ScaffoldGraph.Alignment>> prefEdges = new HashMap<>();
				HashMap<String, ArrayList<ScaffoldGraph.Alignment>> suffEdges = new HashMap<>();
				for(ScaffoldGraph.Alignment a : als)
				{
					if(a.myContigPrefix)
					{
						ReadUtils.addToMap(prefEdges, a.to, a);
					}
					else
					{
						ReadUtils.addToMap(suffEdges, a.to, a);
					}
				}
				
				for(String to : prefEdges.keySet())
				{
					ArrayList<ScaffoldGraph.Alignment> al = prefEdges.get(to);
					ArrayList<ScaffoldGraph.Alignment> valid = new ArrayList<>();
					double totalWeight = 0;
					ArrayList<ScaffoldGraph.ReadInterval> intervals = new ArrayList<ScaffoldGraph.ReadInterval>();
					for(ScaffoldGraph.Alignment a : al)
					{
						totalWeight += a.weight;
						valid.add(a);
						a.from = s;
						ScaffoldGraph.ReadInterval ri = new ScaffoldGraph.ReadInterval(a);
						if(Settings.VERBOSE)
						{
							System.err.println("Adding read interval (from prefix): " + ri.readName
									+ " " + ri.from + " " + ri.to + " " + ri.start + " " + ri.end);
						}
						intervals.add(ri);
					}
					
					if(valid.size() > 0)
					{
						ScaffoldGraph.Alignment toAdd = valid.get(0);
						toAdd.from = s;
						toAdd.weight = totalWeight;
						toAdd.allReads = intervals;
						if(Settings.VERBOSE)
						{
							System.err.println("Adding consensus edge to graph: ");
							System.err.println(" from=" + toAdd.from + ", myContigPrefix=" + toAdd.myContigPrefix +
									", to=" + toAdd.to + ", theirContigPrefix=" + toAdd.theirContigPrefix);
						}
						res.add(toAdd);
					}
					
				}
				for(String to : suffEdges.keySet())
				{
					ArrayList<ScaffoldGraph.Alignment> valid = new ArrayList<>();
					ArrayList<ScaffoldGraph.Alignment> al = suffEdges.get(to);
					double totalWeight = 0;
					ArrayList<ScaffoldGraph.ReadInterval> intervals = new ArrayList<ScaffoldGraph.ReadInterval>();
					for(ScaffoldGraph.Alignment a : al)
					{
						totalWeight += a.weight;
						valid.add(a);
						a.from = s;
						ScaffoldGraph.ReadInterval ri = new ScaffoldGraph.ReadInterval(a);
						if(Settings.VERBOSE)
						{
							System.err.println("Adding read interval (from suffix): " + ri.readName
									+ " " + ri.from + " " + ri.to + " " + ri.start + " " + ri.end);
						}
						intervals.add(ri);
					}
					if(valid.size() > 0)
					{
						ScaffoldGraph.Alignment toAdd = valid.get(0);
						toAdd.from = s;
						toAdd.weight = totalWeight;
						toAdd.allReads = intervals;
						if(Settings.VERBOSE)
						{
							System.err.println("Adding consensus edge to graph: ");
							System.err.println(" from=" + toAdd.from + ", myContigPrefix=" + toAdd.myContigPrefix +
									", to=" + toAdd.to + ", theirContigPrefix=" + toAdd.theirContigPrefix);
						}
						res.add(toAdd);
					}
				}
			}
		}
		return res;
	}
	/*
	 * Gets the best alignment of another contig to follow a given contig
	 */
	/*
	static ScaffoldGraph.Alignment consensus(String from, ArrayList<ScaffoldGraph.Alignment> als, HashSet<String> usedContigs, HashMap<String, ArrayDeque<ScaffoldGraph.Alignment>> scaffoldEdges, HashMap<String, String> lastToFirst)
	{
		HashMap<String, ArrayList<ScaffoldGraph.Alignment>> prefEdges = new HashMap<>();
		HashMap<String, ArrayList<ScaffoldGraph.Alignment>> suffEdges = new HashMap<>();
		for(ScaffoldGraph.Alignment a : als)
		{
			if(a.myContigPrefix)
			{
				ReadUtils.addToMap(prefEdges, a.to, a);
			}
			else
			{
				ReadUtils.addToMap(suffEdges, a.to, a);
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
		
		if(best.size() < Settings.MIN_READ_SUPPORT || bestTotalWeight < Settings.MIN_WEIGHT_SUPPORT) return null;
		
		
		ScaffoldGraph.Alignment res = best.get(0);
		
		return res;
	}
	*/
}
