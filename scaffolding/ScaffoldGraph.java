package scaffolding;

import java.util.*;

/*
 * Map contig name -> contig name -> list of alignments
 */
public class ScaffoldGraph {
	HashMap<String,ArrayList<Alignment>[]> adj;
	ScaffoldGraph()
	{
		adj =  new HashMap<>();
	}
	@SuppressWarnings("unchecked")
	void addEdge(String from, String to, String readName, int fromEnd, int toStart, int readLength, boolean fromPrefix, boolean toPrefix, double weight)
	{
		String[] contigs = new String[] {from, to};
		for(String key : contigs)
		{
			if(!adj.containsKey(key))
			{
				adj.put(key, new ArrayList[2]);
				for(int i = 0; i<adj.get(key).length; i++)
				{
					adj.get(key)[i] = new ArrayList<>();
				}
			}
		}
				
		// Add forward edge
		adj.get(from)[0].add(new Alignment(to, readName, fromEnd, toStart, fromPrefix, toPrefix, 0, weight));
		
		// Add reverse edge
		adj.get(to)[1].add(new Alignment(from, readName, readLength - toStart, readLength - fromEnd, toPrefix, fromPrefix, 1, weight));
	}

static class Alignment
{
	String to;
	String read;
	int myReadEnd;
	int theirReadStart;
	boolean myContigPrefix;
	boolean theirContigPrefix;
	int strand;
	double weight;
	Alignment(String tt, String rr, int mre, int trs, boolean mcp, boolean tcp, int ss, double ww)
	{
		to = tt;
		read = rr;
		myReadEnd = mre;
		theirReadStart = trs;
		myContigPrefix = mcp;
		theirContigPrefix = tcp;
		strand = ss;
		weight = ww;
	}
}
}
