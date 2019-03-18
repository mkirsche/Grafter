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
		Alignment al = new Alignment(to, readName, fromEnd, toStart, fromPrefix, toPrefix, 0, weight, readLength);
		adj.get(from)[0].add(al);
		
		// Add reverse edge
		adj.get(to)[1].add(al.reverse(from));
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
	int readLength;
	Alignment(String tt, String rr, int mre, int trs, boolean mcp, boolean tcp, int ss, double ww, int rl)
	{
		to = tt;
		read = rr;
		myReadEnd = mre;
		theirReadStart = trs;
		myContigPrefix = mcp;
		theirContigPrefix = tcp;
		strand = ss;
		weight = ww;
		this.readLength = rl;
	}
	
	Alignment reverse(String from)
	{
		return new Alignment(from, read, readLength - theirReadStart, readLength - myReadEnd, 
				theirContigPrefix, myContigPrefix, 1 - strand, weight, readLength);
	}
}
}
