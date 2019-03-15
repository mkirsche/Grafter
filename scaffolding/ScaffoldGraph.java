package scaffolding;

import java.util.*;

/*
 * Map contig name -> contig name -> list of alignments
 */
public class ScaffoldGraph {
	HashMap<String,ArrayList<Alignment>> adj;
	ScaffoldGraph()
	{
		adj =  new HashMap<>();
	}
	void addEdge(String from, String to, String readName, int fromEnd, int toStart, boolean fromPrefix, boolean toPrefix)
	{
		if(!adj.containsKey(from))
		{
			adj.put(from, new ArrayList<Alignment>());
		}
		adj.get(from).add(new Alignment(to, readName, fromEnd, toStart, fromPrefix, toPrefix));
	}

static class Alignment
{
	String to;
	String read;
	int myReadEnd;
	int theirReadStart;
	boolean myContigPrefix;
	boolean theirContigPrefix;
	Alignment(String tt, String rr, int mre, int trs, boolean mcp, boolean tcp)
	{
		to = tt;
		read = rr;
		myReadEnd = mre;
		theirReadStart = trs;
		myContigPrefix = mcp;
		theirContigPrefix = tcp;
	}
}
}
