package scaffolding;

import java.util.*;

/*
 * Map contig name -> contig name -> list of alignments
 */
public class ScaffoldGraph {
	HashMap<String, HashMap<String, ArrayList<Alignment>>> adj;
	ScaffoldGraph()
	{
		adj =  new HashMap<>();
	}
	void addEdge(String from, String to, String readName, int fromEnd, int toStart, boolean fromPrefix, boolean toPrefix)
	{
		if(!adj.containsKey(from))
		{
			adj.put(from, new HashMap<String, ArrayList<Alignment>>());
		}
		HashMap<String, ArrayList<Alignment>> fromEdgeList = adj.get(from);
		if(!fromEdgeList.containsKey(to))
		{
			fromEdgeList.put(to, new ArrayList<Alignment>());
		}
		fromEdgeList.get(to).add(new Alignment(readName, fromEnd, toStart, fromPrefix, toPrefix));
	}

static class Alignment
{
	String read;
	int myReadEnd;
	int theirReadStart;
	boolean myContigPrefix;
	boolean theirContigPrefix;
	Alignment(String rr, int mre, int trs, boolean mcp, boolean tcp)
	{
		read = rr;
		myReadEnd = mre;
		theirReadStart = trs;
		myContigPrefix = mcp;
		theirContigPrefix = tcp;
	}
}
}
