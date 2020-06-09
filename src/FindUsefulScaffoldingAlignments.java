/*
Scaffolding assemblies with long-read alignments
 */
import java.util.*;
import java.io.*;
public class FindUsefulScaffoldingAlignments {
static String pafFn = "", outFn = "";
static int minQual = 0;
static String graphFn = "";

/*
 * Parse the command line arguments
 */
static void parseArgs(String[] args)
{
	for(String arg : args)
	{
		int equalsIdx = arg.indexOf("=");
		if(equalsIdx == -1)
		{

		}
		else
		{
			String field = arg.substring(0, equalsIdx);
			String val = arg.substring(1 + equalsIdx);
			if(field.equalsIgnoreCase("aln_fn"))
			{
				pafFn = val;
			}
			if(field.equalsIgnoreCase("out_file"))
			{
				outFn = val;
			}
			if(field.equalsIgnoreCase("graph_fn"))
			{
				graphFn = val;
			}
			if(field.equalsIgnoreCase("minq"))
			{
				minQual = Integer.parseInt(val);
			}
		}
	}
	if(pafFn.length() == 0)
	{
		usage();
		System.exit(1);
	}
	if(outFn.length() == 0)
	{
		usage();
		System.exit(1);
	}
}

/*
 * Print usage menu
 */
static void usage()
{
	System.out.println();
	System.out.println("Ultralong Scaffolding - module for finding useful alignments");
	System.out.println("Usage: java -cp src FindUsefulScaffoldingArguments [args]");
	System.out.println("  Example: java -cp src FindUsefulScaffoldingArguments aln_fn=aln.paf out_file=out.paf");
	System.out.println();
	System.out.println("Required args:");
	System.out.println("  aln_fn   (String) - a file containing the alignments of ultralong reads to contigs");
	System.out.println("  out_file (String) - the name of the file to output the useful alignments to");
	System.out.println();
	System.out.println("Optional args:");
	System.out.println("  minq     (int)    - the minimum quality of alignments needed to be kept");
	System.out.println("  graph_fn (String) - a GFA file containing an assembly graph, causing only alignments which are validated by the graph to be kept");
	System.out.println();
}

public static void  main(String[] args) throws Exception
{
	parseArgs(args);
	
	Scanner input = new Scanner(new FileInputStream(new File(pafFn)));
	PrintWriter out = new PrintWriter(new File(outFn));
	HashMap<String, ArrayList<PafAlignment>> alignmentsPerRead = new HashMap<>();
	
	// Read in the alignments and store them binned by read
	while(input.hasNext())
	{
		String line = input.nextLine();
		PafAlignment cur = new PafAlignment(line);
		if(cur.mapq < minQual)
		{
			continue;
		}
		String readName = cur.readName;
		if(!alignmentsPerRead.containsKey(readName)) alignmentsPerRead.put(readName, new ArrayList<PafAlignment>());
		alignmentsPerRead.get(readName).add(cur);
	}
	
	// If filtering based on an assembly graph, read in the GFA file
	AssemblyGraph graph = null;
	
	if(graphFn.length() > 0)
	{
		graph = new AssemblyGraph(graphFn);
	}
	
	// Now perform the actual filtering
	int usefulAlignmentCount = 0;
	for(String s : alignmentsPerRead.keySet()) 
	{
		// For now look at only reads which align to exactly two contigs - extend this later to catch larger paths
		if(alignmentsPerRead.get(s).size() == 2)
		{
			PafAlignment first = alignmentsPerRead.get(s).get(0);
			PafAlignment second = alignmentsPerRead.get(s).get(1);
			if(first.readStart > second.readStart)
			{
				PafAlignment tmp = first;
				first = second;
				second = tmp;
			}
			if(first.contigName.equals(second.contigName))
			{
				continue;
			}
			if(Math.abs(first.readEnd - first.readStart) < 10000 || Math.abs(second.readEnd - second.readStart) < 10000)
			{
				continue;
			}
			
			int maxHanging = 100;
			if(first.contigStart > maxHanging && first.contigEnd < first.contigLength - maxHanging)
			{
				continue;
			}
			if(second.contigStart > maxHanging && second.contigEnd < second.contigLength - maxHanging)
			{
				continue;
			}
			if(first.readEnd < second.readStart)
			{
				continue;
			}
			if(graphFn.length() > 0 && (!graph.edges.containsKey(first.contigName) || !graph.edges.get(first.contigName).contains(second.contigName)))
			{
				continue;
			}
			
			out.println(first.line + "\n" + second.line);
			usefulAlignmentCount++;
		}
	}
	out.close();
	System.err.println("Number of useful alignments: " + usefulAlignmentCount);
}
/*
 * Takes in paf line assuming minimap2 call was:
 *     minimap2 <contigs> <reads>
 */
static class PafAlignment
{
	String readName, contigName;
	int readLength, readStart, readEnd;
	int contigLength, contigStart, contigEnd;
	char strand;
	String line;
	int mapq;
	// Call with a second parameter to denote that read and contig were flipped when calling minimap2
	PafAlignment(String line, int backwards)
	{
		this.line = line;
		String[] ss = line.split("\t");
		contigName = ss[0];
		contigLength = Integer.parseInt(ss[1]);
		contigStart = Integer.parseInt(ss[2]);
		contigEnd = Integer.parseInt(ss[3]);
		strand = ss[4].charAt(0);
		readName = ss[5];
		readLength = Integer.parseInt(ss[6]);
		readStart = Integer.parseInt(ss[7]);
		readEnd = Integer.parseInt(ss[8]);
		mapq = Integer.parseInt(ss[11]);
	}
	PafAlignment(String line)
	{
		this.line = line;
		String[] ss = line.split("\t");
		readName = ss[0];
		readLength = Integer.parseInt(ss[1]);
		readStart = Integer.parseInt(ss[2]);
		readEnd = Integer.parseInt(ss[3]);
		strand = ss[4].charAt(0);
		contigName = ss[5];
		contigLength = Integer.parseInt(ss[6]);
		contigStart = Integer.parseInt(ss[7]);
		contigEnd = Integer.parseInt(ss[8]);
		mapq = Integer.parseInt(ss[11]);
	}
	PafAlignment()
	{
		
	}
}

}
