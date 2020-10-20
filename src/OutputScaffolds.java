import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.HashMap;

public class OutputScaffolds {
	
	/*
	 * Output the overlap of contigs 
	 */
	static void outputGfa(String ofn, ScaffoldGraph graph, HashMap<String, String> contigSequences) throws Exception
	{
		PrintWriter out = new PrintWriter(new File(ofn));
		out.println("H\t1.0");
		for(String s : contigSequences.keySet())
		{
			out.println("S\t" + s + "\t*\tLN:" + contigSequences.get(s).length());
		}
		for(String from: graph.adj.keySet())
		{
			ArrayList<ScaffoldGraph.Alignment> alnList = graph.adj.get(from)[0];
			for(ScaffoldGraph.Alignment aln : alnList)
			{
				String to = aln.to;
				char fromStrand = aln.myContigPrefix ? '-' : '+';
				char toStrand = aln.theirContigPrefix ? '+' : '-';
				out.printf("%s\t%s\t%s\t%s\t%s\t%s\n", "L", from, fromStrand, to, toStrand, "*");
			}
		}
		out.close();
	}
	
	/*
	 * Prints the orientations of scaffold edges to a txt file
	 */
	static void printOrientations(HashMap<String, ArrayDeque<ScaffoldGraph.Alignment>> als, CorrectMisassemblies.ContigBreaker splitter) throws IOException
	{
		PrintWriter out = new PrintWriter(new File("orientations.txt"));
		for(String s : als.keySet())
		{
			ArrayDeque<ScaffoldGraph.Alignment> cur = als.get(s);
			String contigName = cur.peekFirst().from;
			String oldName = splitter.sourceMap.containsKey(contigName) ? splitter.sourceMap.get(contigName) : contigName;
			out.print(oldName + " " + (cur.peekFirst().myContigPrefix ? '-' : '+'));
			for(ScaffoldGraph.Alignment spa : cur)
			{
				contigName = spa.to;
				oldName = splitter.sourceMap.containsKey(contigName) ? splitter.sourceMap.get(contigName) : contigName;
				out.print(" " + oldName + " " + (cur.peekFirst().theirContigPrefix ? '+' : '-'));
			}
			out.println();
		}
		out.close();
	}

	/*
	 * Create a Fasta header line for a scaffold based on the names of contigs which make it up
	 * format is >scaffold(index)_grafter contig1 contig2 contig3 ...
	 */
	static String createHeaderLine(int index, ArrayDeque<String> contigs, CorrectMisassemblies.ContigBreaker splitter)
	{
		StringBuilder res = new StringBuilder("");
		res.append(">scaffold" + index + "_grafter");
		
		for(String s : contigs)
		{
			if(splitter.sourceMap.containsKey(s))
			{
				res.append(" " + splitter.sourceMap.get(s));
			}
			else
			{
				res.append(" " + s);
			}
		}
		return res.toString();
	}
}
