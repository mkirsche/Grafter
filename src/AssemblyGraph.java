/*
 * Stores an assembly graph parsed from a GFA file
 */

import java.io.File;
import java.io.FileInputStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;

public class AssemblyGraph {
	HashMap<String, HashSet<String>> edges;
	HashMap<String, Integer> contigLengths;
	
	AssemblyGraph(String fn) throws Exception
	{
		edges = new HashMap<String, HashSet<String>>();
		contigLengths = new HashMap<String, Integer>();
		Scanner input = new Scanner(new FileInputStream(new File(fn)));
		while(input.hasNext())
		{
			String line = input.nextLine();
			String[] tokens = line.split("\t");
			if(tokens[0].equalsIgnoreCase("H"))
			{
				continue;
			}
			else if(tokens[0].equalsIgnoreCase("S"))
			{
				String contig = tokens[1];
				int length = Integer.parseInt(tokens[3].substring(3));
				contigLengths.put(contig, length);
			}
			else if(tokens[0].equalsIgnoreCase("L"))
			{
				String contig1 = tokens[1], contig2 = tokens[2];
				if(!edges.containsKey(contig1))
				{
					edges.put(contig1, new HashSet<String>());
				}
				if(!edges.containsKey(contig2))
				{
					edges.put(contig2, new HashSet<String>());
				}
				edges.get(contig1).add(contig2);
				edges.get(contig2).add(contig1);
			}
		}
	}
}
