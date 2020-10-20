import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Scanner;

public class TsvToAgp {
	static String metadataFn = "";
	static String piecesFn = "";
	static String ofn = "";
	static String readsFn = "";
	static String fastaOfn = "";
	static String scaffoldsFn = "";
	
	/*
	 * Parse command line arguments
	 */
	static void parseArgs(String[] args)
	{
		for(String arg : args)
		{
			int equalsIdx = arg.indexOf('=');
			if(equalsIdx == -1)
			{
				
			}
			else
			{
				String key = arg.substring(0, equalsIdx);
				String val = arg.substring(1 + equalsIdx);
				if(key.equalsIgnoreCase("metadata_file"))
				{
					metadataFn = val;
				}
				else if(key.equalsIgnoreCase("out_file"))
				{
					ofn = val;
				}
				else if(key.equalsIgnoreCase("contigs_file"))
				{
					piecesFn = val;
				}
				else if(key.equalsIgnoreCase("scaffolds_file"))
				{
					scaffoldsFn = val;
				}
				else if(key.equalsIgnoreCase("reads_file"))
				{
					readsFn = val;
				}
				else if(key.equalsIgnoreCase("fasta_out"))
				{
					fastaOfn = val;
				}
			}
		}
		
		if(metadataFn.length() == 0 || ofn.length() == 0 || piecesFn.length() == 0 || scaffoldsFn.length() == 0)
		{
			usage();
			System.exit(1);
		}
		
		if(fastaOfn.length() > 0 && readsFn.length() == 0)
		{
			usage();
			System.exit(1);
		}
	}
	
	/*
	 * Print a usage menu
	 */
	static void usage()
	{
		System.out.println();
		System.out.println("Usage: java -cp src TsvToAgp [args]");
		System.out.println("  Example: java -cp src TsvToAgp metadata_file=joins.tsv out_file=out.agp");
		System.out.println("    	contigs_file=contigs.fasta scaffolds_file=scaffolds.fasta");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  metadata_file    (String) - a table with joins from Grafter and their supporting reads");
		System.out.println("  contigs_file     (String) - the pre-Grafter contigs in FASTA format");
		System.out.println("  scaffolds_file   (String) - the post-Grafter scaffolds in FASTA format");
		System.out.println("  out_file         (String) - where to output the AGP file");
		System.out.println();
	}
	
	public static void main(String[] args) throws Exception
	{
		//metadataFn = "/home/mkirsche/physalis/readmeta_ragtag.tsv";
		//piecesFn = "/home/mkirsche/physalis/purged.fa";
		//ofn = "test_ragtag.agp";
		parseArgs(args);
		
		buildScaffolds();
	}
	
	/*
	 * Gets a map from piece name to its sequence length
	 */
	static HashMap<String, Integer> getContigLengths() throws Exception
	{
		HashMap<String, Integer> res = new HashMap<String, Integer>();
		Scanner input = new Scanner(new FileInputStream(new File(piecesFn)));
		String pieceName = "";
		StringBuilder sb = new StringBuilder("");
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.startsWith(">"))
			{
				if(pieceName.length() > 0)
				{
					res.put(pieceName, sb.length());
				}
				sb = new StringBuilder();
				pieceName = line.split(" ")[0].substring(1);
			}
			else
			{
				sb.append(line);
			}
		}
		if(pieceName.length() > 0)
		{
			res.put(pieceName, sb.length());
		}
		return res;
	}
	
	static HashMap<String, String> getContigToScaffoldMap() throws Exception
	{
		HashMap<String, String> res = new HashMap<String, String>();
		Scanner input = new Scanner(new FileInputStream(new File(scaffoldsFn)));
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.startsWith(">"))
			{
				String[] tokens = line.split(" ");
				String scaffoldName = tokens[0].substring(1);
				for(int i = 1; i<tokens.length; i++)
				{
					res.put(tokens[i], scaffoldName);
				}
			}
		}
		input.close();
		return res;
	}
	
	@SuppressWarnings("unchecked")
	static void buildScaffolds() throws Exception
	{
		PrintWriter out = new PrintWriter(new File(ofn));
		HashMap<String, Integer> contigLengths = getContigLengths();
		HashMap<String, String> contigToScaffold = getContigToScaffoldMap();
		Scanner input = new Scanner(new FileInputStream(new File(metadataFn)));
		
		ArrayList<Join> joins = new ArrayList<Join>();
		String header = input.nextLine();
		while(input.hasNext())
		{
			String line = input.nextLine();
			Join j = new Join(header, line);
			if(!j.sequenceUsed)
			{
				continue;
			}
			joins.add(j);
		}
		
		HashMap<String, Integer> contigToIdx = new HashMap<String, Integer>();
		for(Join j : joins)
		{
			if(!contigToIdx.containsKey(j.contigStart))
			{
				contigToIdx.put(j.contigStart, contigToIdx.size());
			}
			if(!contigToIdx.containsKey(j.contigEnd))
			{
				contigToIdx.put(j.contigEnd, contigToIdx.size());
			}
		}
		
		int n = contigToIdx.size();
		
		ArrayList<Join>[] graph = new ArrayList[n];
		for(int i = 0; i<n; i++)
		{
			graph[i] = new ArrayList<Join>();
		}
		
		for(Join j : joins)
		{
			int a = contigToIdx.get(j.contigStart);
			int b = contigToIdx.get(j.contigEnd);
			graph[a].add(j);
			graph[b].add(j);
		}
		
		boolean[] visited = new boolean[n];
		
		for(int i = 0; i<n; i++)
		{
			if(visited[i])
			{
				// This contig was already used in another scaffold.
				continue;
			}
			if(graph[i].size() == 0)
			{
				// Ignore contigs which are in no edges (shouldn't actually be possible)
				continue;
			}
			if(graph[i].size() == 2)
			{
				// If a contig has degree 2, it's in the middle of a scaffold and we want to start from an end.
				continue;
			}
						
			Queue<Integer> q = new LinkedList<Integer>();
			Queue<Join> joinq = new LinkedList<Join>();
			Queue<Integer> offsetq = new LinkedList<Integer>();
			q.add(i);
			offsetq.add(0);
			joinq.add(graph[i].get(0));
			visited[i] = true;
			int scaffoldLengthSoFar = 0;
			int partNum = 0;
			while(!q.isEmpty())
			{
				// The last contig
				int at = q.poll();
				int offset = offsetq.poll();
				
				// The last edge
				Join lastj = joinq.poll();
								
				// The name of the current contig
				String contigFrom = contigToIdx.get(lastj.contigStart).equals(at) ? lastj.contigStart : lastj.contigEnd;
				int fromLength = contigLengths.get(contigFrom); 
				String contigTo = contigToIdx.get(lastj.contigStart).equals(at) ? lastj.contigEnd : lastj.contigStart;
				
				// The index of the current contig
				int to = contigToIdx.get(contigTo);
				visited[to] = true;
				
				// Whether the last edge's previous contig used its prefix
				boolean lastPrefix = contigToIdx.get(lastj.contigStart).equals(at) ? lastj.startPrefix : lastj.endPrefix;
				
				String contigName = contigFrom;
				int contigStart = (lastPrefix ? offset : 0) + 1;
				int contigEnd = lastPrefix ? (fromLength - 1) : (fromLength - 1 - offset);
				int sequenceLength = contigEnd - contigStart + 1;
				String objectName = contigToScaffold.get(contigName);
				int objectStart = scaffoldLengthSoFar + 1;
				int objectEnd = scaffoldLengthSoFar + sequenceLength;
				scaffoldLengthSoFar += sequenceLength;
				char orientation = lastPrefix ? '-' : '+';
				partNum++;
				out.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
						objectName,
						objectStart,
						objectEnd,
						partNum,
						"W",
						contigName,
						contigStart,
						contigEnd,
						orientation
				);
				
				int lastReadStart = lastj.start;
				int lastReadEnd = lastj.end;
				if(lastReadStart < lastReadEnd)
				{
					sequenceLength = lastReadEnd - lastReadStart;
					objectStart = scaffoldLengthSoFar + 1;
					objectEnd = scaffoldLengthSoFar + sequenceLength;
					scaffoldLengthSoFar += sequenceLength;
					partNum++;
					out.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
							objectName,
							objectStart,
							objectEnd,
							partNum,
							"W",
							lastj.readName,
							lastReadStart + 1,
							lastReadEnd,
							lastj.strand == 0 ? "+" : "-"
					);
				}
				
				// Whether the last edge used the prefix of the current contig
				boolean curPrefix = contigToIdx.get(lastj.contigStart).equals(to) ? lastj.startPrefix : lastj.endPrefix;
				boolean found = false;
				for(Join j : graph[to])
				{
					boolean nextPrefix = contigToIdx.get(j.contigStart).equals(to) ? j.startPrefix : j.endPrefix;
					
					// Look for an edge whch uses the opposite side of the current contig from what the last edge used
					if(nextPrefix != curPrefix)
					{
						found = true;
						q.add(to);
						offsetq.add(Math.max(0, contigStart - contigEnd));
						joinq.add(j);
						break;
					}
				}
				
				if(!found)
				{
					// This means we're on the last contig, so add that to the scaffold
					boolean nextPrefix = contigToIdx.get(lastj.contigStart).equals(at) ? lastj.endPrefix : lastj.startPrefix;
					
					contigName = contigTo;
					contigStart = (nextPrefix ? offset : 0) + 1;
					int toLength = contigLengths.get(contigTo); 
					contigEnd = nextPrefix ? (toLength - 1) : (toLength - 1 - offset);
					sequenceLength = contigEnd - contigStart + 1;
					objectName = contigToScaffold.get(contigName);
					objectStart = scaffoldLengthSoFar + 1;
					objectEnd = scaffoldLengthSoFar + sequenceLength;
					scaffoldLengthSoFar += sequenceLength;
					orientation = nextPrefix ? '+' : '-';
					partNum++;
					out.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
							objectName,
							objectStart,
							objectEnd,
							partNum,
							"W",
							contigName,
							contigStart,
							contigEnd,
							orientation
					);
				}
			}
		}
		for(String s : contigLengths.keySet())
		{
			if(!contigToIdx.containsKey(s))
			{
				out.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
						s,
						1,
						contigLengths.get(s),
						1,
						"W",
						s,
						1,
						contigLengths.get(s),
						"+"
				);
			}
		}
		out.println();
		out.close();
	}
	
	static class Join
	{
		String readName;
		int start, end;
		String contigStart, contigEnd;
		boolean startPrefix, endPrefix;
		boolean sequenceUsed;
		int strand;
		Join(String header, String line)
		{
			String[] headerTokens = header.split("\t");
			String[] lineTokens = line.split("\t");
			for(int i = 0; i<headerTokens.length; i++)
			{
				if(headerTokens[i].equals("READNAME"))
				{
					readName = lineTokens[i].trim();
				}
				else if(headerTokens[i].equals("START"))
				{
					start = Integer.parseInt(lineTokens[i]);
				}
				else if(headerTokens[i].equals("CONTIG_START"))
				{
					contigStart = lineTokens[i];
				}
				else if(headerTokens[i].equals("START_PREFIX"))
				{
					startPrefix = lineTokens[i].equalsIgnoreCase("yes");
				}
				else if(headerTokens[i].equals("END"))
				{
					end = Integer.parseInt(lineTokens[i]);
				}
				else if(headerTokens[i].equals("CONTIG_END"))
				{
					contigEnd = lineTokens[i];
				}
				else if(headerTokens[i].equals("END_PREFIX"))
				{
					endPrefix = lineTokens[i].equalsIgnoreCase("yes");
				}
				else if(headerTokens[i].equals("SEQUENCE_USED"))
				{
					sequenceUsed = lineTokens[i].equalsIgnoreCase("yes");
				}
				else if(headerTokens[i].equals("STRAND"))
				{
					strand = Integer.parseInt(lineTokens[i]);
				}
			}
		}
	}
}
