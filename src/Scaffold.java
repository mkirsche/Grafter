import java.util.*;

import java.io.*;
public class Scaffold {
	static int maxHanging = 100;
	
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
					alnFn = val;
				}
				if(field.equalsIgnoreCase("fasta_fn"))
				{
					fastaFn = val;
				}
				if(field.equalsIgnoreCase("read_fn"))
				{
					readFn = val;
				}
				if(field.equalsIgnoreCase("read_map_file"))
				{
					readMapFile = val;
				}
				if(field.equalsIgnoreCase("contig_map_file"))
				{
					contigMapFile = val;
				}
				if(field.equalsIgnoreCase("out_file"))
				{
					outFile = val;
				}
				if(field.equalsIgnoreCase("max_hanging"))
				{
					maxHanging = Integer.parseInt(val);
				}
			}
		}
		if(alnFn.length() == 0 || fastaFn.length() == 0 || readFn.length() == 0)
		{
			usage();
			System.exit(1);
		}
		if(readMapFile.length() == 0 || contigMapFile.length() == 0 || outFile.length() == 0)
		{
			usage();
			System.exit(1);
		}
	}
	
	static void usage()
	{
		System.out.println();
		System.out.println("Ultralong Scaffolding");
		System.out.println("Usage: java -cp src Scaffold [args]");
		System.out.println("  Example: java -cp src Scaffold aln_fn=aln.paf fasta_fn=contigs.fasta read_fn=reads.fastq");
		System.out.println("    read_map_file=useful_reads.paf contig_map_file=useful_contigs.paf out_file=out.fasta");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  aln_fn          (String) - a file containing the alignments of ultralong reads to contigs");
		System.out.println("  fasta_fn        (String) - the contigs in FASTA format");
		System.out.println("  read_fn         (String) - the ultralong reads in FASTQ format");
		System.out.println("  read_map_file   (String) - TODO");
		System.out.println("  contig_map_file (String) - TODO");
		System.out.println("  out_file        (String) - the name of the file to output the scaffolded contigs to");
		System.out.println();
		System.out.println("Optional args");
		System.out.println("  max_hanging     (int)    - the maximum amount by which the end of a contig can exceed the alignment and still be joined");
		System.out.println();
	}
	
static String alnFn = "", fastaFn = "", readFn = "", readMapFile = "", contigMapFile = "", outFile = "";
public static void main(String[] args) throws IOException
{
	parseArgs(args);
	
	ArrayList<FindUsefulScaffoldingAlignments.PafAlignment> als = new ArrayList<FindUsefulScaffoldingAlignments.PafAlignment>();
	HashSet<String> readNames = new HashSet<String>();
	HashSet<String> contigNames = new HashSet<String>();
	
	HashMap<String, String> readMap, contigMap;
	Scanner input = new Scanner(new FileInputStream(new File(alnFn)));
	
	/*
	 * Read paf alignments and store names of relevant reads/contigs
	 */
	while(input.hasNext())
	{
		FindUsefulScaffoldingAlignments.PafAlignment cur = new FindUsefulScaffoldingAlignments.PafAlignment(input.nextLine());
		als.add(cur);
		readNames.add(cur.readName);
		contigNames.add(cur.contigName);
	}
	
	System.err.println("Found " + readNames.size() + " relevant reads");
	System.err.println("Found " + contigNames.size() + " relevant contigs");
	/*
	 * Get sequences for relevant reads/sequences
	 */
	if((readMap = readMap(readMapFile)).size() == 0)
	{
		readMap = getFastqMap(readFn, readNames);
		writeMap(readMapFile, readMap);
	}
	if((contigMap = readMap(contigMapFile)).size() == 0)
	{
		contigMap = getFastaMap(fastaFn, contigNames);
		writeMap(contigMapFile, contigMap);
	}
	
	HashMap<String, HashSet<String>> countPrefix = new HashMap<String, HashSet<String>>();
	HashMap<String, HashSet<String>> countSuffix = new HashMap<String, HashSet<String>>();
	
	for(int i = 0; i<als.size(); i += 2)
	{
		FindUsefulScaffoldingAlignments.PafAlignment first = als.get(i), second = als.get(i+1);
		
		if(first.contigStart < maxHanging)
		{
			mapAdd(countPrefix, first.contigName, second.contigName);
		}
		if(first.contigEnd + maxHanging > first.contigLength)
		{
			mapAdd(countSuffix, first.contigName, second.contigName);
		}
		if(second.contigStart < maxHanging)
		{
			mapAdd(countPrefix, second.contigName, first.contigName);
		}
		if(second.contigEnd + maxHanging > second.contigLength)
		{
			mapAdd(countSuffix, second.contigName, first.contigName);
		}
	}
	
	PrintWriter out = new PrintWriter(new File(outFile));
	
	HashSet<String> joinedPref = new HashSet<String>();
	HashSet<String> joinedSuff = new HashSet<String>();
	
	HashMap<String, String> scaffoldStarters = new HashMap<String, String>();
	HashMap<String, String> scaffoldEnders = new HashMap<String, String>();
	
	HashSet<String> scaffoldNames = new HashSet<String>();
	
	HashMap<String, String> scaffoldID = new HashMap<String, String>();
	
	int joins = 0;
	
	for(int i = 0; i<als.size(); i += 2)
	{
		FindUsefulScaffoldingAlignments.PafAlignment first = als.get(i), second = als.get(i+1);
		boolean firstPref = first.contigStart < maxHanging;
		boolean firstSuff = first.contigEnd + maxHanging > first.contigLength;
		boolean secondPref = second.contigStart < maxHanging;
		boolean secondSuff = second.contigEnd + maxHanging > second.contigLength;
		
		boolean fp = firstPref;
		boolean sp = secondPref;
		
		if(firstPref && firstSuff) continue;
		if(secondPref && secondSuff) continue;
		if(firstPref && joinedPref.contains(first.contigName)) continue;
		if(firstSuff && joinedSuff.contains(first.contigName)) continue;
		if(secondPref && joinedPref.contains(second.contigName)) continue;
		if(secondSuff && joinedSuff.contains(second.contigName)) continue;
		
		String firstName = first.contigName;
		String secondName = second.contigName;
		
		if(scaffoldID.containsKey(firstName) && scaffoldID.containsKey(secondName) 
				&& scaffoldID.get(firstName).equals(scaffoldID.get(secondName)))
		{
			continue;
		}
		
		System.err.println("Candidate join between " + firstName + " and " + secondName + " " +firstPref+" "+secondPref);
		
		if(!contigMap.containsKey(first.contigName))
		{
			if(scaffoldStarters.containsKey(first.contigName))
			{
				firstPref = true;
				first.contigName = scaffoldStarters.get(first.contigName);
			}
			else if(scaffoldEnders.containsKey(first.contigName))
			{
				firstPref = false;
				first.contigName = scaffoldEnders.get(first.contigName);
			}
			else
			{
				continue;
			}
		}
		
		if(!contigMap.containsKey(second.contigName))
		{
			if(scaffoldStarters.containsKey(second.contigName))
			{
				secondPref = true;
				second.contigName = scaffoldStarters.get(second.contigName);
			}
			else if(scaffoldEnders.containsKey(second.contigName))
			{
				secondPref = false;
				second.contigName = scaffoldEnders.get(second.contigName);
			}
			else
			{
				continue;
			}
		}
		
		System.err.println("Joining");
		
		if(scaffoldStarters.containsKey(firstName))
		{
			scaffoldStarters.remove(firstName);
		}
		if(scaffoldStarters.containsKey(secondName))
		{
			scaffoldStarters.remove(secondName);
		}
		if(scaffoldEnders.containsKey(firstName))
		{
			scaffoldEnders.remove(firstName);
		}
		if(scaffoldEnders.containsKey(secondName))
		{
			scaffoldEnders.remove(secondName);
		}
		if(fp) joinedPref.add(firstName);
		else joinedSuff.add(firstName);
		
		if(sp) joinedPref.add(secondName);
		else joinedSuff.add(secondName);
		
		String firstSeq = contigMap.get(first.contigName);
		String secondSeq = contigMap.get(second.contigName);
		
		String seq = firstPref ? reverseComplement(firstSeq) : firstSeq;
		String seq2 = secondPref ? secondSeq : reverseComplement(secondSeq);
		
		int overlap1 = first.contigEnd - first.contigStart;
		int overlap2 = second.contigEnd - second.contigStart;
		int extraAligned = first.readEnd - second.readStart;
		
		String nname = first.contigName + "&" + second.contigName;
		String[] subcontigs = nname.split("&");
		for(String s : subcontigs) scaffoldID.put(s, nname);
		scaffoldStarters.put(subcontigs[0], nname);
		scaffoldEnders.put(subcontigs[subcontigs.length-1], nname);
		if(scaffoldNames.contains(first.contigName))
		{
			scaffoldNames.remove(first.contigName);
		}
		if(scaffoldNames.contains(second.contigName))
		{
			scaffoldNames.remove(second.contigName);
		}
		scaffoldNames.add(nname);
		String stitched = stitch(seq, seq2, overlap1, overlap2, extraAligned);
		joins++;
		
		contigMap.remove(first.contigName);
		contigMap.remove(second.contigName);
		contigMap.put(nname, stitched);
	}
	
	for(String s : scaffoldNames)
	{
		out.println(">" + s);
		out.println(contigMap.get(s));
	}
	
	
	System.err.println("Joins performed: " + joins);
	out.close();
}
static String stitch(String seq1, String seq2, int overlap1, int overlap2, int extra)
{
	return seq1 + seq2.substring(extra);
}
static String reverseComplement(String s)
{
	int n = s.length();
	char[] res = new char[n];
	for(int i = 0; i<n; i++)
	{
		char c = s.charAt(n-1-i);
		if(c == 'a' || c == 'A') res[i] = (char)(c + 'T' - 'A');
		if(c == 'c' || c == 'C') res[i] = (char)(c + 'G' - 'C');
		if(c == 'g' || c == 'G') res[i] = (char)(c + 'C' - 'G');
		if(c == 't' || c == 'T') res[i] = (char)(c + 'A' - 'T');
	}
	return new String(res);
}
static void mapAdd(HashMap<String, HashSet<String>> map, String s, String t)
{
	if(!map.containsKey(s)) map.put(s, new HashSet<String>());
	map.get(s).add(t);
}
static HashMap<String, String> readMap(String fn) throws IOException
{
	System.err.println("Reading map from " + fn);
	try {
		Scanner input = new Scanner(new FileInputStream(new File(fn)));
		HashMap<String, String> res = new HashMap<String, String>();
		while(input.hasNext())
		{
			String line = input.nextLine();
			int idx = line.indexOf(' ');
			String key = line.substring(0, idx);
			String val = line.substring(idx+1);
			res.put(key, val);
		}
		return res;
	} catch(Exception e) {
		return new HashMap<String, String>();
	}
}
static void writeMap(String fn, HashMap<String, String> map) throws IOException
{
	PrintWriter out = new PrintWriter(new File(fn));
	for(String s : map.keySet())
	{
		out.println(s+" "+map.get(s));
	}
	out.close();
}
@SuppressWarnings("resource")
static HashMap<String, String> getFastqMap(String fn, HashSet<String> names)  throws IOException
{
	HashMap<String, String> res = new HashMap<String, String>();
	BufferedReader br = new BufferedReader(new FileReader(new File(fn)));
	while(true)
	{
		try {
			String name = br.readLine().substring(1).split(" ")[0];
			if(names.contains(name)) res.put(name, br.readLine());
			else br.readLine();
			for(int i = 0; i<2; i++) br.readLine();
		} catch(Exception e) {
			break;
		}
	}
	return res;
}
@SuppressWarnings("resource")
static HashMap<String, String> getFastaMap(String fn, HashSet<String> names) throws IOException
{
	HashMap<String, String> res = new HashMap<String, String>();
	BufferedReader br = new BufferedReader(new FileReader(new File(fn)));
	String readName = br.readLine().split(" ")[0].substring(1);
	StringBuilder seq = new StringBuilder("");
	boolean useful = names.contains(readName);
	while(true)
	{
		try {
			String line = br.readLine();
			if(line.charAt(0) == '>')
			{
				// process last read
				if(useful)
				{
					res.put(readName, seq.toString());
					seq = new StringBuilder("");
				}
				// new read name
				readName = line.split(" ")[0].substring(1);
				useful = names.contains(readName);
			}
			else if(useful)
			{
				seq.append(line);
			}
		} catch(Exception e) {
			break;
		}
	}
	if(useful)
	{
		res.put(readName, seq.toString());
	}
	return res;
}
}
