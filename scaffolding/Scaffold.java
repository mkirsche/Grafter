package scaffolding;

import java.util.*;

import scaffolding.FindUsefulScaffoldingAlignments.PafAlignment;

import java.io.*;
public class Scaffold {
public static void main(String[] args) throws IOException
{
	String fn = "/scratch/groups/mschatz1/mkirsche/ultralong/ccs/rel2_200kplus_ccs_useful.paf";
	String fastaFn = "/scratch/groups/mschatz1/mkirsche/ultralong/ccs/paternal_and_unknown.contigs.mmpoa.fa";
	String readFn = "/scratch/groups/mschatz1/mkirsche/ultralong/rel2_200kplus.fastq";
	String readMapFile = "/scratch/groups/mschatz1/mkirsche/ultralong/ccs/readmap_paternal.txt";
	String contigMapFile = "/scratch/groups/mschatz1/mkirsche/ultralong/ccs/contigmap_paternal.txt";
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	ArrayList<PafAlignment> als = new ArrayList<PafAlignment>();
	HashSet<String> readNames = new HashSet<String>();
	HashSet<String> contigNames = new HashSet<String>();
	HashMap<String, String> readMap, contigMap;
	
	/*
	 * Read paf alignments and store names of relevant reads/contigs
	 */
	while(input.hasNext())
	{
		PafAlignment cur = new PafAlignment(input.nextLine());
		als.add(cur);
		readNames.add(cur.readName);
		contigNames.add(cur.contigName);
	}
	
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
	
	for(int i = 0; i<als.size(); i += 2)
	{
		PafAlignment first = als.get(i), second = als.get(i+1);
		System.out.println(first+" "+second);
	}
}
static HashMap<String, String> readMap(String fn) throws IOException
{
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
static HashMap<String, String> getFastqMap(String fn, HashSet<String> names)  throws IOException
{
	HashMap<String, String> res = new HashMap<String, String>();
	BufferedReader br = new BufferedReader(new FileReader(new File(fn)));
	while(true)
	{
		try {
			String name = br.readLine().substring(1);
			if(names.contains(name)) res.put(name, br.readLine());
			else br.readLine();
			for(int i = 0; i<2; i++) br.readLine();
		} catch(Exception e) {
			break;
		}
	}
	return res;
}
static HashMap<String, String> getFastaMap(String fn, HashSet<String> names) throws IOException
{
	HashMap<String, String> res = new HashMap<String, String>();
	BufferedReader br = new BufferedReader(new FileReader(new File(fn)));
	boolean keep = false;
	String readName = br.readLine().substring(1);
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
				readName = line.substring(1);
				useful = names.contains(readName);
			}
			else
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
