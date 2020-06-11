import java.util.*;

import java.io.*;
public class ReadUtils {

/*
 * Add a key, value pair to a map from string to list, but initialize list if it's not already there
 */
static <T> void addToMap(HashMap<String, ArrayList<T>> map, String key, T val)
{
	if(!map.containsKey(key)) map.put(key, new ArrayList<T>());
	map.get(key).add(val);
}
/*
 * Get the reverse complement of a genomic sequence 
 */
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
/*
 * Reads a file consisting of pairs of string and produces
 * a map from the first string in each pair to the second
 */
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
	} catch(Exception e){
		return new HashMap<String, String>();
	}
}
/*
 * Write a map from string to string to a file as key-value pairs
 */
static void writeMap(String fn, HashMap<String, String> map) throws IOException
{
	PrintWriter out = new PrintWriter(new File(fn));
	for(String s : map.keySet())
	{
		out.println(s+" "+map.get(s));
	}
	out.close();
}
/*
 * Reads a file in FASTQ format and maps read names to their sequences
 */
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
/*
 * Reads a file in FASTA format and maps read names to their sequences
 */
@SuppressWarnings("resource")
static HashMap<String, String> getFastaMap(String fn, HashSet<String> names) throws IOException
{
	HashMap<String, String> res = new HashMap<String, String>();
	BufferedReader br = new BufferedReader(new FileReader(new File(fn)));
	String readName = br.readLine().split(" ")[0].substring(1);
	StringBuilder seq = new StringBuilder("");
	boolean useful = names.contains(readName);
	int curLength = 0;
	long totLength = 0;
	ArrayList<Integer> contigLengths = new ArrayList<Integer>();
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
				contigLengths.add(curLength);
				totLength += curLength;
				curLength = 0;
			}
			else
			{
				curLength += line.length();
				if(useful)
				{
					seq.append(line);
				}
			}
		} catch(Exception e) {
			break;
		}
	}
	
	contigLengths.add(curLength);
	totLength += curLength;
	
	if(useful)
	{
		res.put(readName, seq.toString());
	}
	
	assemblyStats(contigLengths, totLength);
	
	return res;
}

static void assemblyStats(ArrayList<Integer> contigLengths, long totLength)
{
	Collections.sort(contigLengths);
	System.err.println("Assembly Stats");
	System.err.println("Number of contigs: " + contigLengths.size());
	System.err.println("Total assembly length: " + totLength);
	
	int n50 = 0, nn50 = 0;
	
	long cumulativeLength = 0;
	for(int i = contigLengths.size() - 1; i >= 0; i--)
	{
		int curLength = contigLengths.get(i);
		cumulativeLength += curLength;
		n50 = curLength;
		nn50++;
		if(cumulativeLength * 2 >= totLength)
		{
			break;
		}
	}
	
	System.err.println("Contig N50: " + n50);
	System.err.println("Contig NN50: " + nn50);
	
	System.err.println("Longest contig: " +  contigLengths.get(contigLengths.size() - 1));
}
}
