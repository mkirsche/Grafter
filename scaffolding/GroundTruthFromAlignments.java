package scaffolding;

import java.util.*;
import java.io.*;

public class GroundTruthFromAlignments {
@SuppressWarnings("resource")
public static void main(String[] args) throws IOException
{
	String fn = "maternaltoref.paf";
	
	if(args.length > 0)
	{
		fn = args[0];
	}
	
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	
	HashMap<String, Chromosome> chrMap = new HashMap<>();
	
	HashMap<String, Mapping> contigNameToMapping = new HashMap<>();
	
	HashMap<String, Integer> contigToLength = new HashMap<>();
	
	while(input.hasNext())
	{
		String[] line = input.nextLine().split("\t");
		String contigName = line[0];
		boolean rev = line[4].equals("-");
		String chrName = line[5];
		
		if(chrName.contains("v"))
		{
			continue;
		}
		
		int chrStart = Integer.parseInt(line[7]);
		int chrEnd = Integer.parseInt(line[8]);
		
		contigToLength.put(contigName, Integer.parseInt(line[1]));
		
		int contigLength = Integer.parseInt(line[3]) - Integer.parseInt(line[2]);
		
		if(contigLength < 5000)
		{
			continue;
		}
		
		boolean primary = false;
		for(int i = 12; i<line.length; i++)
		{
			if(line[i].equals("tp:A:P"))
			{
				primary = true;
			}
		}
		
		if(!primary)
		{
			continue;
		}
		
		if(!chrMap.containsKey(chrName))
		{
			chrMap.put(chrName, new Chromosome(chrName));
		}
		
		//System.out.println(contigLength); 
		
		int oldLength = contigNameToMapping.containsKey(contigName) ? contigNameToMapping.get(contigName).contigLength : 0;
		
		if(contigLength > oldLength)
		{
			contigNameToMapping.put(contigName, new Mapping(contigName, chrName, chrStart, chrEnd, contigLength, rev));
		}
	}
	
	for(String contigName : contigNameToMapping.keySet())
	{
		Mapping cur = contigNameToMapping.get(contigName);
		chrMap.get(cur.chrName).mappings.add(cur);
	}
	
	PrintWriter out = new PrintWriter(System.out);
	
	int joins = 0;
	
	HashMap<String, String> contigToChr = new HashMap<>();
	HashMap<String, Integer> contigToPos = new HashMap<>();
	
	for(String chrName : chrMap.keySet())
	{
		Mapping last = null;
		int index = 0;
		for(Mapping m : chrMap.get(chrName).mappings)
		{
			contigToChr.put(m.contigName, chrName);
			contigToPos.put(m.contigName, index);
			if(last != null)
			{
				// Process adjacency between m and last
				//out.println(last.contigName+" "+m.contigName);
				joins++;
			}
			last = m;
			index++;
		}
	}
	
	System.out.println("Reference-based joins: " + joins);
	
	String newContigsFn = "new_contigs.fa";//"/home/mkirsche/eclipse-workspace/Ultralong/new_contigs.fa";
	
	input = new Scanner(new FileInputStream(new File(newContigsFn)));
	
	int chimeras = 0;
	
	HashSet<Misjoin> misjoins = new HashSet<>();
	
	while(input.hasNext())
	{
		String[] line = input.nextLine().split(" ");
		
		//Ignore sequence lines
		if(line.length == 0 || line[0].charAt(0) != '>') continue;
		
		//System.out.println(Arrays.toString(line));
		
		for(int i = 1; i<line.length-1; i++)
		{
			String chrA = contigToChr.get(line[i]);
			String chrB = contigToChr.get(line[i+1]);
			
			if(chrA == null || chrB == null)
			{
				continue;
			}
			
			if(!chrA.equals(chrB))
			{
				misjoins.add(new Misjoin(line[i], line[i+1]));
				chimeras++;
				System.out.println("Chimeric join: " + line[i] + " (length " + contigToLength.get(line[i]) + ") on " + chrA
						+" to "+line[i+1]+ " (length " + contigToLength.get(line[i+1]) + ") on "+chrB);
			}
		}
	}
	
	String pafFn = "alignments.paf";
	
	System.out.println("Chimera joins: " + chimeras);
	
	System.out.println("Building alignment index");
	
	AlignmentIndex ai = new AlignmentIndex(pafFn);
	
	for(Misjoin mj : misjoins)
	{
		ai.debugPrint(mj.c1, mj.c2);
	}
	
	out.close();
}
static class Chromosome
{
	String name;
	TreeSet<Mapping> mappings;
	Chromosome(String nn)
	{
		name = nn;
		mappings = new TreeSet<>();
	}
}
/*
 * Represents a misjoin between contigs c1 and c2
 */
static class Misjoin implements Comparable<Misjoin>
{
	String c1, c2;
	Misjoin(String cc1, String cc2)
	{
		c1 = cc1;
		c2 = cc2;
		if(c1.compareTo(c2) > 0)
		{
			String tmp = c1;
			c1 = c2;
			c2 = tmp;
		}
	}
	@Override
	public int compareTo(Misjoin o) {
		if(!c1.equals(o.c1)) return c1.compareTo(o.c1);
		return c2.compareTo(o.c2);
	}
}
/*
 * Represents a PAF file in a way which makes different types of queries easier
 */
static class AlignmentIndex
{
	HashMap<String, HashSet<Integer>> readToAlignments;
	HashMap<String, HashSet<Integer>> contigToAlignments;
	ArrayList<String> reads;
	ArrayList<String> contigs;
	ArrayList<String> lines;
	int n;
	AlignmentIndex(String fn) throws IOException
	{
		readToAlignments = new HashMap<>();
		contigToAlignments = new HashMap<>();
		reads =  new ArrayList<>();
		contigs  = new ArrayList<>();
		lines = new ArrayList<>();
		Scanner input = new Scanner(new FileInputStream(new File(fn)));
		while(input.hasNext())
		{
			addAlignment(input.nextLine());
		}
	}
	void addAlignment(String line)
	{
		if(line.length() == 0)
		{
			return;
		}
		IncludeContained.SortablePafAlignment spa = new IncludeContained.SortablePafAlignment(line);
		addInit(contigToAlignments, spa.contigName, n);
		addInit(readToAlignments, spa.readName, n);
		reads.add(spa.readName);
		contigs.add(spa.contigName);
		lines.add(line);
		n++;
	}
	static <T> void addInit(HashMap<String, HashSet<T>> map, String key, T val)
	{
		if(!map.containsKey(key)) map.put(key, new HashSet<T>());
		map.get(key).add(val);
	}
	static <T> HashSet<T> intersect(HashSet<T> a, HashSet<T> b)
	{
		HashSet<T> res = new HashSet<T>();
		for(T x : a)
		{
			if(b.contains(x))
			{
				res.add(x);
			}
		}
		return res;
	}
	HashSet<String> alnsToReads(HashSet<Integer> alns)
	{
		HashSet<String> res = new HashSet<>();
		for(int x : alns)
		{
			res.add(reads.get(x));
		}
		return res;
	}
	void printAlignmentsFromRead(String readName)
	{
		System.out.println("Alignments for " + readName + ":");
		HashSet<Integer> alignments = readToAlignments.get(readName);
		ArrayList<Integer> toPrint = new ArrayList<>();
		toPrint.addAll(alignments);
		Collections.sort(toPrint, compareAlignmentsByField(3, true, false));
		for(int x : toPrint)
		{
			System.out.println(lines.get(x));
		}
	}
	Comparator<Integer> compareAlignmentsByField(int i, boolean intCompare, boolean reverse)
	{
		return new Comparator<Integer>() {

			@Override
			public int compare(Integer a, Integer b) {
				String[] aSplit = lines.get(a).split("\t");
				String[] bSplit = lines.get(b).split("\t");
				
				if(intCompare)
				{
					int res = Integer.parseInt(aSplit[i-1]) - Integer.parseInt(bSplit[i-1]);
					return reverse ? -res : res;
				}
				else
				{
					int res = aSplit[i-1].compareTo(bSplit[i-1]);
					return reverse ? -res : res;
				}
			}};
	}
	void debugPrint(String c1, String c2)
	{
		System.out.println("Alignments including " + c1 + " and " + c2);
		HashSet<String> c1Reads = alnsToReads(contigToAlignments.get(c1));
		HashSet<String> c2Reads = alnsToReads(contigToAlignments.get(c2));
		HashSet<String> intersection = intersect(c1Reads, c2Reads);
		for(String x : intersection)
		{
			printAlignmentsFromRead(x);
		}
		System.out.println();
	}
}
static class Mapping implements Comparable<Mapping>
{
	String contigName;
	String chrName;
	int chrStart, chrEnd;
	int contigLength;
	boolean reversed;
	Mapping(String nn, String chnn, int ss, int ee, int ll, boolean rr)
	{
		contigName = nn;
		chrName = chnn;
		chrStart = ss;
		chrEnd = ee;
		contigLength = ll;
		reversed = rr;
	}
	@Override
	public int compareTo(Mapping o) {
		// TODO Auto-generated method stub
		if(chrStart != o.chrStart)
		{
			return chrStart - o.chrStart;
		}
		
		if(chrEnd != o.chrEnd)
		{
			return chrEnd - o.chrEnd;
		}
		
		return contigName.compareTo(o.contigName);
	}
}
}
