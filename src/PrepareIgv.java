import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Scanner;

public class PrepareIgv {
	static String agpFile = "", scaffoldFile = "", tdfFile = "", ontTdfFile = "";
	static String contigFile = "";
	static String outPrefix = "igv";
	
	/*
	 * Parses command line arguments
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
				if(key.equalsIgnoreCase("agp_file"))
				{
					agpFile = val;
				}
				else if(key.equalsIgnoreCase("scaffold_file"))
				{
					scaffoldFile = val;
				}
				else if(key.equalsIgnoreCase("ill_tdf_file"))
				{
					tdfFile = val;
				}
				else if(key.equalsIgnoreCase("ont_tdf_file"))
				{
					ontTdfFile = val;
				}
				else if(key.equalsIgnoreCase("contig_file"))
				{
					contigFile = val;
				}
				else if(key.equalsIgnoreCase("out_prefix"))
				{
					outPrefix = val;
				}
			}
		}
		
		if(agpFile.length() == 0 || scaffoldFile.length() == 0 || contigFile.length() == 0)
		{
			usage();
			System.exit(1);
		}
	}
	
	/*
	 * Prints a usage menu
	 */
	static void usage()
	{
		System.out.println();
		System.out.println("Usage: java -cp src PrepareIgv [args]");
		System.out.println("  Example: java -cp src PrepareIgv agp_file=layout.agp scaffold_file=scaffolds.fa.txt");
		System.out.println("    contig_file=contigs.fa ill_tdf_file=ill.tdf ont_tdf_file=ont.tdf");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  agp_file      (String)       - the AGP format file listing scaffold layouts");
		System.out.println("  scaffold_file (String)       - the sequences of scaffoldsto which Illumina and ONT data were aligned");
		System.out.println("  contig_file   (String)       - the sequences of the contigs which were scaffolded by Grafter");
		System.out.println();
		System.out.println("Optional args");
		System.out.println("  ill_tdf_file  (String)       - the coverage of Illumina data across the scaffolds");
		System.out.println("  ont_tdf_file  (String)       - the coverage of ONT data across the scaffolds");
		System.out.println("  out_prefix    (String) [igv] - where the IGV screenshots should go");
		System.out.println();
	}
	
	
public static void main(String[] args) throws Exception
{
	parseArgs(args);
	contigFile = "/home/mkirsche/physalis/purged.fa";
	agpFile = "/home/mkirsche/eclipse-workspace/Grafter/test_ragtag.agp";
	scaffoldFile = "/home/mkirsche/physalis/scaffolds_ragtag.fa";
	tdfFile = "/home/mkirsche/physalis/M82_grafter_self.sr.tdf";
	ontTdfFile = "/home/mkirsche/physalis/M82_grafter_self.ont.tdf";
	
	HashMap<String, Integer> contigLengths = getLengths();
	
	HashMap<String, String> contigToScaffold = new HashMap<String, String>();
	Scanner input = new Scanner(new FileInputStream(new File(scaffoldFile + ".fai")));
	while(input.hasNext())
	{
		String line = input.nextLine();
		String scaffold = line.split("\t")[0];
		String[] contigs = scaffold.split("&");
		for(String c : contigs)
		{
			contigToScaffold.put(c, scaffold);
		}
	}
	input.close();
	
	HashMap<String, ArrayList<Position>> positions = getPositions(contigLengths, contigToScaffold);
	
	Path currentRelativePath = Paths.get("");
	String outDir = currentRelativePath.toAbsolutePath().toString() + "/" + outPrefix;
	File outDirFile = new File(outDir);
	if(outDirFile.isDirectory())
	{
		final File[] files = outDirFile.listFiles();
		for (File f: files) f.delete();
		outDirFile.delete();
	}
	outDirFile.mkdir();
	String ofn = outDir + "/" + outPrefix + ".bat";
	PrintWriter out = new PrintWriter(new File(ofn));
	out.println("new");
	out.println("genome " + (scaffoldFile.startsWith("/") ? 
			scaffoldFile : (currentRelativePath.toAbsolutePath().toString() + "/" + scaffoldFile)));
	out.println("load " + (tdfFile.startsWith("/") ? 
			tdfFile : (currentRelativePath.toAbsolutePath().toString() + "/" + tdfFile)));
	out.println("load " + (ontTdfFile.startsWith("/") ? 
			ontTdfFile : (currentRelativePath.toAbsolutePath().toString() + "/" + ontTdfFile)));
	
	PrintWriter bedOut = new PrintWriter(new File(outDir + "/" + "out.bed"));
	PrintWriter bedReadOut = new PrintWriter(outDir + "/" + new File("outreads.bed"));
	//PrintWriter bedOverlapOut = new PrintWriter(outDir + "/" + new File("outoverlaps.bed"));
	
	for(String s : positions.keySet())
	{
		//System.out.println(s);
		for(int i = 0; i < positions.get(s).size(); i++)
		{
			Position p = positions.get(s).get(i);
			
			if(contigLengths.containsKey(p.contig))
			{
				bedOut.println(s+"\t"+p.start+"\t"+p.end+"\t"+p.contig);
			}
			else
			{
				bedReadOut.println(s+"\t"+p.start+"\t"+p.end+"\t"+p.contig);
			}
			
			//if(i < positions.get(s).size() - 1 && positions.get(s).get(i+1).contig.startsWith("tig"))
			//{
			//	if(positions.get(s).get(i).end > positions.get(s).get(i+1).start)
			//}
		}
	}
	bedOut.close();
	bedReadOut.close();
	//bedOverlapOut.close();
	
	out.println("load " + (outDir + "/out.bed"));
	out.println("load " + (outDir + "/outreads.bed"));
	//out.println("load " + (outDir + "/outoverlaps.bed"));

	
	out.println("snapshotDirectory " + outDir);

	
	for(String s : positions.keySet())
	{
		//System.out.println(s);
		for(int i = 0; i < positions.get(s).size()-1; i++)
		{
			Position p = positions.get(s).get(i);
			if(!p.contig.startsWith("tig"))
			{
				continue;
			}
			int nexti = i+1;
			while(nexti < positions.get(s).size() && !positions.get(s).get(nexti).contig.startsWith("tig"))
			{
				nexti++;
			}
			if(nexti == positions.get(s).size())
			{
				break;
			}
			Position next = positions.get(s).get(nexti);
			out.println("goto " + s + ":" + (p.end - 5000) + "-" + (next.start + 5000));
			out.println("snapshot " + p.contig + "_joinedwith_" + next.contig + (".png"));
			//System.out.println(p.contig+" "+p.start+" "+p.end);
		}
	}
	
	out.println("exit");
	
	out.close();
}

static HashMap<String, ArrayList<Position>> getPositions(HashMap<String, Integer> lengths, HashMap<String, String> contigToScaffold) throws Exception
{
	Scanner input = new Scanner(new FileInputStream(new File(agpFile)));
	
	// For each scaffold in the AGP file, map its name to a list of which contigs are in it and which intervals they make up
	HashMap<String, ArrayList<Position>> componentList = new HashMap<String, ArrayList<Position>>();
	while(input.hasNext())
	{
		String line = input.nextLine();
		String[] tokens = line.split("\t");
		String scaffoldName = tokens[0];
		String contigName = tokens[5];
		if(!contigToScaffold.containsKey(contigName))
		{
			//continue;
		}
		if(!componentList.containsKey(scaffoldName))
		{
			componentList.put(scaffoldName, new ArrayList<Position>());
		}
		int scaffoldStart = Integer.parseInt(tokens[1]), scaffoldEnd = Integer.parseInt(tokens[2]);
		componentList.get(scaffoldName).add(new Position(contigName, scaffoldStart, scaffoldEnd));
	}
	
	// Integrate with the FASTA scaffolds
	HashMap<String, ArrayList<Position>> res = new HashMap<String, ArrayList<Position>>();
	for(String s : componentList.keySet())
	{
		ArrayList<Position> pieces = componentList.get(s);
		if(pieces.size() == 1)
		{
			continue;
		}
				
		String scaffoldName = contigToScaffold.get(pieces.get(0).contig);
		
		res.put(scaffoldName, new ArrayList<Position>());

		int length = pieces.get(pieces.size() - 1).end;
		
		if(scaffoldName.startsWith(pieces.get(0).contig))
		{
			for(Position p : pieces)
			{
				res.get(scaffoldName).add(p);
			}
		}
		else
		{
			for(Position p : pieces)
			{
				res.get(scaffoldName).add(new Position(p.contig, length - p.end, length - p.start));
			}
			Collections.reverse(res.get(scaffoldName));
		}
	}
	
	return res;
}

static class Piece
{
	String line;
	Piece(String line)
	{
		this.line = line;
		String[] tokens = line.split("\t");
	}
}

static class Position
{
	String contig;
	int start, end;
	Position(String c, int s, int e)
	{
		contig = c;
		start = s;
		end = e;
	}
}

static HashMap<String, Integer> getLengths() throws Exception
{
	HashMap<String, Integer> res = new HashMap<String, Integer>();
	Scanner input = new Scanner(new FileInputStream(new File(contigFile + ".fai")));
	while(input.hasNext())
	{
		String line = input.nextLine();
		String[] tokens = line.split("\t");
		res.put(tokens[0], Integer.parseInt(tokens[1]));
	}
	return res;
}
}
