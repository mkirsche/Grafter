import java.util.*;
import java.io.*;

public class StitchFasta {
static void usage()
{
	System.out.println("java -cp src StitchFasta [contigs] [scaffolds] [outputfasta]");
}
@SuppressWarnings("resource")
public static void main(String[] args) throws IOException
{
	String fastaFn = "", extraFastaFn = "", outFn = "";
	if(args.length >= 3)
	{
		fastaFn = args[0];
		extraFastaFn = args[1];
		outFn = args[2];
	}
	
	else
	{
		usage();
		System.exit(1);
	}
	
	PrintWriter out = new PrintWriter(new File(outFn));
	
	HashSet<String> used = new HashSet<String>();
	
	BufferedReader br = new BufferedReader(new FileReader(new File(extraFastaFn)));
	
	ArrayList<Integer> contigLengths = new ArrayList<Integer>();
	long totLength = 0;
	
	while(true)
	{
		try {
			String[] line = br.readLine().split(" ");
			String name = line[0].substring(1);
			String[] subcontigs = new String[0];
			if(line.length > 1)
			{
				subcontigs = line;
			}
			else
			{
				subcontigs = name.split("&");
			}
			for(String sc : subcontigs)
			{
				used.add(sc);
				
				/*
				 * Handle contigs which were broken in a previous run
				 * This will occur if the contigs are realigned after being broken in misassembly correction
				 */
				if(sc.indexOf('_') != -1)
				{
					used.add(sc.substring(0, sc.lastIndexOf('_')));
				}
			}
			String contig = br.readLine();
			out.println(">" + name);
			contigLengths.add(contig.length());
			totLength += contig.length();
			printContig(out, contig);
		} catch(Exception e) {
			break;
		}
	}
	
	br = new BufferedReader(new FileReader(new File(fastaFn)));
	boolean print = false;
	int curLength = 0;
	while(true)
	{
		try {
			String line = br.readLine();
			if(line.length() == 0) continue;
			if(line.charAt(0) == '>')
			{
				if(print)
				{
					totLength += curLength;
					contigLengths.add(curLength);
				}
				curLength = 0;
				String name = line.substring(1).split(" ")[0];
				if(used.contains(name))
				{
					print = false;
				}
				else
				{
					print = true;
					out.println(line);
				}
			}
			else
			{
				curLength += line.length();
				if(print) out.println(line);
			}
		} catch(Exception e) {
			break;
		}
	}
	out.close();
	
	ReadUtils.assemblyStats(contigLengths, totLength);
}

/*
 * Prints out a contig with 80 characters per line
 */
static void printContig(PrintWriter out, String contig)
{
	int n = contig.length();
	for(int i = 0; i<n; i+= 80)
		out.println(contig.substring(i, Math.min(i+80, n)));
}
}
