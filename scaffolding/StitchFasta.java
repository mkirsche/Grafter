package scaffolding;

import java.util.*;
import java.io.*;

public class StitchFasta {
public static void main(String[] args) throws IOException
{
	String fastaFn = "/scratch/groups/mschatz1/mkirsche/ultralong/ccs/paternal_and_unknown.contigs.mmpoa.fa";
	String extraFastaFn = "/scratch/groups/mschatz1/mkirsche/ultralong/ccs/paternal_newcontigs.fa";
	String outfn = "/scratch/groups/mschatz1/mkirsche/ultralong/ccs/paternal_and_unknown.contigs.mmpoa.scaffolds.fa";
	
	if(args.length >= 3)
	{
		fastaFn = args[0];
		extraFastaFn = args[1];
		outfn = args[2];
	}
	
	PrintWriter out = new PrintWriter(new File(outfn));
	
	HashSet<String> used = new HashSet<String>();
	
	BufferedReader br = new BufferedReader(new FileReader(new File(extraFastaFn)));
	
	while(true)
	{
		try {
			String line = br.readLine();
			String name = line.substring(1);
			String[] subcontigs = name.split("&");
			for(String sc : subcontigs)
			{
				used.add(sc);
			}
			String contig = br.readLine();
			out.println(">" + name);
			printContig(out, contig);
		} catch(Exception e) {
			break;
		}
	}
	
	br = new BufferedReader(new FileReader(new File(fastaFn)));
	boolean print = false;
	while(true)
	{
		try {
			String line = br.readLine();
			if(line.length() == 0) continue;
			if(line.charAt(0) == '>')
			{
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
				if(print) out.println(line);
			}
		} catch(Exception e) {
			break;
		}
	}
	out.close();
}
static void printContig(PrintWriter out, String contig)
{
	int n = contig.length();
	for(int i = 0; i<n; i+= 80)
		out.println(contig.substring(i, Math.min(i+80, n)));
}
}
