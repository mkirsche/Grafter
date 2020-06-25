/*
 * Script for building scaffolds given joins as well as a known layout of pieces
 * Only grafter joins between ends of pieces will be kept
 */

import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;

public class BuildFromLayout {
	static String metadataFn = "";
	static String layoutFn = "";
	static String readsFn = "";
	static String piecesFn = "";
	static String ofn = "";
	static String bedOfn = "";
	static String layoutOfn = "";
	
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
				else if(key.equalsIgnoreCase("layout_file"))
				{
					layoutFn = val;
				}
				else if(key.equalsIgnoreCase("out_file"))
				{
					ofn = val;
				}
				else if(key.equalsIgnoreCase("bed_out"))
				{
					bedOfn = val;
				}
				else if(key.equalsIgnoreCase("layout_out"))
				{
					layoutOfn = val;
				}
				else if(key.equalsIgnoreCase("reads_file"))
				{
					readsFn = val;
				}
				else if(key.equalsIgnoreCase("pieces_file"))
				{
					piecesFn = val;
				}
			}
		}
		
		if(metadataFn.length() == 0 || layoutFn.length() == 0 || ofn.length() == 0 || bedOfn.length() == 0 || readsFn.length() == 0 || piecesFn.length() == 0 || layoutOfn.length() == 0)
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
		System.out.println("Usage: java -cp src BuildFromLayout [args]");
		System.out.println("  Example: java -cp src BuildFromLayout metadata_file=joins.tsv layout_file=layout.txt out_file=out.fasta");
		System.out.println("    reads_file=reads.fastq pieces_file=pieces.fasta out_bed=joins.bed layout_out=layout.bed");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  metadata_file    (String) - a table with joins from Grafter and their supporting reads");
		System.out.println("  layout_file      (String) - a file with the layout of contigs making up each piece");
		System.out.println("  reads_file       (String) - the reads in FASTQ format");
		System.out.println("  pieces_file      (String) - the pieces in FASTA format");
		System.out.println("  out_file         (String) - where to output the scaffolds in FASTA format");
		System.out.println("  bed_out          (String) - where to output the BED file with supporting reads");
		System.out.println("  bed_out          (String) - where to output the BED file with layout of the scaffold");
		System.out.println();
	}
	
	public static void main(String[] args) throws Exception
	{
		parseArgs(args);
		
		buildScaffolds();
	}
	
	/*
	 * Build up scaffolds and produces FASTA as well as BED with supporting reads
	 */
	static void buildScaffolds() throws Exception
	{
		System.err.println("Parsing Grafter joins");
		ArrayList<Join> joins = parseMetadata();
		System.err.println("Found " + joins.size() + " supporting read entries");
		
		System.err.println();
		
		System.err.println("Getting read sequences");
		HashMap<String, String> readSeqs = getRelevantReadSequences(joins);
		System.err.println("Found " + readSeqs.size() + " relevant read sequences");
		
		System.err.println();
		
		System.err.println("Getting piece sequences");
		HashMap<String, String> pieceSeqs = getPieceSeqs();
		System.err.println("Found " + pieceSeqs.size() + " piece sequences");
		
		System.err.println();
		
		System.err.println("Parsing piece layouts");
		ArrayList<Layout> layouts = getLayouts();
		System.err.println("Found layouts for " + layouts.size() + " pieces");
		
		int numPieces = layouts.size();
		
		System.err.println("Building scaffolds");
		
		PrintWriter out = new PrintWriter(new File(ofn));
		PrintWriter bedOut = new PrintWriter(new File(bedOfn));
		PrintWriter layoutOut = new PrintWriter(new File(layoutOfn));
		
		StringBuilder curSeq = new StringBuilder("");
		int scaffolds = 0;
		
		boolean startNewScaffold = true;
		
		// Whether or not the current piece needs to be flipped
		boolean revComp = false;
		
		for(int i = 0; i<numPieces; i++)
		{
			// If starting a new scaffold print the last one and reset variables
			if(startNewScaffold)
			{
				// Print last scaffold
				if(scaffolds != 0)
				{
					out.println(">scaffold" + String.format("%03d", scaffolds));
					out.println(curSeq.toString());
				}
				
				// Move onto the next scaffold and reset sequence
				scaffolds++;
				curSeq = new StringBuilder("");
			}
			
			Layout curLayout = layouts.get(i);

			startNewScaffold = true;
			
			// Try to join with the next sequence
			if(i != numPieces-1)
			{
				Layout nextLayout = layouts.get(i+1);
				
				String lastContig = curLayout.contigNames[curLayout.contigNames.length - 1];
				char lastStrand = curLayout.strands[curLayout.strands.length - 1];
				
				String nextContig = nextLayout.contigNames[0];
				char nextStrand = nextLayout.strands[0];
				
				ArrayList<String> supportingReads = new ArrayList<String>();
				
				String mainRead = "";
				String patchSeq = "";
				
				boolean nextRevComp = false;
				
				boolean readSeqFlipped = false;
				int readStart = -1, readEnd = -1;
								
				for(Join j : joins)
				{
					if(j.contigStart.equals(lastContig) && j.contigEnd.equals(nextContig))
					{
						if(startNewScaffold || j.startPrefix == ((lastStrand == '-') ^ revComp))
						{
							if(i == 0 && j.startPrefix != ((lastStrand == '-') ^ revComp))
							{
								revComp =  true;
							}
							nextRevComp |= j.endPrefix == (nextStrand == '-');
							if(j.sequenceUsed)
							{
								startNewScaffold = false;
								mainRead = j.readName;
								String readSeq = readSeqs.get(j.readName);
								if(j.start > j.end)
								{
									readSeqFlipped = true;
								}
								readStart = Math.min(j.start, j.end);
								readEnd = Math.max(j.start, j.end);
								patchSeq = j.start <= j.end ? readSeq.substring(j.start, j.end) : reverseComplement(readSeq.substring(j.end, j.start));
							}
							else
							{
								supportingReads.add(j.readName);
							}
						}
					}
					
					else if(j.contigEnd.equals(lastContig) && j.contigStart.equals(nextContig))
					{
						if(startNewScaffold || j.endPrefix == ((lastStrand == '-') ^ revComp))
						{
							if(startNewScaffold && j.endPrefix != ((lastStrand == '-') ^ revComp))
							{
								revComp =  true;
							}
							nextRevComp |= j.startPrefix == (nextStrand == '-');
							if(j.sequenceUsed)
							{
								startNewScaffold = false;
								mainRead = j.readName;
								String readSeq = readSeqs.get(j.readName);
								if(j.start > j.end)
								{
									readSeqFlipped = true;
								}
								readStart = Math.min(j.start, j.end);
								readEnd = Math.max(j.start, j.end);
								patchSeq = j.start <= j.end ? reverseComplement(readSeq.substring(j.start, j.end)) : reverseComplement(readSeq.substring(j.end, j.start));
							}
							else
							{
								supportingReads.add(j.readName);
							}
						}
					}
				}
				
				// Add this piece's sequence to the scaffold
				String pieceSeq = pieceSeqs.get(curLayout.pieceName);
				
				String layoutToPrint = curLayout.pieceName + "\t" + 0 + "\t" + pieceSeq.length() + "\t.\t0\t" + (revComp ? '-' : '+');
				
				if(revComp)
				{
					System.err.println("Flipping " + curLayout.pieceName);
				}
				
				curSeq.append(revComp ? reverseComplement(pieceSeq) : pieceSeq);
				revComp = nextRevComp;
				layoutOut.println(layoutToPrint);
				
				if(patchSeq.length() > 0)
				{
					int patchStart = curSeq.length();
					int patchEnd = patchStart + patchSeq.length();
					//int patchLength = patchSeq.length();
					StringBuilder supportingReadField = new StringBuilder(mainRead);
					for(int j = 0; j<supportingReads.size(); j++)
					{
						supportingReadField.append("," + supportingReads.get(j));
					}
					if(supportingReadField.length() == 0)
					{
						supportingReadField = new StringBuilder(".");
					}
					layoutOut.println(mainRead + "\t" + readStart + "\t" + readEnd + "\t.\t0\t" + (readSeqFlipped ? '-' : '+'));
					bedOut.println("scaffold" + String.format("%03d", scaffolds) + "\t" + patchStart + "\t" + patchEnd + "\t" + supportingReadField);
					
					curSeq.append(patchSeq);
				}
				
			}
		}
		
		if(curSeq.length() > 0)
		{
			out.println(">scaffold" + String.format("%03d", scaffolds));
			out.println(curSeq.toString());
		}
		
		out.close();
		bedOut.close();
		layoutOut.close();
		
		System.err.println("Built " + scaffolds + " scaffolds");
	}
	
	/*
	 * Gets the joins and supporting reads from a Grafter metadata table
	 */
	static ArrayList<Join> parseMetadata() throws Exception
	{
		ArrayList<Join> res = new ArrayList<Join>();
		Scanner input = new Scanner(new FileInputStream(new File(metadataFn)));
		String header = input.nextLine();
		while(input.hasNext())
		{
			String line = input.nextLine();
			res.add(new Join(header, line));
		}
		
		return res;
	}
	
	/*
	 * Gets the sequences of reads which are used to fill gaps
	 */
	static HashMap<String, String> getRelevantReadSequences(ArrayList<Join> joins) throws Exception
	{
		HashSet<String> readsNeeded = new HashSet<String>();
		
		for(Join j : joins)
		{
			if(j.sequenceUsed)
			{
				readsNeeded.add(j.readName);
			}
		}
		
		HashMap<String, String> res = new HashMap<String, String>();
		
		Scanner input = new Scanner(new FileInputStream(new File(readsFn)));
		String[] lines = new String[4];
		while(input.hasNext())
		{
			for(int i = 0; i<4; i++)
			{
				lines[i] = input.nextLine();
			}
			String readName = lines[0].split(" ")[0].substring(1).trim();
			
			if(readsNeeded.contains(readName))
			{
				res.put(readName, lines[1]);
			}
		}
		
		return res;
	}
	
	/*
	 * Gets the layout of each piece
	 */
	static ArrayList<Layout> getLayouts() throws Exception
	{
		ArrayList<Layout> res = new ArrayList<Layout>();
		Scanner input = new Scanner(new FileInputStream(new File(layoutFn)));
		while(input.hasNext())
		{
			String line = input.nextLine().trim();
			if(line.length() == 0)
			{
				continue;
			}
			res.add(new Layout(line));
		}
		return res;
	}
	
	/*
	 * Gets a map from piece name to its sequence
	 */
	static HashMap<String, String> getPieceSeqs() throws Exception
	{
		HashMap<String, String> res = new HashMap<String, String>();
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
					res.put(pieceName, sb.toString());
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
			res.put(pieceName, sb.toString());
		}
		return res;
	}

	
	static class Join
	{
		String readName;
		int start, end;
		String contigStart, contigEnd;
		boolean startPrefix, endPrefix;
		boolean sequenceUsed;
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
			}
		}
	}
	
	static class Layout
	{
		String pieceName;
		String[] contigNames;
		char[] strands;
		Layout(String line)
		{
			int spaceIndex = line.indexOf(' ');
			pieceName = line.substring(0, spaceIndex);
			String[] tokens = line.substring(1 + spaceIndex).trim().split(",");
			int n = tokens.length;
			contigNames = new String[n];
			strands = new char[n];
			for(int i = 0; i<n; i++)
			{
				contigNames[i] = tokens[i].substring(0, tokens[i].length() - 1);
				strands[i] = tokens[i].charAt(tokens[i].length() - 1);
			}
		}
	}
	
	/*
	 * Computes the reverse complement of a string
	 */
	static String reverseComplement(String s)
	{
		int n = s.length();
		char[] res = new char[n];
		for(int i = 0; i<n; i++)
		{
			char c = s.charAt(n - 1 - i);
			res[i] = c;
			if(c == 'A') res[i] = 'T';
			if(c == 'C') res[i] = 'G';
			if(c == 'G') res[i] = 'C';
			if(c == 'T') res[i] = 'A';
			if(c == 'a') res[i] = 't';
			if(c == 'c') res[i] = 'g';
			if(c == 'g') res[i] = 'c';
			if(c == 't') res[i] = 'a';
		}
		return new String(res);
	}
}
