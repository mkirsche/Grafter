import java.util.*;

import java.io.*;
public class TigsFromLayout {
	
	static String layoutFn, ofn;
public static void main(String[] args) throws Exception
{
	if(args.length != 2)
	{
		System.out.println("Usage: java TigsFromLayout <layoutfile> <outfile>");
		System.exit(1);
	}
	layoutFn = args[0];
	ofn = args[1];
	ArrayList<Layout> layouts = getLayouts();
	PrintWriter out = new PrintWriter(new File(ofn));
	for(Layout lay : layouts)
	{
		for(String s : lay.contigNames)
		{
			out.println(s);
		}
	}
	out.close();
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
}
