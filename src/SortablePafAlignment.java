/*
 * Alignment of a contig to an ultralong read - sortable by start position in the read 
 */
public class SortablePafAlignment implements Comparable<SortablePafAlignment>
{
	String readName, contigName;
	int readLength, readStart, readEnd;
	int contigLength, contigStart, contigEnd;
	int mapq;
	char strand;
	String line;
	// Call with a second parameter to denote that read and contig were flipped when calling minimap2
	SortablePafAlignment(String line, int backwards)
	{
		this.line = line;
		String[] ss = line.split("\t");
		contigName = ss[0];
		contigLength = Integer.parseInt(ss[1]);
		contigStart = Integer.parseInt(ss[2]);
		contigEnd = Integer.parseInt(ss[3]);
		strand = ss[4].charAt(0);
		readName = ss[5];
		readLength = Integer.parseInt(ss[6]);
		readStart = Integer.parseInt(ss[7]);
		readEnd = Integer.parseInt(ss[8]);
		mapq = Integer.parseInt(ss[11]);
	}
	SortablePafAlignment(String line)
	{
		this.line = line;
		String[] ss = line.split("\t");
		readName = ss[0];
		readLength = Integer.parseInt(ss[1]);
		readStart = Integer.parseInt(ss[2]);
		readEnd = Integer.parseInt(ss[3]);
		strand = ss[4].charAt(0);
		contigName = ss[5];
		contigLength = Integer.parseInt(ss[6]);
		contigStart = Integer.parseInt(ss[7]);
		contigEnd = Integer.parseInt(ss[8]);
		mapq = Integer.parseInt(ss[11]);
	}

	public int compareTo(SortablePafAlignment o) {
		if(readStart != o.readStart)
		{
			return readStart - o.readStart;
		}
		return readEnd - o.readEnd;
	}
	SortablePafAlignment copy()
	{
		SortablePafAlignment res = new SortablePafAlignment(line);
		if(!res.contigName.equals(contigName)) res.contigName = contigName;
		if(res.contigStart != contigStart) res.contigStart = contigStart;
		if(res.contigLength != contigLength) res.contigLength = contigLength;
		if(res.contigEnd != contigEnd) res.contigEnd = contigEnd;
		return res;
	}
	
}