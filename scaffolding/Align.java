/*
 * TODO fill in aligner with dynamic programming
 */

package scaffolding;

public class Align {

	static String joinContigs(String lastSeq, String curSeq, int overlap, boolean lastReversed, boolean curReversed)
	{
		String s = lastSeq;
		if(lastReversed)
		{
			s = Scaffold.reverseComplement(s);
		}
		String t = curSeq;
		if(curReversed)
		{
			t = Scaffold.reverseComplement(t);
		}
		
		return t.substring(overlap);
	}
	
	static String joinReadContigs(String lastSeq, String curSeq, int overlap, boolean lastReversed, boolean curReversed, String readSeq, boolean readReversed)
	{
		return "";
	}
}
