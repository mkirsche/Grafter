import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

public class AlignmentGatherer {
	
	/*
	 * Compresses the alignments to a given read by combining alignments of the same contig into one
	 * Also, filters out invalid alignments
	 */
	static ArrayList<SortablePafAlignment> compress(ArrayList<SortablePafAlignment> alignments, boolean filterInvalid)
	{
		int n = alignments.size();

		Comparator<SortablePafAlignment> byContigName = new Comparator<SortablePafAlignment>() {

			@Override
			public int compare(SortablePafAlignment a, SortablePafAlignment b) {
				if(a.contigName.equals(b.contigName))
				{
					return a.compareTo(b);
				}
				return a.contigName.compareTo(b.contigName);
			}
		};
		
		// Sort by contig name and break ties by read start position
		Collections.sort(alignments, byContigName);
		ArrayList<SortablePafAlignment> filtered = new ArrayList<>();
		for(int i = 0; i<n; i++)
		{
			// Find the end of the run of alignments of the current contig
			int j = i+1;
			while(j< n && alignments.get(i).contigName.equals(alignments.get(j).contigName))
			{
				j++;
			}
			
			// Now alignments[i:j) has all the alignments of this contig - combine or remove them
			boolean[] rse = new boolean[2], cse = new boolean[2];
			boolean gapFree = true;
			int lastReadEndPosition = alignments.get(i).readEnd;
			int lastContigEndPosition = alignments.get(i).contigEnd;
			for(int k = i; k<j; k++)
			{
				SortablePafAlignment cur = alignments.get(k);
				
				if(k > i && alignments.get(k-1).strand != alignments.get(k).strand)
				{
					j = k;
					break;
				}
				
				// Check for a gap between this alignment and the last one in either the read or contig
				if(cur.contigStart - lastContigEndPosition > Settings.MAX_GAP)
				{
					gapFree = false;
					break;
				}
				if(cur.readStart - lastReadEndPosition > Settings.MAX_GAP)
				{
					gapFree = false;
					break;
				}
				
				lastContigEndPosition = cur.contigEnd;
				lastReadEndPosition = cur.readEnd;
				boolean[] crse = readStartEnd(cur);
				boolean[] ccse = contigStartEnd(cur);
				for(int q = 0; q<2; q++)
				{
					rse[q] |= crse[q];
					cse[q] |= ccse[q];
				}
			}
			
			/*
			 * If the set of alignments had a large gap, ignore it 
			 */
			if(!gapFree)
			{
				i = j - 1;
				continue;
			}
			
			/*
			 * We have whether the alignment set covers the start/end of contig/read, so check that it's valid
			 */
			if(filterInvalid && !cse[0] && !cse[1])
			{
				/*
				 * Middle portion of contig aligns somewhere on read but neither end of it
				 * Only possible case is read contained in contig, making alignment useless
				 * Throw out the alignment and update i
				 */
				i = j - 1;
				continue;
			}
			else if(filterInvalid && !rse[0] && !rse[1])
			{
				/*
				 * Neither end of the read is involved, so contig must be contained in the read
				 * Filter out cases which don't reflect this
				 */
				if(!cse[0] || !cse[1])
				{
					i = j - 1;
					continue;
				}
			}
			else
			{
				/*
				 * Create consensus of all of the alignments by taking earliest start and latest end
				 */
				SortablePafAlignment total = alignments.get(i).copy();
				for(int k = i+1; k<j; k++)
				{
					SortablePafAlignment cur = alignments.get(k);
					total.contigStart = Math.min(total.contigStart, cur.contigStart);
					total.contigEnd = Math.max(total.contigEnd, cur.contigEnd);
					total.readStart = Math.min(total.readStart, cur.readStart);
					total.readEnd = Math.max(total.readEnd, cur.readEnd);
				}
				i = j - 1;
				filtered.add(total);
			}
		}
		
		Collections.sort(filtered);
		if(filtered.size() > 0)
		{
		System.err.println("Chain " + alignments.get(0).readName);
		for(SortablePafAlignment spa : filtered)
		{
			System.err.println("  " + spa.readStart + " " + spa.readEnd + " " + spa.contigName + " " + spa.contigStart +" "+spa.contigEnd);
		}
		}
		return filtered;
	}

	/*
	 * Gets chains of unique matches to a read given the list of all of the alignments to it
	 */
	@SuppressWarnings("unchecked")
	static ArrayList<ArrayList<SortablePafAlignment>> getUniqueMatches(ArrayList<SortablePafAlignment> alignments)
	{
		/*
		 * Sort by start point
		 */
		Collections.sort(alignments);
		
		/*
		 * Compress all alignments of the same contig and remove invalid alignments
		 */
		alignments = compress(alignments, true);
		
		/*
		 * List of chains of alignments
		 */
		ArrayList<ArrayList<SortablePafAlignment>> res = new ArrayList<>();
		
		/*
		 * The list of alignments in the current chain
		 */
		ArrayList<SortablePafAlignment> cur = new ArrayList<>();
		for(int i = 0 ; i<alignments.size(); i++)
		{
			
			SortablePafAlignment a = alignments.get(i);
			
			/*
			 * Cases: 
			 *   1.) Contained in a previous alignment -> Ignore this alignment
			 *   2.) Overlaps last two alignments -> End chain here
			 *   3.) Valid continuation of chain
			 */
			if(cur.size() >= 1 && cur.get(cur.size()-1).readEnd >= a.readEnd)
			{
				// Contained in a previous alignment
				continue;
			}
			else if(cur.size() >= 2 && cur.get(cur.size() - 2).readEnd > a.readStart)
			{
				// Overlaps last two alignments
				cur.remove(cur.size()-1);
				if(cur.size() >= 2)
				{
					res.add((ArrayList<SortablePafAlignment>) cur.clone());
				}
				cur.clear();
			}
			else if(cur.size() >= 1 && cur.get(cur.size() - 1).readEnd + 100000 < a.readStart)
			{
				if(cur.size() >= 2)
				{
					res.add((ArrayList<SortablePafAlignment>) cur.clone());
				}
				cur.clear();
				cur.add(a);
			}
			else
			{
				// Valid continuation of chain
				cur.add(a);
			}
		}
		
		// Add leftover chain
		if(cur.size() >= 2)
		{
			res.add(cur);
		}
		
		return res;
	}
	
	/*
	 * General end checking for alignments
	 */
	static boolean[] startEnd(int startPos, int endPos, int length)
	{
		double curMaxHanging = Math.min(Settings.MAX_HANGING_PROP*length, Settings.MAX_HANGING);
		return new boolean[] {startPos < curMaxHanging, endPos + curMaxHanging >= length};
	}

	/*
	 * Whether or not an alignment contains the start/end of the contig involved
	 */
	static boolean[] contigStartEnd(SortablePafAlignment pa)
	{
		return startEnd(pa.contigStart, pa.contigEnd, pa.contigLength);
	}

	/*
	 * Whether or not an alignment contains the start/end of the read involved
	 */
	static boolean[] readStartEnd(SortablePafAlignment pa)
	{
		return startEnd(pa.readStart, pa.readEnd, pa.readLength);
	}
}
