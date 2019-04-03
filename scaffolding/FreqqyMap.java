package scaffolding;

import java.util.*;

public class FreqqyMap {

	// Map contig name to cumulative sum array of its kmer frequencies
	HashMap<String, long[]> contigToFreqSum;
	
	// Map kmer code to its frequency
	int[] kmerFrequencies;
	
	HashMap<String, Integer> contigLengths;
	
	// Length of kmers to use (default 13)
	int k;
	static int defaultK = 13;
	
	int samplingFrequency = 100;
	
	int totalLength = 0;
	
	FreqqyMap(HashMap<String, String> seqMap, int k)
	{
		this.k = k;
		kmerFrequencies = new int[1<<(2*k)];
		contigToFreqSum = new HashMap<>();
		contigLengths = new HashMap<>();
		countKmers(seqMap);
		buildSumArrays(seqMap);
	}
	
	FreqqyMap(HashMap<String, String> seqMap)
	{
		this(seqMap, defaultK);
	}
	
	FreqqyMap()
	{
		k = defaultK;
		kmerFrequencies = new int[1<<(2*k)];
		contigToFreqSum = new HashMap<>();
		contigLengths = new HashMap<>();
	}
	
	double getAverageFrequency(String name, int start, int end)
	{
		System.out.println(name+" "+contigLengths.get(name)+" "+start+" "+end);
		long[] csum = contigToFreqSum.get(name);
		int a = Math.max(0, start - k + 1);
		int b = Math.min(contigLengths.get(name)-k-1, end);
		System.out.println(a+" "+b);
		double totalFreq = csumQuery(csum, b) - csumQuery(csum, a-1);
		return totalFreq / (b - a + 1);
	}
	
	void addSumArray(String key, String s)
	{
		int[] f = new int[Math.max(0, s.length()-k+1)];
		int kmer = 0;
		for(int i = 0; i<s.length(); i++)
		{
			kmer <<= 2;
			kmer &= (1 << (2*k)) - 1;
			kmer |= charToInt(s.charAt(i));
			if(i >= k-1)
			{
				f[i-k+1] = kmerFrequencies[kmer];
			}
		}
		long[] sampleF = sample(f);
		for(int i = 1; i<sampleF.length; i++)
		{
			sampleF[i] += sampleF[i-1];
		}
		contigToFreqSum.put(key, sampleF);
	}
	
	void buildSumArrays(HashMap<String, String> seqMap)
	{
		for(String key : seqMap.keySet())
		{
			addSumArray(key, seqMap.get(key));
		}
	}
	
	long[] sample(int[] f)
	{
		long[] res = new long[(f.length + samplingFrequency - 1) / samplingFrequency];
		for(int i = 0; i<f.length; i+=samplingFrequency)
		{
			for(int j = i; j<i+samplingFrequency && j<f.length; j++)
			{
				res[i/samplingFrequency] += f[j];
			}
		}
		return res;
	}
	
	double csumQuery(long[] sample, int a)
	{
		if(a < 0) return 0;
		int entireBlocks = (a + 1) / samplingFrequency;
		double res = entireBlocks == 0 ? 0 : sample[entireBlocks-1];
		int leftover = (a+1) % samplingFrequency;
		res += (sample[entireBlocks] - (entireBlocks == 0 ? 0 : sample[entireBlocks-1])) * 1. * leftover / samplingFrequency;
		return res;
	}
	
	void addKmerCount(String key, String s)
	{
		totalLength += s.length();
		int kmer = 0;
		for(int i = 0; i<s.length(); i++)
		{
			kmer <<= 2;
			kmer &= (1 << (2*k)) - 1;
			kmer |= charToInt(s.charAt(i));
			if(i >= k-1)
			{
				addKmer(kmer);
			}
		}
		contigLengths.put(key, s.length());
	}
	
	void countKmers(HashMap<String, String> seqMap)
	{
		for(String key : seqMap.keySet())
		{
			String s =  seqMap.get(key);
			addKmerCount(key, s);
		}
	}
	
	static int charToInt(char c)
	{
		if(c == 'A' || c == 'a') return 0;
		else if(c == 'C' || c == 'c') return 1;
		else if(c == 'G' || c == 'g') return 2;
		return 3;
	}
	
	void addKmer(int kmer)
	{
		kmerFrequencies[kmer]++;
	}
}
