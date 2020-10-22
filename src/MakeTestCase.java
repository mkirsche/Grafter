import java.io.File;
import java.io.PrintWriter;
import java.util.Random;

public class MakeTestCase {
	public static void main(String[] args) throws Exception
	{
		makeTest(1, 100000, 80000, 70000, 60000, 5000);
		makeTest(2, 100000, 80000, 70000, 60000, 5000);
		makeOverlapTest(3, 100000, 200000, 70000, 60000, 1000);
	}
	
	static void makeTest(int testNum, int contigLen1, int contigLen2, int readC1, int readC2, int readOnlyLen) throws Exception
	{
		String contigsFn = "testcontigs" + testNum + ".fasta";
		String readsFn = "testreads" + testNum + ".fasta";
		String alnFn = "testaln" + testNum + ".paf";
		String truthFn = "testtruth" + testNum + ".fasta";

		String ctg1 = randomSeq(contigLen1);
		String ctg2 = randomSeq(contigLen2);
		String ctg3 = randomSeq(1000);
		String readOnly = randomSeq(readOnlyLen);
		
		PrintWriter out = new PrintWriter(new File(contigsFn));
		out.println(">contig1");
		out.println(ctg1);
		out.println(">contig2");
		out.println(ctg2);
		out.println(">contig3");
		out.println(ctg3);
		out.close();
		
		String readSeq = ctg1.substring(ctg1.length() - readC1) + readOnly + ctg2.substring(0, readC2);
		out = new PrintWriter(new File(readsFn));
		out.println(">read1");
		if(testNum%2 == 1)
		{
			out.println(readSeq);
		}
		else
		{
			out.println(ReadUtils.reverseComplement(readSeq));
		}
		out.close();
		
		String genome = ctg1 + readOnly + ctg2;
		out = new PrintWriter(new File(truthFn));
		out.println(">scaffold1");
		out.println(genome);
		out.close();
		
		out = new PrintWriter(new File(alnFn));
		
		if(testNum == 1)
		{
			out.printf("%s\t%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", 
					"read1",
					readSeq.length(),
					0,
					readC1,
					'+',
					"contig1",
					contigLen1,
					ctg1.length() - readC1,
					ctg1.length(),
					0,
					0,
					60
			);
			out.printf("%s\t%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", 
					"read1",
					readSeq.length(),
					readSeq.length() - readC2,
					readSeq.length(),
					'+',
					"contig2",
					contigLen2,
					0,
					readC2,
					0,
					0,
					60
			);
		}
		else if(testNum == 2)
		{
			out.printf("%s\t%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", 
					"read1",
					readSeq.length(),
					0,
					readC2,
					'-',
					"contig2",
					contigLen2,
					0,
					readC2,
					0,
					0,
					60
			);
			out.printf("%s\t%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", 
					"read1",
					readSeq.length(),
					readSeq.length() - readC1,
					readSeq.length(),
					'-',
					"contig1",
					contigLen1,
					contigLen1 - readC1,
					contigLen1,
					0,
					0,
					60
			);	
		}
		out.close();
		
	}
	
	static void makeOverlapTest(int testNum, int contigLen1, int contigLen2, int readC1, int readC2, int overlap) throws Exception
	{
		String contigsFn = "testcontigs" + testNum + ".fasta";
		String readsFn = "testreads" + testNum + ".fasta";
		String alnFn = "testaln" + testNum + ".paf";
		String truthFn = "testtruth" + testNum + ".fasta";

		String overlapSeq = randomSeq(overlap);
		String overlapSeq2 = randomSeq(overlap);
		String ctg1 = randomSeq(contigLen1 - overlap) + overlapSeq;
		String ctg2 = overlapSeq + randomSeq(contigLen2 - 2*overlap) + overlapSeq2;
		String ctg3 = overlapSeq2 + randomSeq(contigLen1 - overlap);
		
		PrintWriter out = new PrintWriter(new File(contigsFn));
		out.println(">contig1");
		out.println(ctg1);
		out.println(">contig2");
		out.println(ctg2);
		out.println(">contig3");
		out.println(ctg3);
		out.close();
		
		String readSeq = ctg1.substring(ctg1.length() - readC1) + ctg2.substring(overlap) + ctg3.substring(overlap, readC2);
		out = new PrintWriter(new File(readsFn));
		out.println(">read1");
		out.println(readSeq);
		out.close();
		
		String genome = ctg1 + ctg2.substring(overlap) + ctg3.substring(overlap);
		out = new PrintWriter(new File(truthFn));
		out.println(">scaffold1");
		out.println(genome);
		out.close();
		
		out = new PrintWriter(new File(alnFn));

		out.printf("%s\t%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", 
				"read1",
				readSeq.length(),
				0,
				readC1,
				'+',
				"contig1",
				contigLen1,
				ctg1.length() - readC1,
				ctg1.length(),
				0,
				0,
				60
		);
		out.printf("%s\t%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", 
				"read1",
				readSeq.length(),
				readC1 - overlap,
				readSeq.length() - readC2 + overlap,
				'+',
				"contig2",
				contigLen2,
				0,
				contigLen2,
				0,
				0,
				60
		);
		out.printf("%s\t%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", 
				"read1",
				readSeq.length(),
				readSeq.length() - readC2,
				readSeq.length(),
				'+',
				"contig3",
				contigLen1,
				0,
				readC2,
				0,
				0,
				60
		);
		out.close();
		
	}
		
	static String randomSeq(int len)
	{
		Random r = new Random();
		char[] res = new char[len];
		for(int i = 0; i<len; i++)
		{
			int x = r.nextInt(4);
			if(x == 0) res[i] = 'A';
			else if(x == 1) res[i] = 'C';
			else if(x == 2) res[i] = 'G';
			else if(x == 3) res[i] = 'T';
		}
		return new String(res);
	}
}
