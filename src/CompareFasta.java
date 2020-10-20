import java.io.File;
import java.io.FileInputStream;
import java.util.Scanner;

public class CompareFasta {
public static void main(String[] args) throws Exception
{
	String fn1 = args[0], fn2 = args[1];
	Scanner input = new Scanner(new FileInputStream(new File(fn1)));
	input.nextLine();
	String seq1 = input.nextLine();
	input.close();
	input = new Scanner(new FileInputStream(new File(fn2)));
	input.nextLine();
	String seq2 = input.nextLine();
	int res1 = compare(seq1,seq2);
	int res2 = compare(seq1, rc(seq2));
	System.out.println(res1+" "+res2);
	int res = Math.min(res1, res2);
	System.out.println(res);
	input.close();
}
static String rc(String s)
{
	int n = s.length();
	char[] res = new char[n];
	for(int i = 0; i<n; i++)
	{
		char c = s.charAt(i);
		char rc = 'A';
		if(c == 'a' || c == 'A') rc = 'T';
		if(c == 'c' || c == 'C') rc = 'G';
		if(c == 'g' || c == 'G') rc = 'C';
		if(c == 't' || c == 'T') rc = 'A';
		res[n-1-i] = rc;
	}
	return new String(res);
}
static int compare(String a, String b)
{
	int n = a.length(), m = b.length();
	int res = Math.abs(n - m);
	for(int i = 0; i<Math.min(n, m); i++)
	{
		if(a.charAt(i) != b.charAt(i))
		{
			res++;
		}
	}
	return res;
}
}
