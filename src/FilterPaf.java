import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Scanner;

public class FilterPaf {
public static void main(String[] args) throws Exception
{
	String pafFn = args[0];
	String fastaFn = args[1];
	String ofn = args[2];
	
	HashSet<String> contigNames = new HashSet<String>();
	Scanner input = new Scanner(new FileInputStream(new File(fastaFn)));
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.startsWith(">"))
		{
			String name = line.split(" ")[0].substring(1);
			contigNames.add(name);
		}
	}
	input.close();
	
	input = new Scanner(new FileInputStream(new File(pafFn)));
	PrintWriter out = new PrintWriter(new File(ofn));
	while(input.hasNext())
	{
		String line = input.nextLine();
		String[] tokens = line.split("\t");
		String contigName = tokens[5];
		if(contigNames.contains(contigName))
		{
			out.println(line);
		}
	}
	input.close();
	out.close();
}
}
