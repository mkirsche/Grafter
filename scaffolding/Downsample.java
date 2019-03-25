package scaffolding;

import java.util.*;
import java.io.*;

public class Downsample {
@SuppressWarnings("resource")
public static void main(String[] args) throws IOException
{
	String fn = args[0];
	double proportion = Double.parseDouble(args[1]);
	PrintWriter out = new PrintWriter(System.out);
	BufferedReader br = new BufferedReader(new FileReader(new File(fn)));
	Random r = new Random(987654321);
	while(true)
	{
		try {
			boolean keep = r.nextDouble() < proportion;
			if(!keep)
			{
				for(int i = 0; i<4; i++)
				{
					br.readLine();
				}
			}
			else
			{
				for(int i = 0; i<4; i++)
				{
					out.println(br.readLine());
				}
			}
		} catch(Exception e) {
			break;
		}
	}
	out.close();
}
}
