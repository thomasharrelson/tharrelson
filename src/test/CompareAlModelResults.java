package test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.TreeSet;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;

import org.apache.log4j.Logger;

import broad.core.parser.CommandLineParser;

public class CompareAlModelResults {
	private static Logger logger = Logger.getLogger(CompareAlModelResults.class.getName());
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		ArrayList<Double> diff = new ArrayList<Double>();
		ArrayList<Double> sqDiff = new ArrayList<Double>();
		Collection<Annotation> incorrectWindows = new TreeSet<Annotation>();
		double totalDiff = 0;
		double variance = 0;
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-r1", "first test result file name", true);
		p.addStringArg("-r2", "second test result file name", true);
		
		p.parse(args);
		
		String result1 = p.getStringArg("-r1");
		String result2 = p.getStringArg("-r2");
		
		BufferedReader br1 = new BufferedReader(new FileReader(result1));
		BufferedReader br2 = new BufferedReader(new FileReader(result2));
		
		//skip header
		br1.readLine();
		br2.readLine();
		
		while (br1.ready()&&br2.ready()) {
			String[] line1 = br1.readLine().split("\t");
			String[] line2 = br2.readLine().split("\t");
			
			double d = Double.parseDouble(line1[5])-Double.parseDouble(line2[5]);
			if (d!=0) {
				incorrectWindows.add(new BasicAnnotation(line1[2],Integer.parseInt(line1[3]),Integer.parseInt(line1[4])));
			}
			diff.add(Double.parseDouble(line1[5])-Double.parseDouble(line2[5]));
			totalDiff = totalDiff + Double.parseDouble(line1[5])-Double.parseDouble(line2[5]);
			sqDiff.add(Math.pow(Double.parseDouble(line1[5])-Double.parseDouble(line2[5]),2));
			variance = variance + Math.pow(Double.parseDouble(line1[5])-Double.parseDouble(line2[5]),2);
		}
		logger.info("sum of differences between the two files:\t" + totalDiff);
		logger.info("");
		logger.info("variance between the two files:\t" + variance);
		logger.info("");
		logger.info(incorrectWindows);
	}

}
