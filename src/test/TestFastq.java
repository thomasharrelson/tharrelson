package test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.logging.Logger;



import broad.core.parser.CommandLineParser;

public class TestFastq {
	public static Logger logger = Logger.getLogger(TestFastq.class.getName());
	private static String reverseComplement(String bar) {
		String rtrn = "";
		for (int i=bar.length()-1;i>=0;i--) {
			switch (bar.charAt(i)) {
				case 'A': rtrn = rtrn + "T"; break;
				case 'T': rtrn = rtrn + "A"; break;
				case 'G': rtrn = rtrn + "C"; break;
				case 'C': rtrn = rtrn + "G"; break;
			}
		}
		return rtrn;
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-f","fastq file to check",true);
		p.addStringArg("-b","barcode to check for, default is ...", false, "put barcode here");
		p.parse(args);
		String file = p.getStringArg("-f");
		String barcode = p.getStringArg("-b");
		BufferedReader br = new BufferedReader(new FileReader(file));
		String revBarcode = reverseComplement(barcode);
		int forCount = 0;
		int revCount = 0;
		while (br.ready()) {
			String line = br.readLine();
			if (line.contains(barcode)) {
				forCount = forCount+1;
			}
			if (line.contains(revBarcode)) {
				revCount = revCount + 1;
			}
		}
		logger.info("number of times barcode appeared in fastq files: " + forCount);
		logger.info("number of times reverse complement of barcode appeared in fastq file: " + revCount);
		br.close();
	}

}
