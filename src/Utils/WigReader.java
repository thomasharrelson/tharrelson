package Utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Annotation;


public class WigReader {
	private BufferedReader br;
	private String file;
	private TreeMap<String, TreeMap<Integer,Double>> wigData;
	private static Logger logger = Logger.getLogger(WigReader.class.getName());
	
	//Constructor - should construct solely from a wig file name
	public WigReader(String wigFile) throws IOException {
		this.file = wigFile;
		br = new BufferedReader(new FileReader(wigFile));
		this.wigData = new TreeMap<String,TreeMap<Integer,Double>>();
		// initialize loop by finding first chromosome
		String chr = null;
		while (br.ready()) {
			logger.info("Entered initialization of wig file: this message should appear only once");
			String line = br.readLine();
			if (line.contains("variableStep chrom=")) {
				chr = checkChr(line);
				logger.info("got first chromosome: " + chr);
				break;
			}
		}
		while (br.ready()) {	
			TreeMap<String,TreeMap<Integer,Double>> chrCounts = readChrCounts(chr);
			for (String newChr : chrCounts.keySet()) {
				this.wigData.put(chr, chrCounts.get(newChr));
				chr = newChr;
			}
		}
		br.close();
		return;
	}
	
	private TreeMap<String,TreeMap<Integer,Double>> readChrCounts(String chr) throws NumberFormatException, IOException {
		TreeMap<Integer,Double> rtrn = new TreeMap<Integer,Double>();
		TreeMap<String,TreeMap<Integer,Double>> fullRtrn = new TreeMap<String,TreeMap<Integer,Double>>();
		String line = null;
		if (br.ready()) {
			line = br.readLine();
		}
		while (br.ready() && !line.contains("variableStep chrom=")) {
			
			String[] splitLine = line.split("\t");
			if (splitLine.length!=2) { continue;}
			rtrn.put(Integer.parseInt(splitLine[0]), Double.parseDouble(splitLine[1]));	
			line = br.readLine();
		}
		String newChr = "end";
		if (line.contains("variableStep chrom=")) {
			newChr = checkChr(line);
		}
		fullRtrn.put(newChr, rtrn);
		logger.info("About to return newly created wig data map");
		logger.info("keyset size: " + rtrn.keySet().size());
		logger.info("values size: " + rtrn.values().size());
		logger.info("newChr = " + newChr + "\n");
		return fullRtrn;
	}
	
	public Map<Integer,Double> getWigData(String chr) {
		return wigData.get(chr);
	}
	
	public Map<Integer,Double> getWigData(String chr, int start, int end) {
		return wigData.get(chr).subMap(start, end);
	}
	
	public Map<Integer,Double> getWigData(String chr, int pos) {
		return getWigData(chr, pos, (pos+1));
	}
	
	public Map<Integer,Double> getWigData(Annotation region) throws IOException {
		Iterator<? extends Annotation> exons = region.getBlocks().iterator();
		Map<Integer,Double> rtrn = new TreeMap<Integer,Double>();
		while (exons.hasNext()) {
			Annotation exon = exons.next();
			rtrn.putAll(getWigData(exon.getChr(),exon.getStart(),exon.getEnd()));
		}
		return rtrn;
	}
	
	private String checkChr(String line) {
		return line.substring(19);
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// The purpose of this main method is to test the wig reader w/ a variety of examples.
		WigReader wr1 = new WigReader("/seq/lincRNA/Tommy/WigWriterRuns/mergedSamples_1.5K_norm.wig");
		logger.info("reading from full chromosome (No. 1)");
		logger.info("chromosome 1 keyset size (should be number of positions w/ counts>0): " + wr1.getWigData("chr1").keySet().size());
		logger.info("reading from section of chromosome 5 (5505191 - 5505653)\n" + wr1.getWigData("chr5",5505191,5505653));
		logger.info("reading single point from chromosome 1 (3661061): " + wr1.getWigData("chr1",3661061));
	}

}
