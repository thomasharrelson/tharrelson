package secondReadBias;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import broad.core.parser.CommandLineParser;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;
import broad.pda.seq.protection.WigSequenceComposition;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.readers.WigReader;
import nextgen.core.writers.WigWriter;

public class PileUpSequenceComposition {
	private Map<String,Collection<Gene>> chrRegions;
	private WigReader wr;
	private boolean zScoreMode;
	private WigSequenceComposition wsc;
	private Map<String, Sequence> genome;
	private double zScore;
	private Map<String,Map<Integer,Double>> pileUpPositions;
	
	public PileUpSequenceComposition(String bedFile, String wigFile, String genomeFile, boolean z, int shift, double zscore) throws IOException {
		chrRegions = BEDFileParser.loadDataByChr(new File(bedFile));
		wr = new WigReader(wigFile);
		genome = FastaSequenceIO.getChrSequencesFromFasta(genomeFile);
		this.zScoreMode = z;
		this.zScore = zscore;
		this.pileUpPositions = findPileUps();
		wsc = new WigSequenceComposition(pileUpPositions, genome, chrRegions, shift);
	}
	
	public void writePileUps(String outFile) throws IOException {
		wsc.writeBaseCountsToTable(outFile);
	}
	
	public void writePileUpsToWig(String outFile) throws IOException {
		FileWriter w = new FileWriter(outFile);
		Map<String,Map<Integer,Double>> wrData = wr.getAllValues();
		for (String chr : pileUpPositions.keySet()) {
			Map<Integer,Double> chrPileUps = pileUpPositions.get(chr);
			w.write("variableStep chrom=" + chr + "\n");
			for (int pos : chrPileUps.keySet()) {
				if (wrData.get(chr).containsKey(pos)) {
					w.write(WigWriter.coordinateToWigPosition(pos) + "\t" + wrData.get(chr).get(pos) + "\n");
				}
			}
		}
		w.close();
	}
	
	private Map<String,Map<Integer,Double>> findPileUps() {
		Map<String,Map<Integer,Double>> rtrn = new TreeMap<String,Map<Integer,Double>>();
		
		for (String chr : chrRegions.keySet()) {
			Map<Integer,Double> chrPileUps = new HashMap<Integer,Double>();
			for (Gene region : chrRegions.get(chr)) {
				if (zScoreMode) {
					Map<Integer,Double> zScoreMap = getZscores(region);
					if (zScoreMap == null) {continue;}
					for (int pos : zScoreMap.keySet()) {
						if (zScoreMap.get(pos) >= zScore) {
							
							chrPileUps.put(pos, 1.0);
						}
					}
				} else {
					Map<Integer,Double> scoreMap = wr.getValues(region);
					int greatestCountPos = 0;
					int greatestCount = 0;
					for (int pos : scoreMap.keySet()) {
						if (scoreMap.get(pos)>greatestCount) {
							greatestCountPos = pos;
						}
					}
					
					chrPileUps.put(greatestCountPos, 1.0);
				}
			}
			rtrn.put(chr, chrPileUps);
		}
		return rtrn;
	}
	
	private Map<Integer,Double> getZscores(Annotation region) {
		Map<Integer,Double> rtrn = new TreeMap<Integer,Double>();
		Map<Integer, Double> scoreMap = wr.getValues(region);
		if (scoreMap == null) {return null;}
		double avg = getCount(region)/((double) region.size());
		double stdDev = getStdDev(scoreMap, avg);
		for (int pos : scoreMap.keySet()) {
			rtrn.put(pos, (scoreMap.get(pos)-avg)/stdDev);
		}
		return rtrn;
	}
	
	private Double getStdDev(Map<Integer,Double> scoreMap, double avg) {
		double total = 0;
		for (double val : scoreMap.values()) {
			total += Math.pow((val-avg), 2);
		}
		return Math.sqrt(total/((double) scoreMap.values().size()));
	}
	
	private Double getCount(Annotation region) {
		Map<Integer,Double> scoreMap = wr.getValues(region);
		double rtrn = 0;
		for (Double d : scoreMap.values()) {
			rtrn += d.doubleValue();
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
		p.addStringArg("-b", "bed file of either peaks or genes", true);
		p.addStringArg("-w", "wig file for second/first read positions", true);
		p.addStringArg("-g", "genome fasta file", true);
		p.addStringArg("-o", "name of output file", true);
		p.addBooleanArg("-z", "true if you want to use z scores to find pile-ups; otherwise, the position with the highest count will be used", false, false);
		p.addIntArg("-s", "shift from pileup (default is 1; used for read2 bias analysis)", false,1);
		p.addDoubleArg("-zs", "standard score that will be used to define a pileup; any position with a z-score >= 6 will be a peak", false, 6.0);
		
		p.parse(args);
		
		String bedFile = p.getStringArg("-b");
		String wigFile = p.getStringArg("-w");
		String genomeFile = p.getStringArg("-g");
		String outFile = p.getStringArg("-o");
		int shift = p.getIntArg("-s");
		boolean z = p.getBooleanArg("-z");
		double zscore = p.getDoubleArg("-zs");
		
		PileUpSequenceComposition pusc = new PileUpSequenceComposition(bedFile, wigFile, genomeFile, z, shift, zscore);
		pusc.writePileUps(outFile);
	}

}
