package secondReadBias;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.AbstractPairedEndAlignment;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;
import nextgen.core.readFilters.GenomicSpanFilter;

import broad.core.math.Statistics;
import broad.core.parser.CommandLineParser;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.pda.annotation.BEDFileParser;

public class AlignmentWidthBias {
	
	static Logger logger = Logger.getLogger(AlignmentWidthBias.class.getName());
	
	private AlignmentWidthBias() {}
	
	public static char reverseBase(char initBase) {
		switch (initBase) {
		case 'A': case 'a': return 'T';
		case 'G': case 'g': return 'C';
		case 'C': case 'c': return 'G';
		case 'T': case 't': return 'A';
		case 'N': return 'N';
		default: return 'z';
		}
	}
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "bam file", true);
		p.addStringArg("-g", "gene bed file", true);
		p.addStringArg("-f", "genome fasta file", true);
		p.addStringArg("-p", "bed file of peaks to iterate over", true);
		p.addStringArg("-o", "output file name, no extension", false, "alignmentWidthBiasOutput");
		
		p.parse(args);
		
		String bamFile = p.getStringArg("-b");
		String bedFile = p.getStringArg("-g");
		String genomeFile = p.getStringArg("-f");
		String peakFile = p.getStringArg("-p");
		String outFile = p.getStringArg("-o");
		
		Map<String,Sequence> genome = FastaSequenceIO.getChrSequencesFromFasta(genomeFile);
		Map<String,Collection<Gene>> genesByChr = BEDFileParser.loadDataByChr(new File(bedFile));
		Map<String,Collection<Gene>> peaksByChr = BEDFileParser.loadDataByChr(new File(peakFile));
		TranscriptomeSpace tspace = new TranscriptomeSpace(genesByChr);
		AlignmentModel data = new AlignmentModel(bamFile,tspace,true);
		data.addFilter(new GenomicSpanFilter(1000));
		
		AlignmentWidthBias awb = new AlignmentWidthBias();
		
		Map<Integer,Integer> fragmentAsizes = new TreeMap<Integer,Integer>();
		Map<Integer,Integer> fragmentCsizes = new TreeMap<Integer,Integer>();
		Map<Integer,Integer> fragmentGsizes = new TreeMap<Integer,Integer>();
		Map<Integer,Integer> fragmentTsizes = new TreeMap<Integer,Integer>();
		Map<String,Map<Double,Double>> peakAmedians = new TreeMap<String,Map<Double,Double>>();
		Map<String,Map<Double,Double>> peakCmedians = new TreeMap<String,Map<Double,Double>>();
		Map<String,Map<Double,Double>> peakGmedians = new TreeMap<String,Map<Double,Double>>();
		Map<String,Map<Double,Double>> peakTmedians = new TreeMap<String,Map<Double,Double>>();
		ArrayList<PeakStats> Astats = new ArrayList<PeakStats>();
		ArrayList<PeakStats> Cstats = new ArrayList<PeakStats>();
		ArrayList<PeakStats> Gstats = new ArrayList<PeakStats>();
		ArrayList<PeakStats> Tstats = new ArrayList<PeakStats>();
		
		double totalCSize = 0;
		double totalCReads = 0;
		double totalASize = 0;
		double totalAReads = 0;
		double totalGSize = 0;
		double totalGReads = 0;
		double totalTSize = 0;
		double totalTReads = 0;
		
		int pos;
		AbstractPairedEndAlignment read;
		CloseableIterator<Alignment> readItr;
		Sequence chrSeq;
		Collection<Double> Asizes = new TreeSet<Double>();
		Collection<Double> Csizes = new TreeSet<Double>();
		Collection<Double> Gsizes = new TreeSet<Double>();
		Collection<Double> Tsizes = new TreeSet<Double>();
		
		for (String chr : peaksByChr.keySet()) {
			chrSeq = genome.get(chr);
			String chrBases = chrSeq.getSequenceBases();
			logger.info("GETTING FRAGMENTS FOR " + chr);
			for (Gene g : peaksByChr.get(chr)) {
				logger.info("getting fragment sizes for gene : " + g);
				// get iterator over reads
				readItr = data.getOverlappingAnnotations(g, false);
				String geneBases;
				
				/*
				if (g.isNegativeStrand()) {
					geneBases = chrSeq.getRegion(g.getStart(),g.getEnd()).getAntisense().getSequenceBases();
				} else {
					geneBases = chrSeq.getRegion(g.getStart(), g.getEnd()).getSequenceBases();
				}
				*/
				
				while (readItr.hasNext()) {
					read = (AbstractPairedEndAlignment) readItr.next();
					Strand orientation;
					Alignment firstMate = read.getFirstMate();
					Alignment secondMate = read.getSecondMate();
					if (read.toSAMRecord().getSecondOfPairFlag()) {
						orientation = firstMate.getFragmentStrand();
					} else {
						orientation = secondMate.getFragmentStrand();
					}
					char transcribedBase;
					//logger.info("getting base");
					if (orientation == Strand.NEGATIVE) {
						pos = read.getEnd();
						if (pos > chrBases.length()) { continue;}
						
						char tmpBase = chrBases.charAt(pos);
						transcribedBase = reverseBase(tmpBase);
					} else { 
						pos = read.getStart()-1;
						if (pos < 0) { continue;}
						
						transcribedBase = chrBases.charAt(pos);
					}
					//logger.info("got base");
					if (transcribedBase=='A' || transcribedBase=='a') {
						Collection<Integer> tmpSizes = read.getFragmentSize(tspace);
						for (int tmp : tmpSizes) {
							Asizes.add((double) tmp);
						}
						totalASize += read.getSize();
						totalAReads ++;
						/*
						peakCTsizes += read.getSize();
						peakCTreads ++;
						*/
						if (fragmentAsizes.containsKey(read.getSize())) {
							int count = fragmentAsizes.get(read.getSize());
							fragmentAsizes.put(read.getSize(), count+1);
						} else {
							fragmentAsizes.put(read.getSize(), 1);
						}
					} else if (transcribedBase=='C' || transcribedBase=='c'){
						Collection<Integer> tmpSizes = read.getFragmentSize(tspace);
						for (int tmp : tmpSizes) {
							Csizes.add((double) tmp);
						}
						totalCSize += read.getSize();
						totalCReads ++;
						/*
						peakAGsizes += read.getSize();
						peakAGreads ++;
						*/
						if (fragmentCsizes.containsKey(read.getSize())) {
							int count = fragmentCsizes.get(read.getSize());
							fragmentCsizes.put(read.getSize(), count+1);
						} else {
							fragmentCsizes.put(read.getSize(), 1);
						}
					} else if (transcribedBase=='G' || transcribedBase=='g'){
						Collection<Integer> tmpSizes = read.getFragmentSize(tspace);
						for (int tmp : tmpSizes) {
							Gsizes.add((double) tmp);
						}
						totalGSize += read.getSize();
						totalGReads ++;
						/*
						peakAGsizes += read.getSize();
						peakAGreads ++;
						*/
						if (fragmentGsizes.containsKey(read.getSize())) {
							int count = fragmentGsizes.get(read.getSize());
							fragmentGsizes.put(read.getSize(), count+1);
						} else {
							fragmentGsizes.put(read.getSize(), 1);
						}
					} else if (transcribedBase=='T' || transcribedBase=='t'){
						Collection<Integer> tmpSizes = read.getFragmentSize(tspace);
						for (int tmp : tmpSizes) {
							Tsizes.add((double) tmp);
						}
						totalTSize += read.getSize();
						totalTReads ++;
						/*
						peakAGsizes += read.getSize();
						peakAGreads ++;
						*/
						if (fragmentTsizes.containsKey(read.getSize())) {
							int count = fragmentTsizes.get(read.getSize());
							fragmentTsizes.put(read.getSize(), count+1);
						} else {
							fragmentTsizes.put(read.getSize(), 1);
						}
					}
				}
				if (Asizes.size()!=0) {/*
					Map<Double,Double> tmpMap = new TreeMap<Double,Double>();
					tmpMap.put(Statistics.median(Asizes),(double) Asizes.size());
					peakAmedians.put(g.getName(), tmpMap);
					*/
					Astats.add(awb.new PeakStats(g, Statistics.median(Asizes), (double) Asizes.size()));
				} else { /*
					Map<Double,Double> tmpMap = new TreeMap<Double,Double>();
					tmpMap.put(0.0, 0.0);
					peakAmedians.put(g.getName(), tmpMap);
					*/
					Astats.add(awb.new PeakStats(g, 0.0, 0.0));
				}
				
				if (Csizes.size()!=0) {/*
					Map<Double,Double> tmpMap = new TreeMap<Double,Double>();
					tmpMap.put(Statistics.median(Csizes),(double) Csizes.size());
					peakCmedians.put(g.getName(), tmpMap);
					*/
					Cstats.add(awb.new PeakStats(g, Statistics.median(Csizes), (double) Csizes.size()));
				} else { /*
					Map<Double,Double> tmpMap = new TreeMap<Double,Double>();
					tmpMap.put(0.0, 0.0);
					peakCmedians.put(g.getName(), tmpMap);
					*/
					Cstats.add(awb.new PeakStats(g, 0.0, 0.0));
				}
				
				if (Gsizes.size()!=0) {
					/*
					Map<Double,Double> tmpMap = new TreeMap<Double,Double>();
					tmpMap.put(Statistics.median(Gsizes),(double) Gsizes.size());
					peakGmedians.put(g.getName(), tmpMap);
					*/
					Gstats.add(awb.new PeakStats(g, Statistics.median(Gsizes), (double) Gsizes.size()));
				} else { /*
					Map<Double,Double> tmpMap = new TreeMap<Double,Double>();
					tmpMap.put(0.0, 0.0);
					peakGmedians.put(g.getName(), tmpMap);
					*/
					Gstats.add(awb.new PeakStats(g, 0.0, 0.0));
				}
				
				if (Tsizes.size()!=0) {/*
					Map<Double,Double> tmpMap = new TreeMap<Double,Double>();
					tmpMap.put(Statistics.median(Tsizes),(double) Tsizes.size());
					peakTmedians.put(g.getName(), tmpMap);
					*/
					Tstats.add(awb.new PeakStats(g, Statistics.median(Tsizes), (double) Tsizes.size()));
				} else {/* 
					Map<Double,Double> tmpMap = new TreeMap<Double,Double>();
					tmpMap.put(0.0, 0.0);
					peakTmedians.put(g.getName(), tmpMap);
					*/
					Tstats.add(awb.new PeakStats(g, 0.0, 0.0));
				}
				
				Asizes.clear();
				Csizes.clear();
				Gsizes.clear();
				Tsizes.clear();
			}
		}
		/*
		logger.info("Average fragment size ending in C/T = " + (totalCTSize/totalCTReads));
		logger.info("Average fragment size ending in A/G = " + (totalAGSize/totalAGReads));
		*/
		
		FileWriter wDistA = new FileWriter(outFile + "_fragSizeDist_A.txt");
		FileWriter wDistC = new FileWriter(outFile + "_fragSizeDist_C.txt");
		FileWriter wDistG = new FileWriter(outFile + "_fragSizeDist_G.txt");
		FileWriter wDistT = new FileWriter(outFile + "_fragSizeDist_T.txt");
		FileWriter wMedA = new FileWriter(outFile + "_peakSizeMedians_A.txt");
		FileWriter wMedC = new FileWriter(outFile + "_peakSizeMedians_C.txt");
		FileWriter wMedG = new FileWriter(outFile + "_peakSizeMedians_G.txt");
		FileWriter wMedT = new FileWriter(outFile + "_peakSizeMedians_T.txt");
		
		for (Integer size : fragmentAsizes.keySet()) {
			wDistA.write(size + "\t" + fragmentAsizes.get(size) + "\n");
		}
		for (Integer size : fragmentCsizes.keySet()) {
			wDistC.write(size + "\t" + fragmentCsizes.get(size) + "\n");
		}
		for (Integer size : fragmentGsizes.keySet()) {
			wDistG.write(size + "\t" + fragmentGsizes.get(size) + "\n");
		}
		for (Integer size : fragmentTsizes.keySet()) {
			wDistT.write(size + "\t" + fragmentTsizes.get(size) + "\n");
		}
		for (PeakStats stats : Astats) {
			wMedA.write(stats.getPeak().getName() + "\t" + stats.getMedian() + "\t" + stats.getTotal() + "\n");
		}
		for (PeakStats stats : Cstats) {
			wMedC.write(stats.getPeak().getName() + "\t" + stats.getMedian() + "\t" + stats.getTotal() + "\n");
		}
		for (PeakStats stats : Gstats) {
			wMedG.write(stats.getPeak().getName() + "\t" + stats.getMedian() + "\t" + stats.getTotal() + "\n");
		}
		for (PeakStats stats : Tstats) {
			wMedT.write(stats.getPeak().getName() + "\t" + stats.getMedian() + "\t" + stats.getTotal() + "\n");
		}
		
		wDistA.close();
		wDistC.close();
		wDistG.close();
		wDistT.close();
		
		wMedA.close();
		wMedC.close();
		wMedG.close();
		wMedT.close();
	}
	
	public class PeakStats {
		double median;
		double total;
		Gene peak;
		public PeakStats(Gene g, double median, double total) {
			this.peak = g;
			this.median = median;
			this.total = total;
		}
		
		private double getMedian() {
			return this.median;
		}
		
		private double getTotal() {
			return this.total;
		}
		
		private Gene getPeak() {
			return peak;
		}
		
	}

}
