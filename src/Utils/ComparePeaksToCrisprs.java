package Utils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.parser.CommandLineParser;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;

public class ComparePeaksToCrisprs {
	static Logger logger = Logger.getLogger(ComparePeaksToCrisprs.class.getName());
	
	private static Annotation getCutSite(Annotation crispr) {
		Annotation rtrn;
		if (crispr.isNegativeStrand()) {
			rtrn = new BasicAnnotation(crispr.getChr(), crispr.getStart()+3, crispr.getStart()+9, crispr.getOrientation());
		} else {
			rtrn = new BasicAnnotation(crispr.getChr(), crispr.getEnd()-9, crispr.getEnd()-3, crispr.getOrientation());
		}
		return rtrn;
	}
	
	private static String getCrisprSequence(Annotation crispr, Map<String, Sequence> genome) {
		Sequence chrGenome = genome.get(crispr.getChr());
		if (crispr.isNegativeStrand()) {
			return chrGenome.getRegion(crispr.getStart(), crispr.getEnd()).getAntisense().getSequenceBases();
		} else {
			return chrGenome.getRegion(crispr.getStart(), crispr.getEnd()).getSequenceBases();
		}
	}

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub

		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-p", "bed file of peaks",true);
		p.addStringArg("-c", "bed file of known crisprs", true);
		p.addStringArg("-f", "genome fasta file", true);
		p.addStringArg("-o", "output file", true);
		
		p.parse(args);
		
		String peakFile = p.getStringArg("-p");
		String crisprFile = p.getStringArg("-c");
		String genomeFile = p.getStringArg("-f");
		String output = p.getStringArg("-o");
		
		Map<String, Collection<Gene>> peaksByChr = BEDFileParser.loadDataByChr(new File(peakFile));
		Map<String, IntervalTree<Gene>> crisprsByChr = BEDFileParser.loadDataByChrToTree(new File(crisprFile));
		Map<String, Sequence> genome = FastaSequenceIO.getChrSequencesFromFasta(genomeFile);
		
		//Map<Annotation, Collection<Annotation>> peakToCrisprMap = new TreeMap<Annotation, Collection<Annotation>>();
		FileWriter w1 = new FileWriter(output + "_peakToCrisprs.txt");
		FileWriter w2 = new FileWriter(output + "_overlappedCrisprs.txt");
		FileWriter w3 = new FileWriter(output + "_overlappedPeaks.bed");
		
		
		for (String chr : peaksByChr.keySet()) {
			logger.info("comparing peaks to crisprs for " + chr);
			for (Gene peak : peaksByChr.get(chr)) {
				Iterator<Node<Gene>> crisprOverlappers = crisprsByChr.get(chr).overlappers(peak.getStart(), peak.getEnd());
				w1.write(peak.getName());
				while (crisprOverlappers.hasNext()) {
					Gene crispr = crisprOverlappers.next().getValue();
					Annotation crisprCut = getCutSite(crispr);
					if (peak.overlaps(crisprCut, false)) {
						w1.write("\t" + crispr.getName());
						String crisprBases = getCrisprSequence(crispr, genome);
						w2.write(peak.getName() + "\t" + crisprBases + "\n");
						w3.write(peak.toBED() + "\n");
					}
				}
				w1.write("\n");
			}
		}
		w1.close();
		w2.close();
		w3.close();
	}

}
