package secondReadBias;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.parser.CommandLineParser;
import broad.pda.annotation.BEDFileParser;

public class PeakDifferences {
	
	Map<String, Collection<Peak>> peakSets;
	String outFile;

	public PeakDifferences(String Abed, String Cbed, String Gbed, String Tbed, String totalBed, String outFile) throws IOException {
		Map<String, IntervalTree<Gene>> peaksByChrA = BEDFileParser.loadDataByChrToTree(new File(Abed));
		Map<String, IntervalTree<Gene>> peaksByChrC = BEDFileParser.loadDataByChrToTree(new File(Cbed));
		Map<String, IntervalTree<Gene>> peaksByChrG = BEDFileParser.loadDataByChrToTree(new File(Gbed));
		Map<String, IntervalTree<Gene>> peaksByChrT = BEDFileParser.loadDataByChrToTree(new File(Tbed));
		Map<String, Collection<Gene>> peaksByChrTotal = BEDFileParser.loadDataByChr(new File(totalBed));
		peakSets = matchPeaks(peaksByChrA, peaksByChrC, peaksByChrG, peaksByChrT, peaksByChrTotal);
		this.outFile = outFile;
	}
	
	private Map<String, Collection<Peak>> matchPeaks(Map<String, IntervalTree<Gene>> peaksByChrA, Map<String, IntervalTree<Gene>> peaksByChrC, Map<String, IntervalTree<Gene>> peaksByChrG, Map<String, IntervalTree<Gene>> peaksByChrT, Map<String, Collection<Gene>> peaksByChrTotal) {
		Map<String, Collection<Peak>> rtrn = new TreeMap<String, Collection<Peak>>();
		for (String chr : peaksByChrTotal.keySet()) {
			IntervalTree<Gene> chrTreeA = peaksByChrA.get(chr);
			IntervalTree<Gene> chrTreeC = peaksByChrC.get(chr);
			IntervalTree<Gene> chrTreeG = peaksByChrG.get(chr);
			IntervalTree<Gene> chrTreeT = peaksByChrT.get(chr);
			Collection<Peak> chrPeaks = new TreeSet<Peak>();
			for (Gene peak : peaksByChrTotal.get(chr)) {
				
				Gene peakA = matchTrack(peak, chrTreeA);
				Gene peakC = matchTrack(peak, chrTreeC);
				Gene peakG = matchTrack(peak, chrTreeG);
				Gene peakT = matchTrack(peak, chrTreeT);
				Map<String, Annotation> splitPeaks = new TreeMap<String, Annotation>();
				if (peakA!=null) {
					splitPeaks.put("Apeak", peakA);
				}
				if (peakC!=null) {
					splitPeaks.put("Cpeak", peakC);
				}
				if (peakG!=null) {
					splitPeaks.put("Gpeak", peakG);
				}
				if (peakT!=null) {
					splitPeaks.put("Tpeak", peakT);
				}
				
				Peak splitUpPeak = new Peak(peak, splitPeaks);
				chrPeaks.add(splitUpPeak);
			}
			rtrn.put(chr, chrPeaks);
		}
		return rtrn;
	}
	
	private Gene matchTrack(Gene peak, IntervalTree<Gene> chrTree) {
		Gene largePeak = new Gene(peak);
		largePeak.setStart(peak.getStart()-50);
		largePeak.setEnd(peak.getEnd()+50);
		Iterator<Node<Gene>> overlappers = chrTree.overlappers(largePeak.getStart(), largePeak.getEnd());
		Gene rtrn = null;
		while (overlappers.hasNext()) {
			Gene checkGene = overlappers.next().getValue();
			if (peak.overlaps(checkGene)) {
				if (rtrn==null) {
					rtrn = checkGene;
				} else {
					if (peak.intersect(checkGene).size() > peak.intersect(rtrn).size()) {
						rtrn = checkGene;
					}
				}
			}
		}
		if (rtrn == null) {
			while (overlappers.hasNext()) {
				Gene checkGene = overlappers.next().getValue();
				if (largePeak.overlaps(checkGene)) {
					if (rtrn==null) {
						rtrn = checkGene;
					} else {
						if (largePeak.intersect(checkGene).size() > largePeak.intersect(rtrn).size()) {
							rtrn = checkGene;
						}
					}
				}
			}
		}
		return rtrn;
	}
	
	public void writePeakSizeDifferences() throws IOException {
		FileWriter w = new FileWriter(outFile + "_sizeDiff.txt");
		w.write("peak name\tC-A size diff\tC-G size diff\tC-T size diff\n");
		for (String chr : peakSets.keySet()) {
			for (Peak p : peakSets.get(chr)) {
				w.write(p.getName() + "\t" + p.getCAsizeDiff() + "\t" + p.getCGsizeDiff() + "\t" + p.getCTsizeDiff() + "\n");
			}
		}
		w.close();
	}
	
	public void writePeakShifts() throws IOException {
		FileWriter w = new FileWriter(outFile + "_peakShift.txt");
		w.write("peak name\tC-A shift\tC-G shift\tC-T shift\n");
		for (String chr : peakSets.keySet()) {
			for (Peak p : peakSets.get(chr)) {
				w.write(p.getName() + "\t" + p.getCAshift() + "\t" + p.getCGshift() + "\t" + p.getCTshift() + "\n");
			}
		}
		w.close();
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub

		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-o", "output file prefix (no .txt)", true);
		p.addStringArg("-a", "bed file of 'A' peaks", true);
		p.addStringArg("-c", "bed file of 'C' peaks", true);
		p.addStringArg("-g", "bed file of 'G' peaks", true);
		p.addStringArg("-t", "bed file of 'T' peaks", true);
		p.addStringArg("-b", "bed file of peaks not split by sequence", true);
		
		p.parse(args);
		
		String outFile = p.getStringArg("-o");
		String Abed = p.getStringArg("-a");
		String Cbed = p.getStringArg("-c");
		String Gbed = p.getStringArg("-g");
		String Tbed = p.getStringArg("-t");
		String totalBed = p.getStringArg("-b");
		
		PeakDifferences pd = new PeakDifferences(Abed, Cbed, Gbed, Tbed, totalBed, outFile);
		pd.writePeakShifts();
		pd.writePeakSizeDifferences();
	}
	

}
