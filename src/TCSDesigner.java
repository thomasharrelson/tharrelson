import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.feature.GenomeWindow;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;
import broad.core.primer3.Primer3Configuration;
import broad.core.primer3.Primer3ConfigurationFactory;
import broad.core.primer3.Primer3IO;
import broad.core.primer3.Primer3SequenceInputTags;
import broad.core.primer3.Primer3SequenceInputTags.SequenceRegionCoordinates;
import broad.core.primer3.PrimerPair;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.core.util.PipelineUtils;
import broad.pda.annotation.BEDFileParser;


public class TCSDesigner {
	
	/*
	 * Parameters:
	 * genes2cut - collection of gene objects to find cut sites for.
	 * filename - name of the file containing the gene names to be cut.
	 * 
	 */
	private Map<String,Collection<Gene>> genes2cut;
	private static String filename;
	private static Map<String,Sequence> genome;
	private static String genomeFastaFile = "/seq/lincRNA/data/mm9/mm9.fa";
	private static String h3k4me3File = "/seq/lincRNA/Tommy/Data/H3K4me3_peaks_300_150_60.bed";
	private static String h3k4me1File = "/seq/lincRNA/Tommy/Data/H3K4me1_60_150_300_peaks.bed";
	//private static String homologyFile = "/seq/lincRNA/Tommy/Data/mm9_pi_lods.sorted.bed";
	private static String homologyFile = "/seq/lincRNA/Tommy/Data/mm9_pi_lods.sorted.bed";
	private static String repeatFile = "/seq/lincRNA/Tommy/Data/H3K4me1_60_150_300_peaks.bed";
	private static String geneBedFile = "/seq/lincRNA/Pam/Databases/mouse/annotations/RefSeq_LincV3_ChromatinWithNames_nonrandom.bed";

	private Map<Gene,GenomeWindow> broadRegions;
	private Map<Gene,ArrayList<GenomeWindow>> cutRegions;	// ask Pam about which set would be best here.
	private Map<String,IntervalTree<Gene>> chrGenes;
	private Map<String,IntervalTree<Gene>> h3k4me3Track;// = new TreeMap<String,Collection<Gene>>();		// track that finds promoters
	private Map<String,IntervalTree<Gene>> h3k4me1Track;// = new TreeMap<String,Collection<Gene>>();		// track that finds enhancers
	private Map<String,IntervalTree<Gene>> homologyTrack;// = new TreeMap<String,Collection<Gene>>();		// can help find enhancers and promoters but is misleading; however do not cut in highly conserved regions
	private Map<String,IntervalTree<Gene>> repeatTrack;// = new TreeMap<String,Collection<Gene>>();		// don't use cut regions that have a lot of repeats.
	private Collection<TreeMap<String,Collection<Gene>>> trackBundle = new TreeSet<TreeMap<String,Collection<Gene>>>();
	
	private static Map<Gene,ArrayList<String>> flaggedGenes = new TreeMap<Gene,ArrayList<String>>(); // the set of genes that did not yield proper CRISPRs, and associated string describing problem
	
	Logger logger = Logger.getLogger(TCSDesigner.class.getName());
	
	public TCSDesigner(String file, char fileType) throws IOException {
		/* 
		 * fileType will designate the type of input file and thus the type of constructor to be used.
		 */
		filename = file;
		genome = FastaSequenceIO.getChrSequencesFromFasta(genomeFastaFile);
		BufferedReader r1 = new BufferedReader(new FileReader(filename));
		
		logger.info("setting tracks...");
		h3k4me3Track=(BEDFileParser.loadDataByChrToTree(new File(h3k4me3File)));
		h3k4me1Track=(BEDFileParser.loadDataByChrToTree(new File(h3k4me1File)));
		homologyTrack=(BEDFileParser.loadDataByChrToTree(new File(homologyFile)));
		repeatTrack=(BEDFileParser.loadDataByChrToTree(new File(repeatFile)));
		chrGenes = BEDFileParser.loadDataByChrToTree(new File(geneBedFile));
		/*
		trackBundle.add(h3k4me3Track);
		trackBundle.add(h3k4me1Track);
		trackBundle.add(homologyTrack);
		trackBundle.add(repeatTrack);
		*/
		logger.info("Finished setting tracks.");
		
		// construct genes2cut from file
		switch (fileType) {
		case 'n': nameConstruct(r1);break;
		case 'p': posConstruct(r1);break;
		case 'b': bedConstruct(filename);break;
		default: r1.close(); throw new IllegalArgumentException(fileType + " is not a valid file type");
		}
		// use genes2cut to create broadRegions
	}
	
	public void writeCRISPRsWithSurvAssay(String outFile, int numCRISPRs, int numPrimers) throws IOException, InterruptedException {
		
		CRISPRwriter w = new CRISPRwriter(outFile, numCRISPRs, numPrimers);
		//Runtime.getRuntime().exec("use Blat");
		for (String chr : genes2cut.keySet()) {
			
			IntervalTree<Gene> chrMe3Genes = h3k4me3Track.get(chr);//StaticMethods.makeTree(h3k4me3Track.get(chr));
			IntervalTree<Gene> chrMe1Genes = h3k4me1Track.get(chr);//StaticMethods.makeTree(h3k4me1Track.get(chr));
			IntervalTree<Gene> chrHomolGenes = homologyTrack.get(chr);//StaticMethods.makeTree(homologyTrack.get(chr));
			IntervalTree<Gene> chrRepeatGenes = repeatTrack.get(chr);//StaticMethods.makeTree(repeatTrack.get(chr));
			/*
			chrMe3Genes.addAll(h3k4me3Track.get(chr));
			chrMe1Genes.addAll(h3k4me1Track.get(chr));
			chrHomolGenes.addAll(homologyTrack.get(chr));
			chrRepeatGenes.addAll(repeatTrack.get(chr));
			*/
			
			Map<String,IntervalTree<Gene>> chrBundle = new HashMap<String,IntervalTree<Gene>>();
			chrBundle.put("RPT",chrRepeatGenes);
			chrBundle.put("CONS",chrHomolGenes);
			chrBundle.put("ME1",chrMe1Genes);
			chrBundle.put("ME3",chrMe3Genes);
			//chrBundle.put("EXONS",chrGenes.get(chr));
			
			Collection<Gene> genes = genes2cut.get(chr);
			logger.info("Getting CRISPRs for " + chr);
			Map<Gene,Map<String,Map<Integer,SequenceRegion>>> chrCRISPRs = findAllCRISPRs(genes, chrBundle, chrMe3Genes, chrMe1Genes, chrRepeatGenes);
			logger.info("Found all CRISPRs for " + chr);
			
			Sequence chrSeq = genome.get(chr);
			logger.info("Getting primers for " + chr);
			Map<Gene,Map<String,Collection<PrimerPair>>> chrPrimers = findAllPrimers(chrCRISPRs, chrSeq);
			logger.info("Found all primers for " + chr);
			
			logger.info("Writing CRISPRs and primers to files...");
			//Map<Gene,Map<String,Map<Integer,SequenceRegion>>> crisprs = findAllCRISPRs(genes, chrBundle, chrMe3Genes, chrMe1Genes, chrRepeatGenes);
			for (Gene g : chrCRISPRs.keySet()) {
				w.writeGene(g, chrCRISPRs.get(g), numCRISPRs);
			}
			for (Gene g : chrPrimers.keySet()) {
				w.writePrimers(g, chrPrimers.get(g), numPrimers);
			}
			
		}
		w.close();
	}
	
	private Map<Gene,Map<String,Collection<PrimerPair>>> findAllPrimers(Map<Gene,Map<String,Map<Integer,SequenceRegion>>> chrCRISPRs, Sequence chrSeq) throws IOException {
		Map<Gene,Map<String,Collection<PrimerPair>>> rtrn = new TreeMap<Gene,Map<String,Collection<PrimerPair>>>();
		for (Gene g : chrCRISPRs.keySet()) {
			Map<String,Collection<PrimerPair>> genePrimers = new TreeMap<String,Collection<PrimerPair>>();
			for (String id : chrCRISPRs.get(g).keySet()) {
				if (!id.equals("gblock")) {
					Collection<SequenceRegion> cuts = chrCRISPRs.get(g).get(id).values();
					if (!cuts.isEmpty()) {
						Map<String,Collection<PrimerPair>> primers = findPrimers(cuts,chrSeq,id);
						genePrimers.putAll(primers);
					} else {
						//Map<String,Collection<PrimerPair>> dummyPrimers = new TreeMap<String,Collection<PrimerPair>>();
						genePrimers.put(id, new TreeSet<PrimerPair>());
					}
				}
			}
			rtrn.put(g, genePrimers);
		}
		return rtrn;
	}
	
	public Map<Gene,Map<String,Map<Integer,SequenceRegion>>> findAllCRISPRs(Collection<Gene> genes, Map<String,IntervalTree<Gene>> chrBundle, IntervalTree<Gene> chrMe3Regions, IntervalTree<Gene> chrMe1Regions, IntervalTree<Gene> chrRepeatRegions) throws IOException, InterruptedException {
		Map<Gene,Map<String,Map<Integer,SequenceRegion>>> rtrn = new TreeMap<Gene,Map<String,Map<Integer,SequenceRegion>>>();
		
		
		for (Gene g : genes) {
			logger.info("Finding CRISPRs for locus " + g.getName());
			Map<String,Map<Integer,SequenceRegion>> CRISPRs = findCRISPRsites(g, chrBundle, chrMe3Regions, chrMe1Regions, chrRepeatRegions, 4);
			rtrn.put(g,CRISPRs);
		}
		return rtrn;
	}
	
	public Map<String,Map<Integer,SequenceRegion>> findCRISPRsites(Gene locus, Map<String,IntervalTree<Gene>> chrBundle, IntervalTree<Gene> chrMe3Regions, IntervalTree<Gene> chrMe1Regions, IntervalTree<Gene> chrRepeatRegions, int numCRISPRs) throws IOException, InterruptedException {
		Map<String,Map<Integer,SequenceRegion>> rtrn = new TreeMap<String,Map<Integer,SequenceRegion>>();
		// find intron site
		List<? extends Annotation> introns = getFirstIntron(locus);
		Annotation firstIntron = introns.get(0);
		Annotation TSSwindow;
		int PEcutStart = getPEposition(locus, chrMe1Regions, chrMe3Regions);
		boolean error = checkTSSforError(locus,chrMe3Regions);
		IntervalTree<Gene> tmpMe3 = null;
		if (error) {
			TSSwindow = getTSSwindow(locus, chrBundle.get("CONS"));
			tmpMe3 = chrBundle.remove("ME3");
		} else {
			TSSwindow = getTSSwindow(locus, chrBundle.get("ME3"));
		}
		double checkPos=0;
		if (locus.isNegativeStrand()) {
			
			if (TSSwindow.getStart()-firstIntron.getEnd()>470) {
				checkPos = TSSwindow.getStart();
			} else {
				checkPos = (470-(TSSwindow.getStart()-firstIntron.getEnd()))/2 + TSSwindow.getStart();
			}
		} else {
			if (firstIntron.getStart()-TSSwindow.getEnd()>470) {
				checkPos = TSSwindow.getEnd(); 
			} else {
				checkPos = TSSwindow.getEnd() - (470-(firstIntron.getStart()-TSSwindow.getEnd()))/2 - 1;
			}
		}
		Map<Integer,SequenceRegion> INcuts = findINsequences(introns, locus, chrBundle, chrRepeatRegions, 4, checkPos);
		rtrn.put("IN", INcuts);
		logger.info("intron cuts: " + INcuts);
		TreeSet<SequenceRegion> in = new TreeSet<SequenceRegion>();
		in.addAll(INcuts.values());
		
		if (locus.isNegativeStrand()) {
			if (!in.isEmpty()) {
				checkPos = getCutPosition(in.first());
			} else {
				checkPos = firstIntron.getEnd()-25;
			}
		} else {
			if (!in.isEmpty()) {
				checkPos = getCutPosition(in.last());
			} else {
				checkPos = firstIntron.getStart() + 25;
			}
		}
		// check for problems after getting each set of cuts; i.e. not enough CRISPRs returned
		// check if me3 data overlapped the tss for this loci; if it did re-create the tss window using conservation track
		// and remove me3 track from bundle because it will hinder finding small gblocks bt intron and tss cut sites
		// remember to put the me3 track back before designing the next loci
		
		Map<Integer, SequenceRegion> TSScuts = findTSSsequences(TSSwindow, locus, chrBundle, chrRepeatRegions, 4, checkPos);
		rtrn.put("TSS", TSScuts);
		logger.info("TSS cuts: " + TSScuts);
		Map<Integer, SequenceRegion> PEcuts = findPEsequences(PEcutStart,100,locus,chrBundle,chrRepeatRegions,4);
		rtrn.put("PE", PEcuts);
		logger.info("PE cuts " + PEcuts);
		if (INcuts.size()<numCRISPRs) {
			int offset = numCRISPRs - INcuts.size();
			if (flaggedGenes.containsKey(locus)) {
				ArrayList<String> tmp = flaggedGenes.get(locus);
				tmp.add("Missing " + offset + " CRISPRs at the intron cut site.");
				flaggedGenes.put(locus,tmp);
			} else {
				ArrayList<String> tmp = new ArrayList<String>();
				tmp.add("Missing " + offset + " CRISPRs at the intron cut site.");
				flaggedGenes.put(locus,tmp);
			}
		}
		if (TSScuts.size()<numCRISPRs) {
			int offset = numCRISPRs - TSScuts.size();
			if (flaggedGenes.containsKey(locus)) {
				ArrayList<String> tmp = flaggedGenes.get(locus);
				tmp.add("Missing " + offset + " CRISPRs at the window between the promoter and the TSS.");
				flaggedGenes.put(locus,tmp);
			} else {
				ArrayList<String> tmp = new ArrayList<String>();
				tmp.add("Missing " + offset + " CRISPRs at the window between the promoter and the TSS.");
				flaggedGenes.put(locus,tmp);
			}
		}
		if (PEcuts.size()<numCRISPRs) {
			int offset = numCRISPRs - PEcuts.size();
			if (flaggedGenes.containsKey(locus)) {
				ArrayList<String> tmp = flaggedGenes.get(locus);
				tmp.add("Missing " + offset + " CRISPRs at the window upstream from the proximal enhancer.");
				flaggedGenes.put(locus,tmp);
			} else {
				ArrayList<String> tmp = new ArrayList<String>();
				tmp.add("Missing " + offset + " CRISPRs at the window upstream from the proximal enhancer.");
				flaggedGenes.put(locus,tmp);
			}
		}
		double size = Math.min(TSScuts.size(), INcuts.size());
		Map<Integer,SequenceRegion> gBlocks = getGblocks(TSScuts,INcuts, locus.getChr());
		rtrn.put("gblock", gBlocks);
		for (int i=1; i<=size; i++) {
			double TSScut = getCutPosition(TSScuts.get(i));
			double INcut = getCutPosition(INcuts.get(i));
			if (Math.abs(TSScut-INcut)>470) {
				if (flaggedGenes.containsKey(locus)) {
					ArrayList<String> tmp = flaggedGenes.get(locus);
					tmp.add("TSS cut site and IN cut site are more than 500 bp apart for CRISPR #" + i);
					flaggedGenes.put(locus,tmp);
				} else {
					ArrayList<String> tmp = new ArrayList<String>();
					tmp.add("TSS cut site and IN cut site are more than 500 bp apart for CRISPR #" + i);
					flaggedGenes.put(locus,tmp);
				}
			}
		}
		if (error && tmpMe3!=null) {
			chrBundle.put("ME3", tmpMe3);
		}
		return rtrn;
	}
	
	private Map<Integer,SequenceRegion> getGblocks(Map<Integer,SequenceRegion> TSScuts, Map<Integer,SequenceRegion> INcuts, String chr) {
		Map<Integer, SequenceRegion> rtrn = new TreeMap<Integer, SequenceRegion>();
		int numBlocks = Math.min(TSScuts.size(),INcuts.size());
		for (int i=0;i<numBlocks;i++) {
			Strand orientation = INcuts.get(i+1).getOrientation();
			int INcutPos = (int) Math.floor(getCutPosition(INcuts.get(i+1)));
			int TSScutPos = (int) Math.floor(getCutPosition(TSScuts.get(i+1)));
			int start = Math.min(INcutPos,TSScutPos);
			int end = Math.max(INcutPos,TSScutPos);
			if (orientation == Strand.POSITIVE) {
				SequenceRegion tmp = genome.get(chr).getRegion(start,end);
				if (hasBamHI(tmp)) {
					tmp.setSequenceBases("AACCTGCAGG" + tmp.getSequenceBases() + "AGATCTAA");
				} else {
					tmp.setSequenceBases("AACCTGCAGG" + tmp.getSequenceBases() + "GGATCCAA");
				}
				rtrn.put(i+1, tmp);
			} else {
				SequenceRegion tmp = genome.get(chr).getRegion(start,end);
				tmp.reverse();
				if (hasBamHI(tmp)) {
					tmp.setSequenceBases("AACCTGCAGG" + tmp.getSequenceBases() + "AGATCTAA");
				} else {
					tmp.setSequenceBases("AACCTGCAGG" + tmp.getSequenceBases() + "GGATCCAA");
				}
				rtrn.put(i+1, tmp);
			}
			
		}
		return rtrn;
	}
	
	private boolean hasBamHI(Sequence tmp) {
		String bases = tmp.getSequenceBases();
		for (int i=0; i<=bases.length()-6;i++) {
			if (bases.substring(i,i+6).equalsIgnoreCase("GGATCC")) {
				return true;
			}
		}
		return false;
	}
	
	private boolean checkTSSforError(Gene locus, IntervalTree<Gene> me3Regions) {
		Annotation TSSstart;
		if (locus.isNegativeStrand()) {
			TSSstart = new BasicAnnotation(locus.getChr(), locus.getEnd(), locus.getEnd()+1);
		} else {
			TSSstart = new BasicAnnotation(locus.getChr(), locus.getStart()-1, locus.getStart());
		}
		Iterator<Node<Gene>> me3Itr = me3Regions.overlappers(TSSstart.getStart(), TSSstart.getEnd());
		while (me3Itr.hasNext()) {
			Gene me3Region = me3Itr.next().getValue();
			if (me3Region.overlaps(TSSstart)) {
				if (flaggedGenes.containsKey(locus)) {
					ArrayList<String> tmp = flaggedGenes.get(locus);
					tmp.add("The txn start overlaps a promoter region; program will remove me3 track, and use conservation track to find TSS window");
					flaggedGenes.put(locus,tmp);
				} else {
					ArrayList<String> tmp = new ArrayList<String>();
					tmp.add("The txn start overlaps a promoter region; program will remove me3 track, and use conservation track to find TSS window");
					flaggedGenes.put(locus,tmp);
				}
				return true;
			}
		}
		return false;
	}
	
	private Iterator<Annotation> getWindowIterator(Annotation window, int winSize) {
		TreeSet<Annotation> windows = new TreeSet<Annotation>();

		int pos = window.getStart();
		while (pos < window.getEnd()) {
			if (window.getEnd() - pos < winSize) {
				GenomeWindow gw = new GenomeWindow(window.getChr(), pos, window.getEnd());
				gw.setOrientation(window.getOrientation());
				windows.add(gw);
				pos = window.getEnd();
			} else {
				GenomeWindow gw = new GenomeWindow(window.getChr(), pos, pos + winSize);
				gw.setOrientation(window.getOrientation());
				windows.add(gw);
				pos = pos + winSize;
			}
		}
		//logger.info("windows: " + windows);
		if (window.isNegativeStrand()) {
			return windows.descendingIterator();
		}
		return windows.iterator();
	}
	
	private Map<Integer,SequenceRegion> findPEsequences(int PEwindowStart, int winSize, Gene locus, Map<String,IntervalTree<Gene>> chrBundle, IntervalTree<Gene> repeatTrack, int numCRISPRs) throws IOException, InterruptedException {
		Map<Integer, SequenceRegion> rtrn = new TreeMap<Integer, SequenceRegion>();
		// remember to create TSSwindow to be antisense to the actual gene of interest because I have to search backwards from the start site
		Iterator<Annotation> winItr = getWindowIterator(locus.getChr(), PEwindowStart, 100, locus.getStrand().getReverseStrand());
		logger.info("finding PE sequences for " + locus + "using PE start: " + PEwindowStart);
		TreeSet<SequenceRegion> allCutsConsidered = new TreeSet<SequenceRegion>();
		while (winItr.hasNext() && allCutsConsidered.size()<(3*numCRISPRs)) {
			Annotation window = winItr.next();
			logger.info("initial window: " + window);
			Annotation filteredWindow = filterWindow(window, chrBundle);
			logger.info("filtered window: " + filteredWindow);
			if (filteredWindow == null) continue;
			TreeSet<SequenceRegion> potentialCuts = findAllCutSites(filteredWindow, genome.get(window.getChr()));
			logger.info("potential cuts: " + potentialCuts);
			if (potentialCuts.isEmpty()) continue;
			TreeSet<SequenceRegion> filteredCuts = filterSet(potentialCuts, repeatTrack);
			logger.info("filtered cuts: " + filteredCuts);
			if (filteredCuts.isEmpty()) continue;
			TreeSet<SequenceRegion> blatFilteredCuts = filterBLAT(filteredCuts, locus.getName()+"_PE");
			logger.info("blat filtered cuts: " + blatFilteredCuts);
			if (blatFilteredCuts.isEmpty()) continue;
			allCutsConsidered.addAll(blatFilteredCuts);
		}
		if (!allCutsConsidered.isEmpty()) {
				// find ideal cut location
			ArrayList<SequenceRegion> rankedSites = rankCutSites(allCutsConsidered, locus.getStrand().getReverseStrand(),PEwindowStart, numCRISPRs,0);
			logger.info("ranked cut sites: " + rankedSites);
			for (int i = 1; i<=rankedSites.size(); i++) {
				SequenceRegion region = rankedSites.get(i-1);
				String name = region.getName();
				rtrn.put(Integer.parseInt(name.substring(name.length()-1, name.length())), region);
				if (rtrn.size()==numCRISPRs) {
					return rtrn;
				}
			}
		}
		return rtrn;
	}
	
	private Iterator<Annotation> getWindowIterator(String chr, int start, int winSize, Strand direction) {
		// direction should be the order you want to search the windows in.
		TreeSet<Annotation> windows = new TreeSet<Annotation>();
		Annotation tmp;
		// set total search space to be 2000; possibly include this in the parameters
		double numWindows = Math.ceil(2000/winSize);
		if (direction.toString().equals("-")) {
			for (int i = 0; i < numWindows; i++) {
				tmp = new BasicAnnotation(chr, start-winSize+1, start+1);
				start = start-winSize;
				windows.add(tmp);
			}
			return windows.descendingIterator();
		} else {
			for (int i = 0; i < numWindows; i++) {
				tmp = new BasicAnnotation(chr, start, start+winSize);
				start = start+winSize;
				windows.add(tmp);
			}
		}
		return windows.iterator();
	}
	
	private Map<Integer,SequenceRegion> findTSSsequences(Annotation TSSwindow, Gene locus, Map<String,IntervalTree<Gene>> chrBundle, IntervalTree<Gene> repeatTrack, int numCRISPRs, double checkPos) throws IOException, InterruptedException {
		Map<Integer, SequenceRegion> rtrn = new TreeMap<Integer, SequenceRegion>();
		// remember to create TSSwindow to be antisense to the actual gene of interest because I have to search backwards from the start site
		Iterator<Annotation> winItr = getWindowIterator(TSSwindow, 100);
		logger.info("finding TSS sequences for " + locus + "using TSS window: " + TSSwindow);
		int idealCutLocation;
		if (TSSwindow.isNegativeStrand()) {
			idealCutLocation = TSSwindow.getEnd() - 3;
		} else {
			idealCutLocation = TSSwindow.getStart() + 2;
		}
		TreeSet<SequenceRegion> allCutsConsidered = new TreeSet<SequenceRegion>();
		while (winItr.hasNext() && allCutsConsidered.size()<(3*numCRISPRs)) {
			Annotation window = winItr.next();
			logger.info("initial window: " + window);
			Annotation filteredWindow = filterWindow(window, chrBundle);
			logger.info("filtered window: " + filteredWindow);
			if (filteredWindow==null) continue;
			TreeSet<SequenceRegion> potentialCuts = findAllCutSites(filteredWindow, genome.get(window.getChr()));
			logger.info("potential cuts: " + potentialCuts);
			if (potentialCuts.isEmpty()) continue;
			TreeSet<SequenceRegion> filteredCuts = filterSet(potentialCuts, repeatTrack);
			logger.info("filtered cuts: " + filteredCuts);
			if (filteredCuts.isEmpty()) continue;
			TreeSet<SequenceRegion> blatFilteredCuts = filterBLAT(filteredCuts, locus.getName()+"_TSS");
			logger.info("blat filtered cuts: " + blatFilteredCuts);
			if (blatFilteredCuts.isEmpty()) continue;
			allCutsConsidered.addAll(blatFilteredCuts);
		}
		if (!allCutsConsidered.isEmpty()) {
			ArrayList<SequenceRegion> rankedSites = rankCutSites(allCutsConsidered, TSSwindow.getStrand(), idealCutLocation, numCRISPRs, checkPos);
			logger.info("ranked cut sites: " + rankedSites);
			for (int i = 1; i<=rankedSites.size(); i++) {
				SequenceRegion region = rankedSites.get(i-1);
				String name = region.getName();
				rtrn.put(Integer.parseInt(name.substring(name.length()-1, name.length())), region);
				if (rtrn.size()==numCRISPRs) {
					return rtrn;
				}
			}
		}
		return rtrn;
	}
	
	private Map<Integer,SequenceRegion> findINsequences(List<? extends Annotation> introns, Gene locus, Map<String,IntervalTree<Gene>> chrBundle, IntervalTree<Gene> repeatTrack, int numCRISPRs, double checkPos) throws IOException, InterruptedException {
		Annotation window;
		Annotation firstIntron = introns.get(0);
		logger.info("intron: " + firstIntron);
		TreeMap<Integer,SequenceRegion> rtrn = new TreeMap<Integer,SequenceRegion>();
		// need to re-design; set the window to exclude the first and last 25 bp of intron;
		// if intron is very big, cut the window to a reasonable size
		// if intron is less than 50 bp, flag the corresponding locus, and don't continue.
		// if intron has fewer than 4 cuts, check the size to make sure whole thing was checked.
		// if it wasn't, check the remainder of the 
		Iterator<Annotation> winItr;
		logger.info("finding intron sequences for " + locus + "using intron: " + firstIntron);
		
		int idealCutLocation;
		if (locus.isNegativeStrand() && firstIntron.size()>200) {
			idealCutLocation = firstIntron.getEnd()-103;
		} else if (firstIntron.size()>200) {
			idealCutLocation = firstIntron.getStart() + 102;
		} else if (locus.isNegativeStrand() && firstIntron.size()<=200) {
			idealCutLocation = (int) Math.floor(firstIntron.size()/2);
			idealCutLocation = firstIntron.getReferenceCoordinateAtPosition(idealCutLocation);
		} else {
			idealCutLocation = (int) Math.ceil(firstIntron.size()/2);
			idealCutLocation = firstIntron.getReferenceCoordinateAtPosition(idealCutLocation);
		}
		
		//Iterator<? extends Annotation> intronItr = introns.iterator();
		for (Annotation intron : introns) {
			//intron = intronItr.next();
			if (intron.size()>50) {
				Annotation win = new BasicAnnotation(intron.getChr(), intron.getStart() + 25, intron.getEnd()-25,intron.getOrientation());
				winItr = getWindowIterator(win,150);
			} else {
				//TODO: create flag structure, and return null structure
				logger.warn("Intron smaller than 50 bp, returning no cut sites");
				continue;
			}
			logger.info("ideal cut location: " + idealCutLocation);
			TreeSet<SequenceRegion> allCutsConsidered = new TreeSet<SequenceRegion>();
			while (winItr.hasNext() && allCutsConsidered.size()<(3*numCRISPRs)) {
				window = winItr.next();
				logger.info("initial window: " + window);
				Annotation filteredWindow = filterWindow(window, chrBundle);
				logger.info("filteredWindow: " + filteredWindow);
				if (filteredWindow == null) { continue;}
				TreeSet<SequenceRegion> potentialCuts = findAllCutSites(filteredWindow, genome.get(window.getChr()));
				logger.info("potential cuts : " + potentialCuts);
				if (potentialCuts.isEmpty()) continue;
				TreeSet<SequenceRegion> filteredCuts = filterSet(potentialCuts, repeatTrack);
				logger.info("filteredCuts: " + filteredCuts);
				if (filteredCuts.isEmpty()) continue;
				TreeSet<SequenceRegion> blatFilteredCuts = filterBLAT(filteredCuts, locus.getName()+"_IN");
				logger.info("blat filtered cuts: " + blatFilteredCuts);
				if (blatFilteredCuts.isEmpty()) continue;
				allCutsConsidered.addAll(blatFilteredCuts);
			}
			if (!allCutsConsidered.isEmpty()) {
				ArrayList<SequenceRegion> rankedSites = rankCutSites(allCutsConsidered, intron.getStrand(), idealCutLocation, numCRISPRs, checkPos);
				logger.info("ranked cut sites: " + rankedSites);
				for (int i = 1; i<=rankedSites.size(); i++) {
					SequenceRegion region = rankedSites.get(i-1);
					String name = region.getName();
					rtrn.put(Integer.parseInt(name.substring(name.length()-1, name.length())), region);
					if (rtrn.size()==numCRISPRs) {
						return rtrn;
					}
				}
			}
		}
		return rtrn;
	}
	
	private ArrayList<SequenceRegion> rankCutSites(TreeSet<SequenceRegion> blatFilteredCuts, Strand orientation, int refPos, int numCuts, double checkPos) {
		ArrayList<SequenceRegion> rtrn = new ArrayList<SequenceRegion>();
		NavigableMap<Double,Set<SequenceRegion>> scoreMap = new TreeMap<Double,Set<SequenceRegion>>();
		for (SequenceRegion cut : blatFilteredCuts) {
			double score = findScore(cut, orientation, refPos, checkPos);
			if (scoreMap.containsKey(score)) {
				Set<SequenceRegion> tmp = scoreMap.get(score);
				tmp.add(cut);
				scoreMap.put(score, tmp);
			} else {
				TreeSet<SequenceRegion> initList = new TreeSet<SequenceRegion>();
				initList.add(cut);
				scoreMap.put(score, initList);
			}
		}
		// order scores from greatest to least
		scoreMap = scoreMap.descendingMap();
		//initialize counter for number of cut sites returned
		int i = 0;
		for (double key : scoreMap.keySet()) {
			Set<SequenceRegion> tmp = scoreMap.get(key);
			for (SequenceRegion region : tmp) {
				// re-set name of the region to reflect the proper rank of the cut site
				
				String name = region.getName();
				region.setName(name.substring(0, name.length()-2) + "_" + Integer.toString(i+1));
				
				rtrn.add(region);
				i++;
				if (i==numCuts) {break;}
			}
		}
		return rtrn;
	}
	
	private double findScore(SequenceRegion cut, Strand strand, int refPos, double checkPos) {
		//double checkPos;
		double cutPos;
		double score = cut.getScore();
		if (checkPos!=0) {
			if (strand.equals(Strand.NEGATIVE)) {
				cutPos = getCutPosition(cut);
				if (Math.abs(checkPos-cutPos)>470) {
					score = score - 0.5; // apply score reduction that will decrease even a perfect score below the previously minimum allowed score
				}
			} else {
				//checkPos = window.getStart() + relPos + 2.5;
				cutPos = getCutPosition(cut);
				if (Math.abs(checkPos-cutPos)>470) {
					score = score - 0.5; // apply score reduction that will decrease even a perfect score below the previously minimum allowed score
				}
			}
		}
		cutPos = getCutPosition(cut);
		double diff = Math.abs(cutPos-refPos);
		score = score-0.0025*diff;
		
		logger.info("cut score: " + cut.getName() + "\t" + score);
		return score;
	}
	
	private double getCutPosition(SequenceRegion cut) {
		if (cut.isNegativeStrand()) {
			return cut.getStart() + 2.5;
		} else {
			return cut.getEnd() - 3.5;
		}
	}
	
	private TreeSet<SequenceRegion> filterBLAT(TreeSet<SequenceRegion> filteredRegions, String name) throws IOException, InterruptedException {
		String jobID = Long.valueOf(System.currentTimeMillis()).toString();
		Map<String,SequenceRegion> nameMap = writeSequencesToFile(filteredRegions, name, jobID);
		
		ArrayList<String> jobIDs = new ArrayList<String>();
		
		
		String cmmd = "blat " + genomeFastaFile + " " + name + "_" + jobID + "_tmp.fa " + name + "_" + jobID + "_tmp.psl " + "-noHead -ooc=/seq/lincRNA/Tommy/TCSdesign/mm9.11.ooc -minIdentity=80 -minMatch=1 -maxIntron=0 -minScore=15";
		jobIDs.add(jobID);
		PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, name + jobID + ".bsub", "hour", 4);
		
		logger.info("Waiting for blat jobs to finish...");
		PipelineUtils.waitForAllJobs(jobIDs, Runtime.getRuntime());
		logger.info("All jobs finished");
		File blatInput = new File(name + "_" + jobID + "_tmp.fa");
		if (!blatInput.delete()) {
			logger.warn(name + "_" + jobID + "_tmp.fa" + " could not be deleted");
		}
		return readBLAT(name+"_"+jobID+"_tmp.psl", nameMap);
	}
	
	private TreeSet<SequenceRegion> readBLAT(String BLATfile, Map<String, SequenceRegion> nameMap) throws IOException {
		TreeSet<SequenceRegion> rtrn = new TreeSet<SequenceRegion>();
		Map<String,ArrayList<Integer>> blatOutput = new TreeMap<String,ArrayList<Integer>>();
		Map<String, ArrayList<Annotation>> blatAnnotations = new TreeMap<String,ArrayList<Annotation>>();
		BufferedReader br = new BufferedReader(new FileReader(BLATfile));
		
		while (br.ready()) {
			String line = br.readLine();
			//logger.info("blat output line: " + line);
			// TODO: what is the delimiter for psl files?
			// going with //s+ for now; taken from internet forum
			String[] fields = line.split("\t");
			int matches = Integer.parseInt(fields[0]);
			String name = fields[9];
			if (blatOutput.containsKey(name)) {
				ArrayList<Integer> tmp = blatOutput.get(name);
				tmp.add(matches);
				blatOutput.put(name, tmp);
				ArrayList<Annotation> tmpAnnos = blatAnnotations.get(name);
				Annotation anno = new BasicAnnotation(fields[13],Integer.parseInt(fields[15]),Integer.parseInt(fields[16]));
				anno.setOrientation(fields[8].charAt(0));
				tmpAnnos.add(anno);
				blatAnnotations.put(name, tmpAnnos);
			} else {
				ArrayList<Integer> matchList = new ArrayList<Integer>();
				matchList.add(matches);
				blatOutput.put(name, matchList);
				ArrayList<Annotation> annoList = new ArrayList<Annotation>();
				Annotation anno = new BasicAnnotation(fields[13],Integer.parseInt(fields[15]),Integer.parseInt(fields[16]));
				anno.setOrientation(fields[8].charAt(0));
				annoList.add(anno);
				blatAnnotations.put(name, annoList);
			}
		}
		br.close();
		File blatFile = new File(BLATfile);
		if (!blatFile.delete()) {
			logger.warn(BLATfile + " could not be deleted");
		}
		
		for (String cutSite : blatOutput.keySet()) {
			ArrayList<Integer> matchList = blatOutput.get(cutSite);
			ArrayList<Annotation> targets = blatAnnotations.get(cutSite);
			//String name = cutSite.substring(0, cutSite.length()-4);
			//String id = cutSite.substring(cutSite.length()-5,cutSite.length()-4);
			if (matchList.size()==1) {
				rtrn.add(nameMap.get(cutSite));
			} else if (matchList.size()>1) {
				int size = nameMap.get(cutSite).size();
				int i = 0; // number of perfect matches to target 20mer
				for (Integer match : matchList) {
					if ((20-match.intValue()) == 0) {
						i++;
					}
					if (i>1) {break;}
				}
				/*
				 * TODO: a single mismatch is not enough to cause cutting in a different location
				 * If the mismatch resides in the last 12 bases, then the CRISPR will not cut.
				 * If it resides in the first 8, it will cut at around 50-90% efficiency.
				 * TODO: Figure out the minimum number of mismatches required to reduce cutting efficiency to ~0%
				 */
				if (i==1) {
					// only one perfect 20meric match
					// now check 
					boolean checkFlag = true;
					SequenceRegion region = nameMap.get(cutSite);
					ArrayList<Annotation> annoList = blatAnnotations.get(cutSite);
					double totalScore = 1;
					for (int j = 0; j<annoList.size(); j++) {
						Annotation anno = annoList.get(j);
						SequenceRegion test;
						if (anno.isNegativeStrand()) {
							test = genome.get(anno.getChr()).getRegion(anno.getStart(), anno.getStart()+20);
							test.setOrientation('-');
						} else {
							test = genome.get(anno.getChr()).getRegion(anno.getEnd()-20, anno.getEnd());
							test.setOrientation('+');
						}
						int matches = matchList.get(j);
						if (matches==20) { continue;}
						double score = testCRISPRvalidity(region, test, matches);
						if (score == 0.5) {checkFlag = false; break;}
						totalScore = totalScore*score;
					}
					if (checkFlag) {
						SequenceRegion sr = nameMap.get(cutSite);
						sr.setScore(totalScore);
						rtrn.add(sr);
					}
				}
			}
		}
		return rtrn;
	}
	
	private double testCRISPRvalidity(SequenceRegion original, SequenceRegion test, int matches) {
		double rtrn = 1;
		Sequence testSeq = test.getSequence();
		if (test.isNegativeStrand()) {
			testSeq = testSeq.getAntisense();
		}
		Sequence origSeq = original.getSequence();
		if (original.isNegativeStrand()) {
			origSeq = origSeq.getAntisense();
		}
		String origBases = origSeq.getSequenceBases();
		String testBases = testSeq.getSequenceBases();
		if (origBases.substring(8,20).equalsIgnoreCase(testBases.substring(8, 20)) && origBases.length()==testBases.length()) {
			return (1-Math.pow(0.5, origBases.length()-matches));
		}
		return rtrn;
	}
	
	private Map<String, SequenceRegion> writeSequencesToFile(Collection<SequenceRegion> filteredRegions, String outFile, String jobID) throws IOException {
		Map<String,SequenceRegion> nameMap = new TreeMap<String, SequenceRegion>();
		File file = new File(outFile + "_" + jobID + "_tmp.fa");
		FileWriter w = new FileWriter(file);
		int i = 1;
		for (SequenceRegion region : filteredRegions) {
			w.write(">" + outFile + "_" + i + "_tmp" + "\n");
			region.setName(outFile + "_" + Integer.toString(i));
			nameMap.put((outFile + "_" + Integer.toString(i) + "_tmp"), region);
			w.write(region.getSequenceBases() + "\n");
			i++;
		}
		w.close();
		return nameMap;
	}
	
	private TreeSet<SequenceRegion> filterSet(TreeSet<SequenceRegion> unfiltered, IntervalTree<Gene> repeatTree) {
		TreeSet<SequenceRegion> rtrn = new TreeSet<SequenceRegion>();
		//IntervalTree<Gene> repeatTree = StaticMethods.makeTree(repeatTrack);
		for (SequenceRegion unfilteredCutSite : unfiltered) {
			Iterator<Node<Gene>> rptOverlapper = repeatTree.overlappers(unfilteredCutSite.getStart(), unfilteredCutSite.getEnd());
			
			boolean rOverlapFlag = false;
			
			while (rptOverlapper.hasNext()) {
				if (rptOverlapper.next().getValue().overlapsStranded(unfilteredCutSite)) {
					rOverlapFlag = true;
					break;
				}
			}
			// Check that there are no homology and repeat overlaps of current window and add to return struc
			if (!rOverlapFlag) {
				rtrn.add(unfilteredCutSite);
			}
		}
		return rtrn;
	}
	
	private TreeSet<SequenceRegion> findAllCutSites(Annotation window, Sequence chrSeq) {
		TreeSet<SequenceRegion> rtrn = new TreeSet<SequenceRegion>();
		List<? extends Annotation> winBlocks = window.getBlocks(true);
		for (Annotation block : winBlocks) {
			SequenceRegion blockSeq = chrSeq.getRegion(block.getStart(), block.getEnd());
			for (int i=1; i<block.size()-2; i++) {
				SequenceRegion checkSeq1 = blockSeq.getRegion(i, i+2);
				SequenceRegion checkSeq2 = blockSeq.getRegion(i-1, i+1);
				if (checkSeq1.getSequenceBases().equalsIgnoreCase("GG") && i>7) {
					SequenceRegion tmp = chrSeq.getRegion(block.getStart()+i-21, block.getStart()+i-1);
					tmp.setOrientation(Strand.POSITIVE);
					rtrn.add(tmp);
				}
				if (checkSeq2.getSequenceBases().equalsIgnoreCase("CC") && i+7<block.size()) {
					SequenceRegion tmp = chrSeq.getRegion(block.getStart()+i+2,block.getStart()+i+22);
					tmp.reverse();
					tmp.setOrientation(Strand.NEGATIVE);
					rtrn.add(tmp);
				}
			}
		}
		return rtrn;
	}
	
	private Annotation filterWindow(Annotation window, TreeSet<Gene> homolTrack) {
		return filterWindow(window, homolTrack, null, null, null);
	}
	private Annotation filterWindow(Annotation window, TreeSet<Gene> homolTrack, TreeSet<Gene> repeatTrack) {
		return filterWindow(window, homolTrack, repeatTrack, null, null);
	}
	
	private Annotation filterWindow(Annotation window, TreeSet<Gene> homolTrack, TreeSet<Gene> repeatTrack, TreeSet<Gene> me3Track, TreeSet<Gene> me1Track) {
		Gene tmpWindow = new Gene(window);
		Iterator<Gene> homolOverlappers;
		Iterator<Gene> me3Overlappers;
		Iterator<Gene> me1Overlappers;
		Iterator<Gene> repeatOverlappers;
		if (homolTrack != null) {
			homolOverlappers = homolTrack.subSet(homolTrack.lower(tmpWindow), false, homolTrack.higher(tmpWindow), false).iterator();
			while (homolOverlappers.hasNext()) {
				Annotation homolPeak = homolOverlappers.next();
				window = window.minus(homolPeak);
			}
		}
		if (repeatTrack != null) {
			repeatOverlappers = repeatTrack.subSet(repeatTrack.lower(tmpWindow), false, repeatTrack.higher(tmpWindow), false).iterator();
			while (repeatOverlappers.hasNext()) {
				Annotation repeatPeak = repeatOverlappers.next();
				window = window.minus(repeatPeak);
			}
		}
		if (me3Track != null) {
			me3Overlappers = me3Track.subSet(me3Track.lower(tmpWindow), false, me3Track.higher(tmpWindow), false).iterator();
			while (me3Overlappers.hasNext()) {
				Annotation me3Peak = me3Overlappers.next();
				window = window.minus(me3Peak);
			}
		}
		if (me1Track != null) {
			me1Overlappers = me1Track.subSet(me1Track.lower(tmpWindow), false, me1Track.higher(tmpWindow), false).iterator();
			while (me1Overlappers.hasNext()) {
				Annotation me1Peak = me1Overlappers.next();
				window = window.minus(me1Peak);
				
			}
		}
		//Iterator<Gene> repeatOverlappers = repeatTrack.subSet(repeatTrack.lower(tmpWindow), false, repeatTrack.higher(tmpWindow), false).iterator();
		
		return window;
	}
	
	private Annotation filterWindow(Annotation window, Map<String,IntervalTree<Gene>> trackBundle) {
		//Gene tmpWindow = new Gene(window);
		String chr = window.getChr();
		Iterator<Node<Gene>> trackOverlappers;
		for (String trackID : trackBundle.keySet()) {
			IntervalTree<Gene> track = trackBundle.get(trackID);
			//trackOverlappers = track.subSet(track.lower(tmpWindow), false, track.higher(tmpWindow), false).iterator();
			trackOverlappers = track.overlappers(window.getStart(), window.getEnd());
			while (trackOverlappers.hasNext()) {
				Annotation homolPeak = trackOverlappers.next().getValue();
				window = window.minus(homolPeak);
				if (window.getBlocks().size()==0) {
					return null;
				}
			}
		
		}
		//Iterator<Gene> repeatOverlappers = repeatTrack.subSet(repeatTrack.lower(tmpWindow), false, repeatTrack.higher(tmpWindow), false).iterator();
		
		return window;
	}
	
	private int getPEposition(Gene locus, IntervalTree<Gene> chrMe1Regions, IntervalTree<Gene> chrMe3Regions) {
		if (locus.isNegativeStrand()) {
			//Gene TSS = new Gene(locus.getChr(),locus.getEnd(),locus.getEnd()+1);
			Gene promoter = chrMe3Regions.min(locus.getEnd(),locus.getEnd()+1).getValue();
			return chrMe1Regions.min(promoter.getStart(),promoter.getEnd()).getValue().getEnd();
		} else {
			//Gene TSS = new Gene(locus.getChr(),locus.getStart()-1,locus.getStart());
			Gene promoter = chrMe3Regions.max(locus.getStart()-1,locus.getStart()).getValue();
			return chrMe1Regions.max(promoter.getStart(),promoter.getEnd()).getValue().getStart();
		}
	}
	
	private Annotation getTSSwindow(Gene locus, IntervalTree<Gene> chrMe3Regions) {
		// find promoter location from h3k4me3 track; keep orientation in mind
		Annotation rtrn;
		if (locus.isNegativeStrand()) {
			//Gene TSS = new Gene(locus.getChr(),locus.getEnd(),locus.getEnd()+1);
			logger.info("Txn start site: " + locus.getName() + "\nEnd of promoter: " + chrMe3Regions.min(locus.getEnd(),locus.getEnd()+1).getStart());
			rtrn = new BasicAnnotation(locus.getChr(), locus.getEnd(), chrMe3Regions.min(locus.getEnd(),locus.getEnd()+1).getValue().getStart());
		} else {
			//only have to do this because Node.compare() is shitty for IntervalTree.max() but works well for IntervalTree.min()
			logger.info("Txn start site: " + locus.getName() + "\nEnd of promoter: " + chrMe3Regions.max(locus.getStart()-1,locus.getStart()).getStart());
			Annotation TSScheck = new BasicAnnotation(locus.getChr(), locus.getStart()-1, locus.getStart());
			Gene promoter = null;
			Iterator<Node<Gene>> promItr = chrMe3Regions.reverseIterator(locus.getStart()-1, locus.getStart());
			while (promItr.hasNext()) {
				promoter = promItr.next().getValue();
				if (!promoter.overlaps(TSScheck)) {
					break;
				} else {
					promoter = null;
				}
			}
			if (promoter != null) {
				rtrn = new BasicAnnotation(locus.getChr(), promoter.getEnd(), locus.getStart());
			} else {
				rtrn = new BasicAnnotation(locus.getChr(), 0, locus.getStart());
			}
		}
		rtrn.setOrientation(locus.getStrand().getReverseStrand());
		return rtrn;
	}
	
	private int getFirstIntronStartSite(Gene locus) {
		Annotation intron = getFirstIntron(locus).get(0);
		if (locus.isNegativeStrand()) {
			return intron.getEnd()-100; // do not want cut sites in the beginning part of the intron in case there are important splice features there.
		} else {
			return intron.getStart()+100; // ditto above
		}
	}
	
	private List<? extends Annotation> getFirstIntron(Gene gene) {
		Gene intronList = gene.getIntrons();
		if (intronList!=null)
			return intronList.getBlocks(true);
		else {
			List<Annotation> rtrn = new ArrayList<Annotation>();
			if (gene.isNegativeStrand()) {
				rtrn.add(new BasicAnnotation(gene.getChr(), gene.getStart()-1000, gene.getStart(),Strand.NEGATIVE));
				return rtrn;
			} else {
				rtrn.add(new BasicAnnotation(gene.getChr(), gene.getEnd(), gene.getEnd()+1000,Strand.POSITIVE));
				return rtrn;
			}
		}
	}
	
	private void nameConstruct(BufferedReader r1) {
		throw new IllegalArgumentException("Class not implemented for construction with a list of gene names");
	}
	
	//consider making public
	private void setBroadRegions(Collection<Gene> g2c) {
		Iterator<Gene> iter = g2c.iterator();
		TreeSet<GenomeWindow> cutRegions = new TreeSet<GenomeWindow>();
		while (iter.hasNext()) {
			Gene g = iter.next();
			int start;
			int stop;
			GenomeWindow gWin;
			if (g.isNegativeStrand()) {
				start = g.getStart();
				stop = g.getEnd() + 5000;
				
				Strand orientation = g.getStrand();
				// TODO: figure out how to use the enum Strand, so that I can get the proper window
				gWin = new GenomeWindow(new BasicAnnotation(g.getChr(),start,stop,orientation));
				gWin.addSourceAnnotation(gWin);
			} else {
				if (g.getStart()>=5000) {
					start = g.getStart() - 5000;
				} else {
					start = 0;
				}
				stop = g.getEnd();
				Strand orientation = g.getStrand();
				gWin = new GenomeWindow(new BasicAnnotation(g.getChr(),start,stop,orientation));
				gWin.addSourceAnnotation(g);
			}
			broadRegions.put(g, gWin);
		}
	}
	
	/*
	 * Preliminary builder of gene collection
	 * Uses a file that has format "chrX:(start int)-(stop int)"
	 * i.e. chr1:23456-28000
	 */
	private void posConstruct(BufferedReader r1) throws IOException {
		
		String line = r1.readLine();
		while (!line.isEmpty()) {
			String[] splitLine1 = line.split(":");
			String[] splitLine2 = splitLine1[1].split("-");
			Gene gene = new Gene(splitLine1[0],Integer.parseInt(splitLine2[0]),Integer.parseInt(splitLine2[1]));
			if (genes2cut.containsKey(gene.getChr())) {
				Collection<Gene> tmp = genes2cut.get(gene.getChr());
				tmp.add(gene);
				genes2cut.put(gene.getChr(), tmp);
			} else {
				Collection<Gene> chrGenes = new TreeSet<Gene>();
				chrGenes.add(gene);
				genes2cut.put(gene.getChr(), chrGenes);
			}
			line = r1.readLine();
		}
		return;
	}
	
	private void bedConstruct(String file) throws IOException {
		logger.info("Constructing designer from bed file of genes.");
		genes2cut = BEDFileParser.loadDataByChr(new File(file));
		return;
		//throw new IllegalArgumentException("Class not implemented for construction with a bed file of genes");
	}
	
	private Map<String,Collection<PrimerPair>> findPrimers(Collection<SequenceRegion> cuts, Sequence chrSeq, String cName) throws IOException {
		TreeMap<String,Collection<PrimerPair>> rtrn = new TreeMap<String, Collection<PrimerPair>>();
		
		// create proper configuration object
		Primer3Configuration config = Primer3ConfigurationFactory.getLongRangePCRConfiguration();
		config.minMeltingTemp = 60;
		config.maxProductSize = 3000;
		config.minProductSize = 600;
		config.optimalProductSize = 800;
		config.maxNumPrimersToReturn = 15;
		
		// create input for primer3; use whole chromosome seq, just exclude most of it
		Primer3SequenceInputTags input = new Primer3SequenceInputTags(chrSeq);
		// find most 5' cut
		int fivePrimePos = -1;
		int threePrimePos = -1;
		for (SequenceRegion cut : cuts) {
			
			int pos1;
			int pos2;
			if (cut.isNegativeStrand()) {
				pos1 = cut.getStart();
				pos2 = pos1 + 6;
			} else {
				pos1 = cut.getEnd();
				pos2 = pos1 - 6;
			}
			if (fivePrimePos == -1) { fivePrimePos = pos1;}
			if (threePrimePos == -1) { threePrimePos = pos2;}
			if (pos1 < fivePrimePos) { fivePrimePos = pos1;}
			if (pos2 > threePrimePos) { threePrimePos = pos2;}
		}
		logger.info("id = " + cName);
		logger.info("cut sites: " + cuts);
		logger.info("5prime position: " + fivePrimePos);
		logger.info("3prime position: " + threePrimePos);
		List<Primer3SequenceInputTags.SequenceRegionCoordinates> exCoord = new ArrayList<Primer3SequenceInputTags.SequenceRegionCoordinates>();
		exCoord.add(new SequenceRegionCoordinates(fivePrimePos - 100, threePrimePos + 100));
		exCoord.add(new SequenceRegionCoordinates(0, fivePrimePos - 690));
		exCoord.add(new SequenceRegionCoordinates(threePrimePos + 690, chrSeq.getLength()));
		input.addExcludedRegions(exCoord);
		Primer3IO p3io = new Primer3IO();
		p3io.startPrimer3Communications();
		Collection<PrimerPair> pp = p3io.findPrimerPair(input, config);
		p3io.endPrimer3Communications();
		TreeSet<PrimerPair> tmp = new TreeSet<PrimerPair>();
		for (PrimerPair primer : pp) {
			//logger.info("primer = " + primer);
			int leftPos = primer.getLeftPrimerPosition();
			logger.info("left primer position = " + leftPos);
			int rightPos = primer.getRightPrimerPosition();
			logger.info("right primer position = " + rightPos);
			boolean flag = true;
			for (SequenceRegion cut : cuts) {
				logger.info("crispr cut = " + cut);
				int cutPos;
				if (cut.isNegativeStrand()) {
					cutPos = cut.getStart() + 3;
				} else {
					cutPos = cut.getEnd() - 4;
				}
				int leftFrag = cutPos - leftPos;
				int rightFrag = rightPos - cutPos;
				if (Math.abs(leftFrag - rightFrag) < 50) {
					flag = false;
					logger.info("primer didn't pass test");
					break;
				}
			}
			if (flag) {
				logger.info("primer passed test!");
				tmp.add(primer);
			}
			if (tmp.size() == 2) {
				rtrn.put(cName, tmp);
				return rtrn;
			}
		}
		rtrn.put(cName, tmp);
		return rtrn;
	}
	
	private void writeLociErrors(String outFile) throws IOException {
		FileWriter w = new FileWriter(outFile);
		for (Gene locus : flaggedGenes.keySet()) {
			for (String error : flaggedGenes.get(locus)) {
				w.write(locus.getName() + "\t" + locus.toShortBED() + "\t" + error + "\n");
			}
		}
		w.close();
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws IOException, InterruptedException {
		// TODO Auto-generated method stub
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-g", "File name containing all genes", true);
		p.addStringArg("-t", "Type of the file, must be single char ('b'-bed file, 'p'-position file, 'n'-gene name file, 's'-sequence file)", true);
		p.addStringArg("-o", "output file prefix; all files created will be based on this name, only differing in suffixes", true);
		
		p.parse(args);
		
		String f = p.getStringArg("-g");
		char type = p.getStringArg("-t").charAt(0);
		String out = p.getStringArg("-o");
		
		TCSDesigner tcs = new TCSDesigner(f,type);
		tcs.writeCRISPRsWithSurvAssay(out, 4, 2);
		tcs.writeLociErrors(out + "_flaggedLoci");
		
	}
	
	public class CRISPRwriter {
		
		ArrayList<FileWriter> INwriters = new ArrayList<FileWriter>();
		ArrayList<FileWriter> PEwriters = new ArrayList<FileWriter>();
		ArrayList<FileWriter> TSSwriters = new ArrayList<FileWriter>();
		Map<String,ArrayList<FileWriter>> PrimerWriters = new TreeMap<String,ArrayList<FileWriter>>();
		Map<String,ArrayList<FileWriter>> U6primers = new TreeMap<String,ArrayList<FileWriter>>();
		Map<String,ArrayList<FileWriter>> gRNAprimers = new TreeMap<String,ArrayList<FileWriter>>();
		ArrayList<FileWriter> GblockWriters = new ArrayList<FileWriter>();
		FileWriter MasterListWriter;
		
		public CRISPRwriter(String outFile, int numCRISPRs, int numPrimers) throws IOException {
			ArrayList<FileWriter> U6in = new ArrayList<FileWriter>();
			ArrayList<FileWriter> U6pe = new ArrayList<FileWriter>();
			ArrayList<FileWriter> U6tss = new ArrayList<FileWriter>();
			ArrayList<FileWriter> gDNAin = new ArrayList<FileWriter>();
			ArrayList<FileWriter> gDNApe = new ArrayList<FileWriter>();
			ArrayList<FileWriter> gDNAtss = new ArrayList<FileWriter>();
			MasterListWriter = new FileWriter(new File(outFile + "_master_list.txt"));
			final String u6seq = "";
			final String gDNAseq = "";
			
			for (int i=1;i<=numCRISPRs; i++) {
				INwriters.add(new FileWriter(new File(outFile+"_IN_" + i)));
				PEwriters.add(new FileWriter(new File(outFile+"_PE_" + i)));
				TSSwriters.add(new FileWriter(new File(outFile+"_TSS_" + i)));
				GblockWriters.add(new FileWriter(new File(outFile+"_gBlock_" + i)));
				U6in.add(new FileWriter(new File(outFile+"_U6_IN_" + i)));
				U6pe.add(new FileWriter(new File(outFile+"_U6_PE_" + i)));
				U6tss.add(new FileWriter(new File(outFile+"_U6_TSS_" + i)));
				gDNAin.add(new FileWriter(new File(outFile+"_GRNA_IN_" + i)));
				gDNApe.add(new FileWriter(new File(outFile+"_GRNA_PE_" + i)));
				gDNAtss.add(new FileWriter(new File(outFile+"_GRNA_TSS_" + i)));
			}
			ArrayList<FileWriter> INprimers = new ArrayList<FileWriter>();
			ArrayList<FileWriter> TSSprimers = new ArrayList<FileWriter>();
			ArrayList<FileWriter> PEprimers = new ArrayList<FileWriter>();
			for (int i=1; i<=numPrimers; i++) {
				INprimers.add(new FileWriter(new File(outFile + "_surveyorPrimers_IN_" + i)));
				TSSprimers.add(new FileWriter(new File(outFile + "_surveyorPrimers_TSS_" + i)));
				PEprimers.add(new FileWriter(new File(outFile + "_surveyorPrimers_PE_" + i)));
			}
			PrimerWriters.put("IN",INprimers);
			PrimerWriters.put("TSS",TSSprimers);
			PrimerWriters.put("PE",PEprimers);
			U6primers.put("IN",U6in);
			U6primers.put("PE", U6pe);
			U6primers.put("TSS", U6tss);
			gRNAprimers.put("IN", gDNAin);
			gRNAprimers.put("PE", gDNApe);
			gRNAprimers.put("TSS", gDNAtss);
		}
		
		public void writeGene(Gene gene, Map<String,Map<Integer,SequenceRegion>> sequences, int numCRISPRs) throws IOException {
			for (String id : sequences.keySet()) {
				int i = 0;
				FileWriter cw;
				FileWriter u6;
				FileWriter grna;
				
				for (Integer crisprNum : sequences.get(id).keySet()) {
					SequenceRegion region = sequences.get(id).get(crisprNum);
					logger.info("crispr num = " + crisprNum);
					cw = getCrisprWriter(id, crisprNum);
					if (id.equals("gblock")) {
						String bases = region.getSequenceBases();
					}
					cw.write(gene.getName()+"_"+id+"_"+crisprNum+"\t"+region.getSequenceBases()+"\n");
					if (!id.equals("gblock")) {
						u6 = getU6Writer(id, crisprNum);
						//TODO: find how to modify sequence to add primer portions for u6 and gdna
						u6.write(gene.getName()+"_U6_"+id+"_"+crisprNum+"\t"+region.getAntisense().getSequenceBases()+"CGGTGTTTCGTCCTTTCCAC\n");
						grna = getgDNAWriter(id, crisprNum);
						grna.write(gene.getName()+"_gRNA_"+id+"_"+crisprNum+"\t"+region.getSequenceBases()+"GTTTTAGAGCTAGAAATAGCAAGTT\n");
					}
					i++;
				}
				while (i<4) {
					cw = getCrisprWriter(id, i+1);
					cw.write(" \t \n");
					if (!id.equals("gblock")) {
						u6 = getU6Writer(id, i+1);
						grna = getgDNAWriter(id, i+1);
						
						u6.write(" \t \n");
						grna.write(" \t \n");
					}
					i++;
				}
			}
		}
		
		public void writePrimers(Gene gene, Map<String,Collection<PrimerPair>> primerMap, int numPrimers) throws IOException {
			logger.info("primer map keyset = " + primerMap.keySet());
			for (String id : primerMap.keySet()) {
				ArrayList<FileWriter> writers = PrimerWriters.get(id);
				Collection<PrimerPair> prime = primerMap.get(id);
				int i = 0;
				logger.info("primer id = " + id);
				for (PrimerPair primer : prime) {
					logger.info("About to write: " + gene.getName()+"_surveyorPrimer_"+id+"\t"+primer.getLeftPrimer()+"\t"+primer.getRightPrimer()+"\n");
					FileWriter p = writers.get(i);
					p.write(gene.getName()+"_surveyorPrimer_"+id+"_"+(i+1)+"\t"+primer.getLeftPrimer()+"\t"+primer.getRightPrimer()+"\n");
					i++;
				}
				while (i<numPrimers) {
					FileWriter p = writers.get(i);
					p.write(" \t \n");
					i++;
				}
			}
		}
		
		private FileWriter getCrisprWriter(String id, Integer crisprNum) {
			if (id.equals("IN")) {
				return INwriters.get(crisprNum-1);
			} else if (id.equals("TSS")) {
				return TSSwriters.get(crisprNum-1);
			} else if (id.equals("PE")) {
				return PEwriters.get(crisprNum-1);
			} else if (id.equals("gblock")) {
				return GblockWriters.get(crisprNum-1);
			} else {
				throw new IllegalArgumentException("Invalid CRISPR site type (must be IN, TSS, or PE)");
			}
		}
		
		private FileWriter getU6Writer(String id, Integer crisprNum) {
			return U6primers.get(id).get(crisprNum-1);
		}
		
		private FileWriter getgDNAWriter(String id, Integer crisprNum) {
			
			return gRNAprimers.get(id).get(crisprNum-1);
		}
		
		
		public void close() throws IOException {
			for (int i = 0; i<INwriters.size(); i++) {
				INwriters.get(i).close();
				PEwriters.get(i).close();
				TSSwriters.get(i).close();
				for (String key : U6primers.keySet()) {
					U6primers.get(key).get(i).close();
					gRNAprimers.get(key).get(i).close();
				}
			}
			
			for (String key : PrimerWriters.keySet()) {
				//int size = PrimerWriters.get(key).size();
				for (FileWriter w : PrimerWriters.get(key)) {
				//for (int i = 0; i<size; i++) {
					w.close();
				}
			}
		}
	}
	
	public static class PeakReader {
		
		public static IntervalTree<Annotation> getPeaks(String file) throws IOException {
			IntervalTree<Annotation> rtrn = new IntervalTree<Annotation>();
			BufferedReader br = new BufferedReader(new FileReader(file));
			
			while (br.ready()) {
				String[] tokens = StringParser.getTokens(br.readLine());
				Annotation peak = new BasicAnnotation(tokens[0], Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]));
				rtrn.put(Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]), peak);
			}
			br.close();
			return rtrn;
		}
	}

}
