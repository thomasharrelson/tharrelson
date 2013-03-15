import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import Utils.WigReader;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.math.MannWhitney;
import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;
import broad.core.util.PipelineUtils;
import broad.pda.annotation.BEDFileParser;
import broad.pda.seq.protection.MultiSampleBindingSiteCaller;
import broad.pda.seq.protection.MultiSampleBindingSiteCaller.SampleFileParser;
import broad.pda.seq.protection.SampleData;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.feature.GeneWindow;
import nextgen.core.feature.Window;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.CountScore;
import nextgen.core.model.score.WindowScoreIterator;


public class PeakComparator {
	Map<String,Collection<Gene>> barcodePeaks;
	Map<String,Collection<Gene>> nobarcodePeaks;
	Map<String,Collection<Gene>> genesByChr;
	MultiSampleBindingSiteCaller peakCaller;
	ArrayList<SampleData> samples;
	TranscriptomeSpace tspace;
	SampleData barcodeSample;
	SampleData nobarcodeSample;
	
	public AlignmentModel barcData;
	public AlignmentModel noBarcData;

	private boolean peaksCalculated = false;
	private static Logger logger = Logger.getLogger(PeakComparator.class.getName());
	
	public PeakComparator(String listFile, String bedFile, String nbpeakFile, String bpeakFile) throws IOException, InterruptedException {
		if (nbpeakFile==null || bpeakFile==null) {
			this.peakCaller = new MultiSampleBindingSiteCaller(listFile,bedFile,30,1,10000,0.001,false);
			//MultiSampleBindingSiteCaller.SampleFileParser s = peakCaller.new SampleFileParser(listFile);
			this.samples = peakCaller.getAllSamples();
			this.genesByChr = BEDFileParser.loadDataByChr(new File(bedFile));
			Map<String,Collection<Gene>> tmpGenesByChr = new TreeMap<String, Collection<Gene>>();
			for (String chr : genesByChr.keySet()) {
				Collection<Gene> collapsedGenes = StaticMethods.collapseGenes(genesByChr.get(chr));
				
				tmpGenesByChr.put(chr, collapsedGenes);
			}
			this.genesByChr = tmpGenesByChr;
			logger.info("Done collapsing genes");
			this.tspace = new TranscriptomeSpace(genesByChr);
			
			for (SampleData sample: samples) {
				if (sample.getSampleName().contains("nobarcode")) {
					logger.info("getting nobarcode sample and model");
					logger.info(sample.getSampleName());
					this.nobarcodeSample = sample;
					noBarcData = new AlignmentModel(nobarcodeSample.getSampleName() + ".bam",tspace,true);
				} else {
					this.barcodeSample = sample;
					barcData = new AlignmentModel(barcodeSample.getSampleName() + ".bam", tspace, true);;
				}
				Map<String,Collection<Gene>> peaks = getPeaks(sample);
				if (sample.getSampleName().contains("nobarcode")) {
					this.nobarcodePeaks = peaks;
				} else {
					this.barcodePeaks = peaks;
				}
			}
			if (!barcodePeaks.isEmpty() && !nobarcodePeaks.isEmpty()) {
				peaksCalculated = true;
			}
			if (barcData == null || noBarcData == null) {
				throw new IllegalArgumentException("Did not instantiate both barcode and nonbarcode alignment models; check the sample list file");
			}
			
			return;
		} else if (bpeakFile!=null && nbpeakFile!=null) {
			this.peakCaller = new MultiSampleBindingSiteCaller(listFile,bedFile,30,1,10000,0.001,false);
			//MultiSampleBindingSiteCaller.SampleFileParser s = peakCaller.new SampleFileParser(listFile);
			this.samples = peakCaller.getAllSamples();
			this.genesByChr = BEDFileParser.loadDataByChr(new File(bedFile));
			Map<String,Collection<Gene>> tmpGenesByChr = new TreeMap<String, Collection<Gene>>();
			for (String chr : genesByChr.keySet()) {
				Collection<Gene> collapsedGenes = StaticMethods.collapseGenes(genesByChr.get(chr));
				
				tmpGenesByChr.put(chr, collapsedGenes);
			}
			this.genesByChr = tmpGenesByChr;
			logger.info("Done collapsing genes");
			this.tspace = new TranscriptomeSpace(genesByChr);
			
			for (SampleData sample: samples) {
				if (sample.getSampleName().contains("nobarcode")) {
					logger.info("getting nobarcode sample and model");
					logger.info(sample.getSampleName());
					this.nobarcodeSample = sample;
					noBarcData = new AlignmentModel(nobarcodeSample.getSampleName() + ".bam",tspace,true);
				} else {
					this.barcodeSample = sample;
					barcData = new AlignmentModel(barcodeSample.getSampleName() + ".bam", tspace, true);;
				}
				
				if (sample.getSampleName().contains("nobarcode")) {
					Map<String,Collection<Gene>> peaks = BEDFileParser.loadDataByChr(new File(nbpeakFile));
					this.nobarcodePeaks = peaks;
				} else {
					Map<String,Collection<Gene>> peaks = BEDFileParser.loadDataByChr(new File(bpeakFile));
					this.barcodePeaks = peaks;
				}
			}
			if (!barcodePeaks.isEmpty() && !nobarcodePeaks.isEmpty()) {
				peaksCalculated = true;
			}
			if (barcData == null || noBarcData == null) {
				throw new IllegalArgumentException("Did not instantiate both barcode and nonbarcode alignment models; check the sample list file");
			}
		}
		/*
		try {
			calculatePeaks(listFile,bedFile);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}*/
		
	}
	
	public Map<String,Collection<Gene>> getPeaks(SampleData sample) throws IOException, InterruptedException {
		if (peaksCalculated) {
			if (sample.getSampleName().contains("nobarcode")) {
				return nobarcodePeaks;
			} else {
				return barcodePeaks;
			}
		} else {
			Map<String,Collection<Gene>> rtrn = new TreeMap<String,Collection<Gene>>();
			for (String chr : genesByChr.keySet()) {
				Collection<Gene> genes = genesByChr.get(chr);
				for (Gene g : genes) {
					logger.info(sample.getSampleName());
					logger.info("Gene: " + g);
					Collection<Annotation> peaks = peakCaller.getSingleSampleScanPeaks(sample, g);
					Collection<Gene> genePeaks = convert2Genes(peaks);
					logger.info("size of collection pre-collapsing: " + genePeaks.size());
					Collection<Gene> collPeaks = StaticMethods.collapseGenes(genePeaks);
					logger.info("size of collapsed peaks collection: " + collPeaks.size());
					if (rtrn.containsKey(chr)) {
						Collection<Gene> existingData = rtrn.get(chr);
						existingData.addAll(collPeaks);
						rtrn.put(chr, existingData);
					} else {
						rtrn.put(chr,collPeaks);
					}
				}
			}
			return rtrn;
		}
	}
	
	private Collection<Gene> convert2Genes(Collection<Annotation> annotations) {
		Collection<Gene> rtrn = new TreeSet<Gene>();
		for (Annotation anno : annotations) {
			rtrn.add(new Gene(anno));
		}
		return rtrn;
	}
	
	
	
	public SampleData getBarcodeSample() {
		return barcodeSample;
	}
	
	public SampleData getNobarcodeSample() {
		return nobarcodeSample;
	}
	
	public Map<String,Map<Integer,Double>> AggregatePeaks(Map<String,Collection<Gene>> bPeakMap, Map<String,Collection<Gene>> nbPeakMap, String bWig, String nbWig, boolean norm) throws IOException {
		// this only works if passed the nonbarcoded peaks, because we only want to look at peaks associated w/ nonbarcode sample
		Map<String,Map<Integer,Double>> rtrn = new TreeMap<String,Map<Integer,Double>>();
		Map<Integer,Double> totalBarPeakCount = new TreeMap<Integer,Double>();
		Map<Integer,Double> totalNobarPeakCount = new TreeMap<Integer,Double>();
		WigReader bwr = null;
		WigReader nbwr = null;
		if (bWig != null && nbWig != null) {
			bwr = new WigReader(bWig);
			nbwr = new WigReader(nbWig);
		}
		int numNBpeaks=0;
		for (String chr : bPeakMap.keySet()) {
			Collection<Gene> bPeaks = bPeakMap.get(chr);
			Collection<Gene> nbPeaks = nbPeakMap.get(chr);
			IntervalTree<Gene> barTree = StaticMethods.makeTree(bPeaks);
			
			numNBpeaks += nbPeaks.size();
			
			for (Gene nbPeakWindow : nbPeaks) {
				Iterator<Node<Gene>> ovrItr = barTree.overlappers(nbPeakWindow.getStart(), nbPeakWindow.getEnd());
				ArrayList<Gene> overlaps = new ArrayList<Gene>();
				Gene nextPeak;
				while (ovrItr.hasNext()) {
					// what happens if multiple barcoded peaks overlap the nonbarcoded peak in question?
					nextPeak = ovrItr.next().getValue();
					if (nbPeakWindow.overlapsStranded(nextPeak)) {
						overlaps.add(nextPeak);
					}
				}
				if (overlaps.size()>1) {
					// this is a problem
					// should try to define the best match by determining the size of the intersections
					numNBpeaks--;
					continue;
				} else if (overlaps.size() == 0) {
					numNBpeaks--;
					continue;
				}
				Annotation bPeakWindow = overlaps.get(0);
				// find window to generate counts over.
				int startPos;
				boolean negative = false;
				Collection<? extends Window> checkWindow;
				if (nbPeakWindow.isNegativeStrand()) {
					negative = true;
					
					startPos = bPeakWindow.getEnd()+25;
					checkWindow = tspace.getFragment(nbPeakWindow.getChr(),startPos-50,startPos);
				} else {
					startPos = bPeakWindow.getStart()-25;
					checkWindow = tspace.getFragment(nbPeakWindow.getChr(),startPos,startPos+50);
				}
				// gotta find counts using alignment model or wig files
				Peak bPeak;
				Peak nbPeak;
				
				// gotta convert peaks to prob distributions
				int shift;
				if (nbPeakWindow.isNegativeStrand()) {
					shift = bPeakWindow.getEnd() - nbPeakWindow.getEnd();
				} else {
					shift = nbPeakWindow.getStart() - bPeakWindow.getStart();
				}
				if (shift<0) {
					if (!negative) {
						nbPeakWindow = nbPeakWindow.trimAbsolute(bPeakWindow.getStart(), nbPeakWindow.getEnd());
					} else {
						nbPeakWindow = nbPeakWindow.trimAbsolute(nbPeakWindow.getStart(), bPeakWindow.getEnd());
					}
					shift = 0;
				}
				if (bWig == null) {
					bPeak = new Peak(bPeakWindow, barcData);
				} else {
					bPeak = new Peak(bPeakWindow, bwr);
				}
				if (nbWig == null) {
					nbPeak = new Peak(nbPeakWindow, noBarcData);
				} else {
					nbPeak = new Peak(nbPeakWindow, nbwr);
				}
				logger.info("bPeak: " + bPeak);
				logger.info("nbPeak: " + nbPeak);
				Map<Integer,Double> bPeakProbDist = null;
				Map<Integer,Double> nbPeakProbDist = null;
				if (norm) {
					bPeakProbDist = convertCountsToProbDist(bPeak,0,true); // since we are setting the left end of window to the bPeak
					nbPeakProbDist = convertCountsToProbDist(nbPeak,shift,false);
				} else {
					Map<Integer,Double> bPeakPosCountMap = bPeak.getPosCountMap();
					Map<Integer,Double> nbPeakPosCountMap = nbPeak.getPosCountMap();
					bPeakProbDist = new TreeMap<Integer,Double>();
					nbPeakProbDist = new TreeMap<Integer,Double>();
					for (int genomePos : bPeakPosCountMap.keySet()) {
						bPeakProbDist.put(bPeak.getPositionAtReferenceCoordinate(genomePos, false), bPeakPosCountMap.get(genomePos));
					}
					for (int genomePos : nbPeakPosCountMap.keySet()) {
						nbPeakProbDist.put(nbPeak.getPositionAtReferenceCoordinate(genomePos, false)+shift, nbPeakPosCountMap.get(genomePos));
					}
				}
				// Consider turning following for loops into a method
				// They are redundant
				logger.info("bar prob dist size: " + bPeakProbDist.size());
				logger.info("nobar prob dist size: " + nbPeakProbDist.size());
				for (int pos : bPeakProbDist.keySet()) {
					if (totalBarPeakCount.containsKey(pos)) {
						double total = bPeakProbDist.get(pos) + totalBarPeakCount.get(pos);
						totalBarPeakCount.put(pos, total);
						logger.info("added value at pos: " + pos +", value = " + totalBarPeakCount.get(pos));
					} else {
						totalBarPeakCount.put(pos, bPeakProbDist.get(pos));
						logger.info("new value at pos: " + pos +", value = " + totalBarPeakCount.get(pos));
					}
				}
				for (int pos : nbPeakProbDist.keySet()) {
					if (totalNobarPeakCount.containsKey(pos)) {
						double total = nbPeakProbDist.get(pos) + totalNobarPeakCount.get(pos);
						totalNobarPeakCount.put(pos, total);
						logger.info("added value at pos: " + pos +", value = " + totalNobarPeakCount.get(pos));
					} else {
						totalNobarPeakCount.put(pos, nbPeakProbDist.get(pos));
						logger.info("new value at pos: " + pos +", value = " + totalNobarPeakCount.get(pos));
					}
				}
				logger.info("total number of peaks covered so far: " + numNBpeaks);
			}
			
			
		}
		// Complete finding the average for counts
		for(int pos : totalBarPeakCount.keySet()) {
			double curr = totalBarPeakCount.get(pos);
			totalBarPeakCount.put(pos, (curr/numNBpeaks));
		}
		for(int pos : totalNobarPeakCount.keySet()) {
			double curr = totalNobarPeakCount.get(pos);
			totalNobarPeakCount.put(pos, (curr/numNBpeaks));
		}
		rtrn.put("nobarcode",totalNobarPeakCount);
		rtrn.put("barcode",totalBarPeakCount);
		return rtrn;
		
	}
	
	public Map<String,Map<Integer,Double>> AggregatePeaksFixedNBwindow(Map<String,Collection<Gene>> bPeakMap, Map<String,Collection<Gene>> nbPeakMap, String bWig, String nbWig, int relStart, int winSize, boolean norm, int thresh) throws IOException {
		// this only works if passed the nonbarcoded peaks, because we only want to look at peaks associated w/ nonbarcode sample
		/* This method will create a fixed window of selected size and relative position.
		 * It will be fixed to the relative position of the nobarcode peak.
		 * I am initially going to remove the overlap check portion bc I want to iterate over all nb peaks regardless
		 * of what is happening in the barcode sample. The overlap part will be commented out so that it can be added back in later
		 * if needed.
		 */
		Map<String,Map<Integer,Double>> rtrn = new TreeMap<String,Map<Integer,Double>>();
		Map<Integer,Double> totalBarPeakCount = new TreeMap<Integer,Double>();
		Map<Integer,Double> totalNobarPeakCount = new TreeMap<Integer,Double>();
		WigReader bwr = null;
		WigReader nbwr = null;
		if (bWig != null && nbWig != null) {
			bwr = new WigReader(bWig);
			nbwr = new WigReader(nbWig);
		}
		int numNBpeaks=0;
		for (String chr : bPeakMap.keySet()) {
			Collection<Gene> bPeaks = bPeakMap.get(chr);
			Collection<Gene> nbPeaks = nbPeakMap.get(chr);
			//IntervalTree<Gene> barTree = StaticMethods.makeTree(bPeaks);
			
			numNBpeaks += nbPeaks.size();
			
			for (Gene nbPeakWindow : nbPeaks) {
				
				//Annotation bPeakWindow = tspace.;
				// find window to generate counts over.
				int startPos;
				boolean negative = false;
				Collection<? extends Window> checkWindow;
				logger.info("nbPeakWindow: " + nbPeakWindow);
				int peakSize = nbPeakWindow.size();
				if (nbPeakWindow.isNegativeStrand()) {
					negative = true;
					if (peakSize >= relStart+winSize) {
						checkWindow = tspace.getFragment(nbPeakWindow.getChr(),nbPeakWindow.getReferenceCoordinateAtPosition(relStart+winSize-1,false),nbPeakWindow.getEnd()-relStart);
					} else {
						checkWindow = tspace.getFragment(nbPeakWindow.getChr(),nbPeakWindow.getStart() - (relStart+winSize-peakSize),nbPeakWindow.getEnd()-relStart);
					}
				} else {
					if (peakSize >= relStart+winSize) {
						checkWindow = tspace.getFragment(nbPeakWindow.getChr(),nbPeakWindow.getStart()+relStart,nbPeakWindow.getReferenceCoordinateAtPosition(relStart+winSize-1,false)+1);
					} else {
						checkWindow = tspace.getFragment(nbPeakWindow.getChr(),nbPeakWindow.getStart()+relStart,nbPeakWindow.getEnd() + relStart+winSize-peakSize);
					}
				}
				
				if (checkWindow.size() > 1) {
					checkWindow = StaticMethods.checkOrientation(checkWindow, nbPeakWindow); // getFragment doesn't consider orientation... boo
				} else if (checkWindow.size() == 0) {
					logger.warn("checkWindow size = 0; this is due to the peak having no corresponding fragment in tspace");
					continue;
				}
				Annotation newWindow = null;
				if (checkWindow.size()>1) {
					// temp fix
					
					for (Annotation win : checkWindow) {
						if (newWindow==null) {
							if (nbPeakWindow.overlapsStranded(win))
								newWindow = new Gene(win);
						} else {
							// can't take union: this is a bug. need to find which is the real window
							// can be done through analyzing intersections or overlapStranded method
							// overlapStranded works well with a collapsed tspace, but will have problems with a regular tspace
							if (nbPeakWindow.overlapsStranded(win)) {
								if (nbPeakWindow.intersect(win).size()>nbPeakWindow.intersect(newWindow).size()) {
									newWindow = new Gene(win);
								}
							}
						}
					}
					logger.warn("About to throw exception, here is checkWindow: " + checkWindow);
					//throw new IllegalArgumentException("check window collection size is still greater than one. Consider checking whether current nonbarcode peak window is oriented");
				} else if  (checkWindow.size()==1) {
					newWindow = checkWindow.iterator().next();
				}
				// gotta find counts using alignment model or wig files
				Peak bPeak;
				Peak nbPeak;
				
				if (bWig == null) {
					bPeak = new Peak(newWindow, barcData);
				} else {
					bPeak = new Peak(newWindow, bwr);
				}
				if (nbWig == null) {
					nbPeak = new Peak(newWindow, noBarcData);
				} else {
					nbPeak = new Peak(newWindow, nbwr);
				}
				logger.info("bPeak: " + bPeak);
				logger.info("nbPeak: " + nbPeak);
				// check barcode peak distribution count (can change the criteria for a viable count distribution)
				
				if (bPeak.getTotalCount()<=thresh || nbPeak.getTotalCount()<=thresh) {
					numNBpeaks = numNBpeaks - 1;
					logger.info("numNBpeaks = " + numNBpeaks);
					continue;
				}
				
				// need to find shift bc the peaks are found using getFragment in tspace, so they may be trimmed and are not all 50 bp.
				// shift is required to set the 0 point to the start of the nonbarcode peak; this must be true for all peaks
				int zeroPoint;
				int shift=0;
				int fivePrimePeak;

				if (negative) {
					zeroPoint = nbPeakWindow.getEnd();
					fivePrimePeak = bPeak.getEnd();
					shift = zeroPoint - fivePrimePeak;
				} else {
					zeroPoint = nbPeakWindow.getStart();
					fivePrimePeak = bPeak.getStart();
					shift = fivePrimePeak - zeroPoint;
				}
				Map<Integer,Double> bPeakProbDist = null;
				Map<Integer,Double> nbPeakProbDist = null;
				if (norm) {
					bPeakProbDist = convertCountsToProbDist(bPeak,shift,true); 
					nbPeakProbDist = convertCountsToProbDist(nbPeak,shift,false);
				} else {
					Map<Integer,Double> bPeakPosCountMap = bPeak.getPosCountMap();
					Map<Integer,Double> nbPeakPosCountMap = nbPeak.getPosCountMap();
					bPeakProbDist = new TreeMap<Integer,Double>();
					nbPeakProbDist = new TreeMap<Integer,Double>();
					for (int genomePos : bPeakPosCountMap.keySet()) {
						bPeakProbDist.put(bPeak.getPositionAtReferenceCoordinate(genomePos, false)+shift, bPeakPosCountMap.get(genomePos));
						nbPeakProbDist.put(bPeak.getPositionAtReferenceCoordinate(genomePos, false)+shift, nbPeakPosCountMap.get(genomePos));
					}
				}
				// Consider turning following for loops into a method
				// They are redundant
				logger.info("bar prob dist size: " + bPeakProbDist.size());
				logger.info("nobar prob dist size: " + nbPeakProbDist.size());
				logger.info("bar prob dist: " + bPeakProbDist);
				logger.info("nobar prob dist: " + nbPeakProbDist);
				for (int pos : bPeakProbDist.keySet()) {
					if (totalBarPeakCount.containsKey(pos)) {
						double total = bPeakProbDist.get(pos) + totalBarPeakCount.get(pos);
						totalBarPeakCount.put(pos, total);
						logger.info("added value at pos: " + pos +", value = " + totalBarPeakCount.get(pos));
					} else {
						totalBarPeakCount.put(pos, bPeakProbDist.get(pos));
						logger.info("new value at pos: " + pos +", value = " + totalBarPeakCount.get(pos));
					}
				}
				for (int pos : nbPeakProbDist.keySet()) {
					if (totalNobarPeakCount.containsKey(pos)) {
						double total = nbPeakProbDist.get(pos) + totalNobarPeakCount.get(pos);
						totalNobarPeakCount.put(pos, total);
						logger.info("added value at pos: " + pos +", value = " + totalNobarPeakCount.get(pos));
					} else {
						totalNobarPeakCount.put(pos, nbPeakProbDist.get(pos));
						logger.info("new value at pos: " + pos +", value = " + totalNobarPeakCount.get(pos));
					}
				}
				logger.info("total number of peaks covered so far: " + numNBpeaks);
			}
			
		}
		// Complete finding the average for counts
		for(int pos : totalBarPeakCount.keySet()) {
			double curr = totalBarPeakCount.get(pos);
			totalBarPeakCount.put(pos, (curr/numNBpeaks));
		}
		for(int pos : totalNobarPeakCount.keySet()) {
			double curr = totalNobarPeakCount.get(pos);
			totalNobarPeakCount.put(pos, (curr/numNBpeaks));
		}
		rtrn.put("nobarcode",totalNobarPeakCount);
		rtrn.put("barcode",totalBarPeakCount);
		return rtrn;
		
	}
	
	public double getMedian(Peak peak, int shift) {
		Map<Integer,Double> countMap = peak.getPosCountMap();
		double total = peak.getTotalCount();
		double medianNum = (total+1)/2;
		double curr = 0;
		double rtrn = 0;
		for (int i = 0; i<peak.size(); i++) {
			int relPos = peak.getReferenceCoordinateAtPosition(i, false);
			if (countMap.containsKey(relPos)) {
				curr+=countMap.get(relPos);
				if (curr>=Math.floor(medianNum)) {
					if ((total+1) % 2 != 0) {
						if (curr==Math.floor(medianNum)) {
							if (rtrn!=0) {
								rtrn = (rtrn+i + shift)/2;
								break;
							} else {
								rtrn = i + shift;
							}
						} else {
							rtrn = i + shift;
							break;
						}
					} else {
						rtrn = i + shift;
						break;
					}
				}
			}
		}
		return rtrn;
	}
	
	public double getMean(Peak peak, int shift) {
		double rtrn = 0;
		Map<Integer,Double> countMap = peak.getPosCountMap();
		
		for (int i=0; i<peak.size(); i++) {
			int refPos = peak.getReferenceCoordinateAtPosition(i, false);
			if (countMap.containsKey(refPos)) {
				rtrn += (i+shift)*countMap.get(refPos);
			}
		}
		
		return rtrn/peak.getTotalCount();
	}
	
	public double getVariance(Peak peak, double mean, int shift) {
		double rtrn = 0;
		Map<Integer,Double> countMap = peak.getPosCountMap();
		for (int pos : countMap.keySet()) {
			int relPos = peak.getPositionAtReferenceCoordinate(pos,false);
			relPos = relPos + shift;
			double diffSq = Math.pow(relPos - mean, 2);
			rtrn += diffSq*countMap.get(pos);
		}
		return rtrn/peak.getTotalCount();
	}
	
	public double findDiffBtPeaks(Peak peak1, Peak peak2, int shift) {
		double mean1 = getMean(peak1, shift);
		double mean2 = getMean(peak2, shift);
		double var1 = getVariance(peak1, mean1, shift);
		double var2 = getVariance(peak2, mean2, shift);
		
		return (mean1-mean2)/(Math.sqrt(var1+var2));
	}
	
	public int findShift(Gene nbPeak, Peak nbPeakWindow) {
		int zeroPoint;
		int shift=0;
		int fivePrimePeak;
		boolean negative = nbPeak.isNegativeStrand();

		if (negative) {
			zeroPoint = nbPeak.getEnd();
			fivePrimePeak = nbPeakWindow.getEnd();
			shift = zeroPoint - fivePrimePeak;
		} else {
			zeroPoint = nbPeak.getStart();
			fivePrimePeak = nbPeakWindow.getStart();
			shift = fivePrimePeak - zeroPoint;
		}
		return shift;
	}
	
	public double findTotalDiffBtPeaks(Map<String,Collection<Gene>> bPeakMap, Map<String,Collection<Gene>> nbPeakMap, int relStart, int winSize, String nbwig, String bwig, int thresh) throws IOException {
		double rtrn = 0;
		WigReader bwr = null;
		WigReader nbwr = null;
		if (nbwig!=null) {
			nbwr = new WigReader(nbwig);
		}
		if (bwig!=null) {
			bwr = new WigReader(bwig);
		}
		for (String chr : nbPeakMap.keySet()) {
			Collection<Gene> nbPeaks = nbPeakMap.get(chr);
			for (Gene nbPeak : nbPeaks) {
				Peak nbPeakWindow = createPeakWindow(nbPeak, noBarcData, relStart, winSize, nbwr);
				Peak bPeakWindow = createPeakWindow(nbPeak, barcData, relStart, winSize, bwr);
				if (nbPeakWindow==null || bPeakWindow==null) { continue;}
				int shift = findShift(nbPeak,nbPeakWindow);
				if (nbPeakWindow.getTotalCount()>=thresh && bPeakWindow.getTotalCount()>=thresh) {
					rtrn += findDiffBtPeaks(nbPeakWindow, bPeakWindow, shift);
				}
			}
		}
		return rtrn;
	}
	
	public Peak createPeakWindow(Gene nbPeakWindow, AlignmentModel data, int relStart, int winSize, WigReader wr) throws IOException {
		
		Collection<? extends Window> checkWindow;
		boolean negative = false;
		logger.info("nbPeakWindow: " + nbPeakWindow);
		int peakSize = nbPeakWindow.size();
		if (nbPeakWindow.isNegativeStrand()) {
			negative = true;
			if (peakSize >= relStart+winSize) {
				checkWindow = tspace.getFragment(nbPeakWindow.getChr(),nbPeakWindow.getReferenceCoordinateAtPosition(relStart+winSize-1,false),nbPeakWindow.getEnd()-relStart);
			} else {
				checkWindow = tspace.getFragment(nbPeakWindow.getChr(),nbPeakWindow.getStart() - (relStart+winSize-peakSize),nbPeakWindow.getEnd()-relStart);
			}
		} else {
			if (peakSize >= relStart+winSize) {
				checkWindow = tspace.getFragment(nbPeakWindow.getChr(),nbPeakWindow.getStart()+relStart,nbPeakWindow.getReferenceCoordinateAtPosition(relStart+winSize-1,false)+1);
			} else {
				checkWindow = tspace.getFragment(nbPeakWindow.getChr(),nbPeakWindow.getStart()+relStart,nbPeakWindow.getEnd() + relStart+winSize-peakSize);
			}
		}
		
		if (checkWindow.size() > 1) {
			checkWindow = StaticMethods.checkOrientation(checkWindow, nbPeakWindow); // getFragment doesn't consider orientation... boo
		} else if (checkWindow.size() == 0) {
			logger.warn("checkWindow size = 0; this is due to the peak having no corresponding fragment in tspace");
			return null;
		}
		Annotation newWindow = null;
		if (checkWindow.size()>1) {
			for (Annotation win : checkWindow) {
				if (newWindow==null) {
					if (nbPeakWindow.overlapsStranded(win))
						newWindow = new Gene(win);
				} else {
					if (nbPeakWindow.overlapsStranded(win)) {
						if (nbPeakWindow.intersect(win).size()>nbPeakWindow.intersect(newWindow).size()) {
							newWindow = new Gene(win);
						}
					}
				}
			}
			logger.warn("About to throw exception, here is checkWindow: " + checkWindow);
			logger.warn("Here is the new window: " + newWindow);
			//throw new IllegalArgumentException("check window collection size is still greater than one. Consider checking whether current nonbarcode peak window is oriented");
		} else if  (checkWindow.size()==1) {
			newWindow = checkWindow.iterator().next();
		}
		// gotta find counts using alignment model or wig files
		Peak peak;
		
		if (wr == null) {
			peak = new Peak(newWindow, data);
		} else {
			peak = new Peak(newWindow, wr);
		}
		return peak;
	}
	
	public Map<Annotation, Map<Double,Double>> getAllPeakMedianPairs(Map<String,Collection<Gene>> bPeakMap, Map<String,Collection<Gene>> nbPeakMap, int relStart, int winSize, String nbwig, String bwig, int thresh) throws IOException {
		Map<Annotation, Map<Double,Double>> rtrn = new TreeMap<Annotation, Map<Double,Double>>();
		WigReader bwr = null;
		WigReader nbwr = null;
		if (nbwig!=null) {
			nbwr = new WigReader(nbwig);
		}
		if (bwig!=null) {
			bwr = new WigReader(bwig);
		}
		for (String chr : nbPeakMap.keySet()) {
			Collection<Gene> nbPeaks = nbPeakMap.get(chr);
			for (Gene nbPeak : nbPeaks) {
				Peak nbPeakWindow = createPeakWindow(nbPeak, noBarcData, relStart, winSize, nbwr);
				Peak bPeakWindow = createPeakWindow(nbPeak, barcData, relStart, winSize, bwr);
				if (nbPeakWindow==null || bPeakWindow== null) { continue;}
				int shift = findShift(nbPeak, nbPeakWindow);
				if (nbPeakWindow.getTotalCount()>=thresh && bPeakWindow.getTotalCount()>=thresh) {
					double key = getMedian(nbPeakWindow, shift);
					double val = getMedian(bPeakWindow, shift);
					Map<Double,Double> tmp = new TreeMap<Double,Double>();
					tmp.put(key, val);
					rtrn.put(nbPeakWindow, tmp);
				}
			}
		}
		return rtrn;
	}
	
	public Map<String,Double> calculateAvgPeakSize(Map<String,Collection<Gene>> bPeakMap, Map<String,Collection<Gene>> nbPeakMap) {
		Map<String,Double> rtrn = new TreeMap<String,Double>();
		double totalPeaks = 0;
		double nbTotal = 0;
		double bTotal = 0;
		logger.info("About to find average peak size...");
		logger.info("bPeakMap size: " + bPeakMap.size());
		logger.info("nbPeakMap size: " + nbPeakMap.size());
		for (String chr : nbPeakMap.keySet()) {	
			logger.info("Current chromosome: " + chr);
			Collection<Gene> bPeaks = bPeakMap.get(chr);
			logger.info(chr + "bpeak collection size: " + bPeaks.size());
			Collection<Gene> nbPeaks = nbPeakMap.get(chr);
			logger.info(chr + "nbpeak collection size: " + nbPeaks.size());
			IntervalTree<Gene> tree = StaticMethods.makeTree(bPeaks);
			double bSize;
			double nbSize;
			
			totalPeaks += nbPeaks.size();
			
			for (Gene nbPeak : nbPeaks) {
				nbSize = 0;
				bSize = 0;
				Collection<Gene> overlapBpeaks = findOverlapPeaks(nbPeak, tree);
				logger.info("overlap collection size: " + overlapBpeaks.size());
				logger.info("Peak: " + nbPeak.toBED());
				nbSize = nbPeak.size();
				if (overlapBpeaks.size()==1) {
					
					Annotation overlap = overlapBpeaks.iterator().next();
					bSize = overlap.size();
				} else {
					totalPeaks--;
					continue;}
				bTotal += bSize;
				nbTotal += nbSize;
			}
			//bTotal = bTotal/totalPeaks;
			//nbTotal = nbTotal/totalPeaks;	
		}
		nbTotal = nbTotal/totalPeaks;
		bTotal = bTotal/totalPeaks;
		rtrn.put("nobarcode", nbTotal);
		rtrn.put("barcode", bTotal);
		rtrn.put("total number of peaks averaged", totalPeaks);
		return rtrn;
	}
	
	private Collection<Gene> findOverlapPeaks(Gene peak1, Collection<Gene> otherPeaks) {
		TreeSet<Gene> rtrn = new TreeSet<Gene>();
		IntervalTree<Gene> tree = StaticMethods.makeTree(otherPeaks);
		Iterator<Node<Gene>> bPeakItr = tree.overlappers(peak1.getStart(),peak1.getEnd());

		while (bPeakItr.hasNext()) {
			rtrn.add(bPeakItr.next().getValue());
		}
		return rtrn;
	}
	
	private Collection<Gene> findOverlapPeaks(Gene peak1, IntervalTree<Gene> tree) {
		TreeSet<Gene> rtrn = new TreeSet<Gene>();
		
		Iterator<Node<Gene>> bPeakItr = tree.overlappers(peak1.getStart(),peak1.getEnd());

		while (bPeakItr.hasNext()) {
			rtrn.add(bPeakItr.next().getValue());
		}
		return rtrn;
	}
	
	private Map<Integer,Double> normalizeCountsByGeneAvg(Peak peak, int shift, boolean isbarcode) {
		Map<Integer,Double> rtrn = new TreeMap<Integer,Double>();
		Map<Integer,Double> peakMap = peak.getPosCountMap();
		double totalCount = peak.getTotalCount();
		Iterator<? extends Window> ref = tspace.getOverlappingRegion(peak.getChr(), peak.getStart(), peak.getEnd()).iterator();
		
		Annotation refGene = null;
		while (ref.hasNext()) {
			Collection<? extends Annotation> genes = ref.next().getSourceAnnotations();
			refGene = genes.iterator().next();
			if (refGene != null) {
				if (peak.overlaps(refGene)) {
					break; // beware of snoRNA's; I can't figure out a way to filter them out. It's also possible that a peak is called at a snoRNA site which doesn't overlap anything else.
				}
			}
			refGene = null;
		}
		if (refGene == null) {
			throw new IllegalArgumentException("refGene is still null after iterating over overlapping parts; must not have any overlappers");
		}
		double normCount;
		if (isbarcode) {
			normCount = barcData.getCount(refGene,true);
		} else {
			normCount = noBarcData.getCount(refGene,true);
		}
		if (refGene.size()==0 || normCount==0) {
			logger.warn("refGene has no counts; refGene is = " + refGene + "\n normCount = " + normCount + "\n peak is " + peak + "\n peak counts = " + peak.getTotalCount());
			normCount = 1; // so you're not dividing by 0.
			//throw new IllegalArgumentException("either count or size for reference gene is equal to zero: DEBUG");
		}
		logger.info("Checking shift value (make sure it makes sense): " + shift);
		for (int i = 0; i<peak.size();i++) {
			int refPos = peak.getReferenceCoordinateAtPosition(i, false);
			if (peakMap.containsKey(refPos)) {
				rtrn.put(i+shift, (peakMap.get(refPos)/(normCount/refGene.size())));
			}
		}
		return rtrn;
	}
	
	private Map<Integer,Double> convertCountsToProbDist(Peak peak, int shift, boolean isbarcode) {
		Map<Integer,Double> rtrn = new TreeMap<Integer,Double>();
		Map<Integer,Double> peakMap = peak.getPosCountMap();
		double totalCount = peak.getTotalCount();
		/*
		Iterator<? extends Window> ref = tspace.getOverlappingRegion(peak.getChr(), peak.getStart(), peak.getEnd()).iterator();
		
		Annotation refGene = null;
		while (ref.hasNext()) {
			Collection<? extends Annotation> genes = ref.next().getSourceAnnotations();
			refGene = genes.iterator().next();
			if (refGene != null) {
				if (peak.overlaps(refGene)) {
					break; // beware of snoRNA's; I can't figure out a way to filter them out. It's also possible that a peak is called at a snoRNA site which doesn't overlap anything else.
				}
			}
			refGene = null;
		}
		if (refGene == null) {
			throw new IllegalArgumentException("refGene is still null after iterating over overlapping parts; must not have any overlappers");
		}
		double normCount;
		if (isbarcode) {
			normCount = barcData.getCount(refGene,true);
		} else {
			normCount = noBarcData.getCount(refGene,true);
		}
		if (refGene.size()==0 || normCount==0) {
			logger.warn("refGene has no counts; refGene is = " + refGene + "\n normCount = " + normCount + "\n peak is " + peak + "\n peak counts = " + peak.getTotalCount());
			normCount = 1; // so you're not dividing by 0.
			//throw new IllegalArgumentException("either count or size for reference gene is equal to zero: DEBUG");
		}*/
		logger.info("Checking shift value (make sure it makes sense): " + shift);
		for (int i = 0; i<peak.size();i++) {
			int refPos = peak.getReferenceCoordinateAtPosition(i, false);
			if (peakMap.containsKey(refPos)) {
				rtrn.put(i+shift, (peakMap.get(refPos)/(totalCount)));
			}
		}
		return rtrn;
	}
	
	public void writeScoredPeaks(String outfile, Map<String,Collection<Gene>> nbPeakMap, int relStart, int winSize, String bwig, String nbwig, int thresh) throws IOException {
		// this method will use of set of unified peaks identified from the nonbarcode sample (or some other file)
		// thus it doesn't need the barcode peaks, only needs the barcode alignment model.
		WigReader bwr = null;
		WigReader nbwr = null;
		if (bwig!=null) {
			bwr = new WigReader(bwig);
		}
		if (nbwig!=null) {
			nbwr = new WigReader(nbwig);
		}
		TreeMap<Peak,Double> scoreMap = new TreeMap<Peak,Double>();
		for (String chr : nbPeakMap.keySet()) {
			Collection<Gene> nbPeaks = nbPeakMap.get(chr);
			
			for (Gene nbPeak : nbPeaks) {
				Peak nbPeakWindow = createPeakWindow(nbPeak, noBarcData, relStart, winSize, nbwr);
				Peak bPeakWindow = createPeakWindow(nbPeak, barcData, relStart, winSize, bwr);
				if (nbPeakWindow==null || bPeakWindow==null) { continue;}
				int shift = findShift(nbPeak,nbPeakWindow);
				if (nbPeakWindow.getTotalCount()>=thresh && bPeakWindow.getTotalCount()>=thresh) {
					MannWhitney mw = getMannWhitney(bPeakWindow, nbPeakWindow, shift);
					double score = mw.getPvalue();
					nbPeakWindow.setScore(score);
					
					logger.info(nbPeakWindow + "score: " + score);
					logger.info("test statistic: " + mw.getZ());
					Double Z = new Double(mw.getZ());
					if (Z.doubleValue() < 0 && (!Z.isNaN())) {
						if (score!=0) {
							score = -1*score;
						} else {
							score = -1*Math.pow(10, -15);
						}
					}
					scoreMap.put(nbPeakWindow, score);
				}
			}
		}
		FileWriter w = new FileWriter(new File(outfile));
		for (Peak peak : scoreMap.keySet()) {
			double pval = scoreMap.get(peak);
			double color = convertPvalToColor(pval);
			if (pval<0) {
				w.write(peak.toBED(true, 0, 0, (int) color) + "\n");
			} else {
				w.write(peak.toBED(true, (int) color, 0, 0) + "\n");
			}
		}
		w.close();
	}
	
	private double convertPvalToColor(Double pval) {
		if (!pval.isNaN() && pval!=0) {
			double rtrn = (-83*Math.log10(Math.abs(pval)));
			if (rtrn > 254) {
				return 254;}
			else { return rtrn;}
		} else if (pval.isNaN()) {
			return 0.0;
		} else {
			return 254;
		}
	}
	
	public MannWhitney getMannWhitney(Peak bpeak, Peak nbpeak, int shift) {
		double[] nbPositions = nbpeak.recreateCollection(shift);
		double[] bPositions = bpeak.recreateCollection(shift);
		return new MannWhitney(nbPositions,bPositions);
	}
	
	private void writePeaks() throws IOException {
		peakCaller.writeSingleSampleScanPeaksAllSamples();
	}
	
	public class Peak extends GeneWindow{
		private Annotation window;
		private double totalScore;
		private Map<Integer,Double> posCountMap;
		private double totalCounts;
		
		public Peak(Annotation region, AlignmentModel data) {
			super(region);
			totalScore = data.getCount(region, true);
			posCountMap = createPosCountMap(region, data);
			totalCounts = findTotalCounts(posCountMap);
		}
		
		public Peak(Annotation region, WigReader wr) throws IOException {
			super(region);
			//WigReader wr = new WigReader(wigFile);
			posCountMap = wr.getWigData(region);
			totalCounts = findTotalCounts(posCountMap);
		}
		
		private Map<Integer,Double> createPosCountMap(Annotation region, AlignmentModel data) {
			Map<Integer,Double> rtrn = new TreeMap<Integer,Double>();
			int pos;
			for (int i = 0; i < region.size(); i++) {
				int refPos = region.getReferenceCoordinateAtPosition(i, false);
				rtrn.put(refPos, data.getCount(new BasicAnnotation(region.getChr(),refPos,refPos+1)));
			}
			/*
			WindowScoreIterator<CountScore> winItr = data.scan(region,1,0);
			while (winItr.hasNext()) {
				CountScore posScore = winItr.next();
				pos = posScore.getAnnotation().getStart();
				rtrn.put(pos, posScore.getScore());
			}*/
			return rtrn;
		}
		
		private double findTotalCounts(Map<Integer,Double> countMap) {
			double rtrn = 0;
			for (int pos : countMap.keySet()) {
				rtrn = rtrn + countMap.get(pos);
			}
			return rtrn;
		}
		
		public double getTotalCount() {
			return totalCounts;
		}
		
		public Map<Integer,Double> getPosCountMap() {
			return posCountMap;
		}
		
		public Annotation getAnnotation() {
			return window;
		}
		
		public double[] recreateCollection(int shift) {
			double[] rtrn = new double[(int) getTotalCount()];
			int curr = 0;
			for (int refPos : posCountMap.keySet()) {
				int relPos = getPositionAtReferenceCoordinate(refPos,false) + shift;
				double count = posCountMap.get(refPos);
				for (int i = 0; i<Math.floor(count); i++) {
					rtrn[curr] = relPos;
					curr++;
				}
			}
			return rtrn;
		}
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws IOException, InterruptedException {
		// Fields to use
		CommandLineParser p = new CommandLineParser();
		String listFile;
		String bedFile;
		p.addStringArg("-l", "list of sample bam files to compare; try to put the nonbarcoded sample first", true);
		p.addStringArg("-g", "gene bed file used to find peaks", false, null);
		p.addStringArg("-bo", "barcode output file name", true);
		p.addStringArg("-nbo", "nobarcode output file name", true);
		p.addStringArg("-bw", "barcode wig file; can be any type of wig file you want", false, null);
		p.addStringArg("-nbw", "nobarcode wig file; can be any type of wig file",false, null);
		p.addStringArg("-nbg", "non-barcoded peaks in bed file format; use if you've already calculated peaks from MultiSampleBindingSiteCaller", false, null);
		p.addStringArg("-bg", "barcoded peaks in bed file format; use if you've already calculated peaks from MultiSampleBindingSiteCaller", false, null);
		p.addStringArg("-npf", "file of precomputed nobarcode peaks in bed format", false, null);
		p.addStringArg("-bpf", "file of precomputed barcode peaks in bed format", false, null);
		// possibly consider adding wig files as not required inputs, otherwise we will use AlignmentModel.scan to find counts over peaks
		p.parse(args);
		
		bedFile = p.getStringArg("-g");
		listFile = p.getStringArg("-l");
		String barcodeOutFile = p.getStringArg("-bo");
		String nobarcodeOutFile = p.getStringArg("-nbo");
		String nbarFile = p.getStringArg("-nbg");
		String barFile = p.getStringArg("-bg");
		String bwFile = p.getStringArg("-bw");
		String nbwFile = p.getStringArg("-nbw");
		String nbpeakFile = p.getStringArg("-npf");
		String bpeakFile = p.getStringArg("-bpf");
		PeakComparator pc=null;
		if (bedFile != null) {
			pc = new PeakComparator(listFile, bedFile, nbpeakFile, bpeakFile);
		} else if (nbarFile != null && barFile != null) {
			// use this space for constructing w/ bamfile, barcoded peak file, and nonbarcoded peak file
		} else {
			pc = null;
			throw new IllegalArgumentException("Need to couple list file with either (1) a generic gene bed file to find peaks, or (2) both nonbarcode and barcode peak bed files");
		}
		logger.info("created peak comparator");
		//pc.writePeaks();
		logger.info("wrote peaks");
		Map<String,Collection<Gene>> bPeaks = pc.getPeaks(pc.getBarcodeSample());
		logger.info("got barcode peaks");
		Map<String,Collection<Gene>> nbPeaks = pc.getPeaks(pc.getNobarcodeSample());
		logger.info("got nobarcode peaks");
		/*
		Map<String,Double> avgPeakWidths = pc.calculateAvgPeakSize(bPeaks, nbPeaks);
		FileWriter w = new FileWriter(new File("average_peak_widths"));
		for (String sample : avgPeakWidths.keySet()) {
			w.write(sample + "\t" + avgPeakWidths.get(sample).toString() + "\n");
		}
		w.close();
		*/
		Map<String,Map<Integer,Double>> aggregMap = pc.AggregatePeaksFixedNBwindow(bPeaks,nbPeaks,bwFile,nbwFile, -30, 60, true, 100);
		
		//Map<String,Map<Integer,Double>> aggregMap = pc.AggregatePeaks(bPeaks,nbPeaks,bwFile,nbwFile,true);
		
		//Map<Annotation,Map<Double,Double>> medianPairs = pc.getAllPeakMedianPairs(bPeaks, nbPeaks, -30, 60, nbwFile, bwFile, 10);
		
		//double totalZscore = pc.findTotalDiffBtPeaks(bPeaks, nbPeaks, -30, 60, nbwFile, bwFile, 10);
		
		Map<Integer,Double> nbMap = aggregMap.get("nobarcode");
		Map<Integer,Double> bMap = aggregMap.get("barcode");
		
		FileWriter w1 = new FileWriter(new File(nobarcodeOutFile));
		FileWriter w2 = new FileWriter(new File(barcodeOutFile));
		
		for (int pos : nbMap.keySet()) {
			w1.write(Integer.toString(pos) + "\t" + Double.toString(nbMap.get(pos)) + "\n");
		}
		
		/*
		for (Annotation peak : medianPairs.keySet()) {
			for (double pos : medianPairs.get(peak).keySet()) {
				w1.write(Double.toString(pos) + "\t" + Double.toString(medianPairs.get(peak).get(pos)) + "\n");
			}
		}
		*/
		//pc.writeScoredPeaks(barcodeOutFile,nbPeaks, -30, 60, bwFile, nbwFile, 10);
		
		for (int pos : bMap.keySet()) {
			w2.write(Integer.toString(pos) + "\t" + Double.toString(bMap.get(pos)) + "\n");
		}
		
		//logger.info("Total z score for all peaks: " + totalZscore);
		w1.close();
		w2.close();
		
		
	}

}
