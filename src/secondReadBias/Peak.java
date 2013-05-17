package secondReadBias;

import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;

public class Peak extends BasicAnnotation{

	static Logger logger = Logger.getLogger(Peak.class.getName());
	private Map<String, Annotation> splitPeakMap;
	private static String A_PEAK_ID = "Apeak";
	private static String C_PEAK_ID = "Cpeak";
	private static String G_PEAK_ID = "Gpeak";
	private static String T_PEAK_ID = "Tpeak";
	
	public Peak(Annotation region, Map<String, Annotation> splitPeaks) {
		super(region);
		createSplitPeakMap(splitPeaks);
	}
	
	public Peak(String chr, int start, int end, Map<String, Annotation> splitPeaks) {
		super(chr, start, end);
		createSplitPeakMap(splitPeaks);
	}
	
	public Peak(String chr, int start, int end, Strand orientation, Map<String, Annotation> splitPeaks) {
		super(chr, start, end, orientation);
		createSplitPeakMap(splitPeaks);
	}

	private void createSplitPeakMap(Map<String, Annotation> splitPeaks) {
		splitPeakMap = new TreeMap<String, Annotation>();
		for (String peakID : splitPeaks.keySet()) {
			if (!peakID.equals(A_PEAK_ID) && !peakID.equals(C_PEAK_ID) && !peakID.equals(G_PEAK_ID) && !peakID.equals(T_PEAK_ID)) {
				throw new IllegalArgumentException("key from map of split up peaks did not match any of the predetermined peak ID's; must be one of Apeak, Cpeak, Gpeak, or Tpeak");
			} else {
				// valid key
				splitPeakMap.put(peakID, splitPeaks.get(peakID));
			}
		}
	}
	
	public Annotation getPeakA() {
		if (splitPeakMap.containsKey(A_PEAK_ID))
			return splitPeakMap.get(A_PEAK_ID);
		else {
			logger.warn("no 'A' peak for this peak; returning null");
			return null;
		}
	}
	
	public Annotation getPeakC() {
		
		if (splitPeakMap.containsKey(C_PEAK_ID))
			return splitPeakMap.get(C_PEAK_ID);
		else {
			logger.warn("no 'C' peak for this peak; returning null");
			return null;
		}
	}
	
	public Annotation getPeakG() {
		
		if (splitPeakMap.containsKey(G_PEAK_ID))
			return splitPeakMap.get(G_PEAK_ID);
		else {
			logger.warn("no 'G' peak for this peak; returning null");
			return null;
		}
	}
	
	public Annotation getPeakT() {
		
		if (splitPeakMap.containsKey(T_PEAK_ID))
			return splitPeakMap.get(T_PEAK_ID);
		else {
			logger.warn("no 'T' peak for this peak; returning null");
			return null;
		}
	}
	
	public Integer getCAsizeDiff() {
		if (getPeakA()==null || getPeakC()==null) {
			logger.warn("peak 'A' or peak 'C' are non existent; returning null");
			return null;
		}
		return (splitPeakMap.get(C_PEAK_ID).size() - splitPeakMap.get(A_PEAK_ID).size());
	}
	
	public Integer getCGsizeDiff() {
		if (getPeakG()==null || getPeakC()==null) {
			logger.warn("peak 'G' or peak 'C' are non existent; returning null");
			return null;
		}
		return (splitPeakMap.get(C_PEAK_ID).size() - splitPeakMap.get(G_PEAK_ID).size());
	}

	public Integer getCTsizeDiff() {
		if (getPeakT()==null || getPeakC()==null) {
			logger.warn("peak 'T' or peak 'C' are non existent; returning null");
			return null;
		}
		return (splitPeakMap.get(C_PEAK_ID).size() - splitPeakMap.get(T_PEAK_ID).size());
	}
	
	public Double getCAshift() {
		if (getPeakA()==null || getPeakC()==null) {
			logger.warn("peak 'A' or peak 'C' are non existent; returning 0.0 so that writing a file in PeakDifference is still possible");
			return null;
		}
		double middleC = getPeakC().getReferenceCoordinateAtPosition((int) Math.floor(getPeakC().getSize()/2));
		double middleA = getPeakA().getReferenceCoordinateAtPosition((int) Math.floor(getPeakA().getSize()/2));
		if (this.isNegativeStrand()) {
			return middleC - middleA;
		} else {
			return middleA - middleC;
		}
	}
	
	public Double getCGshift() {
		if (getPeakG()==null || getPeakC()==null) {
			logger.warn("peak 'G' or peak 'C' are non existent; returning 0.0 so that writing a file in PeakDifference is still possible");
			return null;
		}
		double middleC = getPeakC().getReferenceCoordinateAtPosition((int) Math.floor(getPeakC().getSize()/2));
		double middleG = getPeakG().getReferenceCoordinateAtPosition((int) Math.floor(getPeakG().getSize()/2));
		if (this.isNegativeStrand()) {
			return middleC - middleG;
		} else {
			return middleG - middleC;
		}
	}

	public Double getCTshift() {
		if (getPeakT()==null || getPeakC()==null) {
			logger.warn("peak 'T' or peak 'C' are non existent; returning 0.0 so that writing a file in PeakDifference is still possible");
			return null;
		}
		double middleC = getPeakC().getReferenceCoordinateAtPosition((int) Math.floor(getPeakC().getSize()/2));
		double middleT = getPeakT().getReferenceCoordinateAtPosition((int) Math.floor(getPeakT().getSize()/2));
		if (this.isNegativeStrand()) {
			return middleC - middleT;
		} else {
			return middleT - middleC;
		}
	}
}
