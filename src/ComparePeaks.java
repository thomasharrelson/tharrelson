import java.io.IOException;
import java.util.Collection;

import broad.core.util.PipelineUtils;
import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Annotation;


public class ComparePeaks {
	Collection<Annotation> peaks;
	BEDFileParser bparser;
	private String bedFile;
	private String listFile;
	private boolean peaksCalculated = false;
	
	public ComparePeaks(String listFile, String bedFile) {
		this.listFile = listFile;
		this.bedFile = bedFile;
		try {
			calculatePeaks(listFile,bedFile);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public Collection<Annotation> getPeaks() throws IOException, InterruptedException {
		if (peaksCalculated) {
			return peaks;
		} else {
			calculatePeaks(listFile, bedFile);
			return peaks;
		}
	}
	
	private void calculatePeaks(String listFile, String bed) throws IOException, InterruptedException {
		
		// farm out another job using Pam's code.
		String jobID = Long.valueOf(System.currentTimeMillis()).toString();
		PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, "java -jar /seq/lincRNA/Pam/Software/ScriptureV2/MultiSampleBindingSiteCaller.java -g " + bed + " -b " + listFile + " -wsp true ", jobID + ".bsub", "hour", 4);
		PipelineUtils.waitForJobs(jobID, Runtime.getRuntime());
		// need to find out what the output file is called... need to test this jar file before using this code.
		// use a bedfile parser to generate a list of 'genes'
	}
	
	public void comparePeakShifts(Collection<Annotation> peaks) {
		
	}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
