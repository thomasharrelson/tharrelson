package secondReadBias;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;

import broad.core.parser.CommandLineParser;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;

public class PeakSequenceComposition {
	private static Map<String, Double> data;
	/*
	 * All I want this class to do is use WigSequenceComposition for Peaks
	 * Then 
	 */
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-p", "Peak bed file", true);
		p.addStringArg("-g", "Genome fasta", true);
		//p.addStringArg("-w", "Wig file", true);
		p.addStringArg("-o", "Output table file", true);
		//p.addIntArg("-s", "Shift from wig count position", false, 0);
		
		p.parse(args);
		
		String peakFile = p.getStringArg("-p");
		String genomeFile = p.getStringArg("-g");
		//String wigFile = p.getStringArg("-w");
		String outFile = p.getStringArg("-o");
		//int shift = p.getIntArg("-s");
		data = new TreeMap<String, Double>();
		Map<String,Collection<Gene>> chrGenes = BEDFileParser.loadDataByChr(new File(peakFile));
		Map<String,Sequence> genome = FastaSequenceIO.getChrSequencesFromFasta(genomeFile);
		// build 'data' structure
		for (String chr : chrGenes.keySet()) {
			for (Gene g : chrGenes.get(chr)) {
				Sequence chrGenome = genome.get(chr);
				
				for (Annotation block : g.getBlocks()) {
					String blockBases = chrGenome.getRegion(block.getStart(),block.getEnd()).getSequenceBases();
					for (int i=0;i<block.size();i++) {
						//int refPos = g.getReferenceCoordinateAtPosition(i);
						
						//String base = chrGenome.getRegion(refPos,refPos+1).getSequenceBases();
						String base = blockBases.substring(i, i+1);
						if (data.containsKey(base)) {
							Double baseCount = data.get(base);
							baseCount++;
							data.put(base,baseCount);
						} else {
							data.put(base,1.0);
						}
					}
				}
			}
		}
		
		FileWriter w = new FileWriter(outFile);
		for (String base : data.keySet()) {
			w.write(base + "\t" + data.get(base) + "\n");
		}
		w.close();
	}

}
