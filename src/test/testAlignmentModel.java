package test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.CountScore;
import nextgen.core.model.score.WindowScoreIterator;

import broad.core.parser.CommandLineParser;
import broad.pda.annotation.BEDFileParser;

public class testAlignmentModel {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		CommandLineParser parser = new CommandLineParser();
		parser.addStringArg("-g", "gene bed file used to compare different AlignmentModel implementations", true);
		parser.addStringArg("-b", "bam file to be tested", true);
		parser.addStringArg("-o", "output file name; default is AlignmentModelTest.out", false, "AlignmentModelTest.out");
		parser.addStringArg("-t", "test gene bed file", true);
		
		parser.parse(args);
		
		String bedFile = parser.getStringArg("-g");
		String bamFile = parser.getStringArg("-b");
		String test = parser.getStringArg("-t");
		Map<String,Collection<Gene>> genesByChr = BEDFileParser.loadDataByChr(new File(bedFile));
		String out = parser.getStringArg("-o");
		FileWriter w = new FileWriter(new File(out));
		TranscriptomeSpace tspace = new TranscriptomeSpace(genesByChr);
		Map<String,Collection<Gene>> testGenes = BEDFileParser.loadDataByChr(new File(test));
		AlignmentModel data = new AlignmentModel(bamFile,tspace);
		
		w.write("Gene name\tWindowSize\tCount\n");
		for (String chr : testGenes.keySet()) {
			for (Gene g : testGenes.get(chr)) {
				WindowScoreIterator<CountScore> winItr = data.scan(g,1,0);
				while (winItr.hasNext()) {
					CountScore score = winItr.next();
					w.write(g.getName() + "\t1\t" + score.getAnnotation().getChr()+"\t"+score.getAnnotation().getStart()+"\t"+score.getAnnotation().getEnd()+"\t"+score.getCount() + "\n");
				}
				WindowScoreIterator<CountScore> winItr2 = data.scan(g, 100, 0);
				while (winItr2.hasNext()) {
					CountScore score = winItr2.next();
					
					w.write(g.getName() + "\t100\t" + score.getAnnotation().getChr()+"\t"+score.getAnnotation().getStart()+"\t"+score.getAnnotation().getEnd()+"\t"+score.getCount() + "\n");
				}
				
			}
		}
		
		w.close();
	}

}
