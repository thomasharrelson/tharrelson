package secondReadBias;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import org.apache.commons.io.output.NullOutputStream;
import org.apache.log4j.Logger;

import net.sf.picard.sam.BuildBamIndex;
import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.BAMIndex;
import net.sf.samtools.BAMRecordCodec;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.AbstractPairedEndAlignment;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;
import nextgen.core.readFilters.GenomicSpanFilter;
import nextgen.core.readFilters.ProperPairFilter;

import broad.core.parser.CommandLineParser;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;

public class SplitBamFileBySequence {

	static Logger logger = Logger.getLogger(SplitBamFileBySequence.class.getName());
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub

		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "bam file to split up", true);
		p.addStringArg("-g", "gene bed file", true);
		p.addStringArg("-f", "genome fasta file", true);
		p.addStringArg("-o", "output file prefix", false, null);
		
		p.parse(args);
		
		String bamFile = p.getStringArg("-b");
		String bedFile = p.getStringArg("-g");
		String genomeFile = p.getStringArg("-f");
		String outFile = p.getStringArg("-o");
		if (outFile == null) {
			outFile = bamFile;
		}
		
		Map<String,Sequence> genome = FastaSequenceIO.getChrSequencesFromFasta(genomeFile);
		Map<String,Collection<Gene>> genesByChr = BEDFileParser.loadDataByChr(new File(bedFile));
		AlignmentModel data = new AlignmentModel(bamFile,new TranscriptomeSpace(genesByChr),true);
		data.addFilter(new ProperPairFilter());
		data.addFilter(new GenomicSpanFilter(200));
		
		BAMFileWriter aWriter = new BAMFileWriter(new File(outFile + "_A.bam"));
		BAMFileWriter cWriter = new BAMFileWriter(new File(outFile + "_C.bam"));
		BAMFileWriter gWriter = new BAMFileWriter(new File(outFile + "_G.bam"));
		BAMFileWriter tWriter = new BAMFileWriter(new File(outFile + "_T.bam"));
		SAMFileReader reader = new SAMFileReader(new File(bamFile));
		SAMFileHeader header = reader.getFileHeader();
		aWriter.setHeader(header);
		cWriter.setHeader(header);
		gWriter.setHeader(header);
		tWriter.setHeader(header);
		
		BAMRecordCodec aCodec = new BAMRecordCodec(header);
		aCodec.setOutputStream(new NullOutputStream());
		BAMRecordCodec cCodec = new BAMRecordCodec(header);
		cCodec.setOutputStream(new NullOutputStream());
		BAMRecordCodec gCodec = new BAMRecordCodec(header);
		gCodec.setOutputStream(new NullOutputStream());
		BAMRecordCodec tCodec = new BAMRecordCodec(header);
		tCodec.setOutputStream(new NullOutputStream());
		
		for (String chr : genesByChr.keySet()) {
			Sequence chrSeq = genome.get(chr);
			String chrBases = chrSeq.getSequenceBases();
			for (Gene gene : genesByChr.get(chr)) {
				CloseableIterator<Alignment> readItr = data.getOverlappingReads(gene, false);
				while (readItr.hasNext()) {
					AbstractPairedEndAlignment read = (AbstractPairedEndAlignment) readItr.next();
					Strand orientation;
					Alignment firstMate = read.getFirstMate();
					Alignment secondMate = read.getSecondMate();
					if (read.toSAMRecord().getSecondOfPairFlag()) {
						orientation = firstMate.getFragmentStrand();
					} else {
						orientation = secondMate.getFragmentStrand();
					}
					char transcribedBase;
					int pos;
					if (orientation == Strand.NEGATIVE) {
						pos = read.getEnd();
						if (pos > chrBases.length()) { continue;}
						logger.info("neg strand");
						char tmpBase = chrBases.charAt(pos);
						transcribedBase = AlignmentWidthBias.reverseBase(tmpBase);
					} else if (orientation == Strand.POSITIVE){ 
						pos = read.getStart()-1;
						if (pos < 0) { continue;}
						logger.info("pos strand");
						transcribedBase = chrBases.charAt(pos);
					} else {continue;}
					if (transcribedBase=='A' || transcribedBase=='a') {
						aWriter.addAlignment(firstMate.toSAMRecord());
						aWriter.addAlignment(secondMate.toSAMRecord());
					} else if (transcribedBase=='C' || transcribedBase=='c') {
						cWriter.addAlignment(firstMate.toSAMRecord());
						cWriter.addAlignment(secondMate.toSAMRecord());
					} else if (transcribedBase=='G' || transcribedBase=='g') {
						gWriter.addAlignment(firstMate.toSAMRecord());
						gWriter.addAlignment(secondMate.toSAMRecord());
					} else if (transcribedBase=='T' || transcribedBase=='t') {
						tWriter.addAlignment(firstMate.toSAMRecord());
						tWriter.addAlignment(secondMate.toSAMRecord());
					}
				}
				readItr.close();
			}
		}
		aWriter.close();
		File transcriptomeBamIdxFile = new File( outFile+"_A.bam" + ".bai");
		if(transcriptomeBamIdxFile.exists()) { transcriptomeBamIdxFile.delete();}
		SAMFileReader reader2 = new SAMFileReader(new File(outFile+"_A.bam"));
		BuildBamIndex.createIndex(reader2,transcriptomeBamIdxFile);
		reader2.close();
		
		cWriter.close();
		transcriptomeBamIdxFile = new File( outFile+"_C.bam" + ".bai");
		if(transcriptomeBamIdxFile.exists()) { transcriptomeBamIdxFile.delete();}
		reader2 = new SAMFileReader(new File(outFile+"_C.bam"));
		BuildBamIndex.createIndex(reader2,transcriptomeBamIdxFile);
		reader2.close();
		
		gWriter.close();
		transcriptomeBamIdxFile = new File( outFile+"_G.bam" + ".bai");
		if(transcriptomeBamIdxFile.exists()) { transcriptomeBamIdxFile.delete();}
		reader2 = new SAMFileReader(new File(outFile+"_G.bam"));
		BuildBamIndex.createIndex(reader2,transcriptomeBamIdxFile);
		reader2.close();
		
		tWriter.close();
		transcriptomeBamIdxFile = new File( outFile+"_T.bam" + ".bai");
		if(transcriptomeBamIdxFile.exists()) { transcriptomeBamIdxFile.delete();}
		reader2 = new SAMFileReader(new File(outFile+"_T.bam"));
		BuildBamIndex.createIndex(reader2,transcriptomeBamIdxFile);
		reader2.close();
	}

}
