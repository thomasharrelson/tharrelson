import java.util.*;
import java.util.logging.Logger;
import java.io.FileWriter;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.IOException;

import broad.core.parser.CommandLineParser;

import net.sf.samtools.util.CloseableIterator;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;
import nextgen.core.feature.GenomeWindow;


public class TalenDesigner {
	
	/*
	 * Parameters:
	 * genes2cut - collection of gene objects to find cut sites for.
	 * filename - name of the file containing the gene names to be cut.
	 * 
	 */
	private Collection<Gene> genes2cut;
	private static String filename;
	private Collection<TreeSet<GenomeWindow>> broadRegions;
	private Collection<TreeSet<GenomeWindow>> cutRegions;	// ask Pam about which set would be best here.
	
	Logger logger = Logger.getLogger(TalenDesigner.class.getName());
	
	public TalenDesigner(String file, char fileType) throws IOException {
		/* 
		 * fileType will designate the type of input file and thus the type of constructor to be used.
		 */
		filename = file;
		BufferedReader r1 = new BufferedReader(new FileReader(filename));
		// construct genes2cut from file
		switch (fileType) {
		case 'n': nameConstruct(r1);break;
		case 'p': posConstruct(r1);break;
		case 'b': bedConstruct(r1);break;
		default: r1.close(); throw new IllegalArgumentException(fileType + " is not a valid file type");
		}
		// use genes2cut to create broadRegions
		setBroadRegions(genes2cut);
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
			int geneStart = g.getStart();
			Strand orientation = g.getStrand();
			// TODO: figure out how to use the enum Strand, so that I can get the proper window
			GenomeWindow gWin = new GenomeWindow(chr,start,stop,orientation);
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
			genes2cut.add(gene);
			line = r1.readLine();
		}
	}
	
	private void bedConstruct(BufferedReader r1) {
		throw new IllegalArgumentException("Class not implemented for construction with a bed file of genes");
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-f", "File name containing all genes", true);
		p.addStringArg("-t", "Type of the file, must be single char ('b'-bed file, 'p'-position file, 'n'-gene name file, 's'-sequence file)", true);
		String f = p.getStringArg("-f");
		char type = p.getStringArg("-t").charAt(0);
		
		TalenDesigner tal = new TalenDesigner(f,type);
		
	}

}