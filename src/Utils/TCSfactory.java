package Utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Gene;

import broad.core.parser.CommandLineParser;
import broad.core.util.PipelineUtils;
import broad.pda.annotation.BEDFileParser;

public class TCSfactory {

	static Logger logger = Logger.getLogger(TCSfactory.class.getName());
	/**
	 * Concatenate all files with the same extension in current directory
	 * @param fileNames
	 * @param extension
	 * @param output
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public static void concatenateFiles(Collection<String> fileNames, String extension, String output) throws IOException, InterruptedException {
		/*
		Runtime run = Runtime.getRuntime();
		String concatCmmd = "";
		for (String file : fileNames) {
			concatCmmd += " " + file + extension;
		}
		Process p = run.exec("cat" + concatCmmd + " > " + output + extension);
		BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
		while (br.ready()) {
			logger.info(br.readLine());
		}
		p.waitFor();
		p.destroy();
		*/
		FileWriter w = new FileWriter(output + extension);
		for (String file : fileNames) {
			BufferedReader br = new BufferedReader(new FileReader(file + extension));
			while (br.ready()) {
				String line = br.readLine();
				w.write(line + "\n");
			}
			br.close();
		}
		w.close();
		/*
		// remove files to organize folder
		Process p2 = run.exec("rm" + concatCmmd);
		p2.waitFor();
		p2.destroy();
		*/
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws IOException, InterruptedException {
		// TODO Auto-generated method stub
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "bed file of lincs to find CRISPRs for", true);
		p.addStringArg("-o", "output file prefix name", true);
		p.addStringArg("-c", "location of TCS jar file to run", true);
		
		p.parse(args);
		
		String bedFile = p.getStringArg("-b");
		String outFile = p.getStringArg("-o");
		String jar = p.getStringArg("-c");
		
		String cmmd = "java -Xmx30g -Xms15g -Xmn10g -jar " + jar + " -t b";
		Collection<String> fileMap = new TreeSet<String>();
		Collection<String> deleteFiles = new TreeSet<String>();
		ArrayList<String> jobs = new ArrayList<String>();
		BufferedReader bedReader = new BufferedReader(new FileReader(bedFile));
		while (bedReader.ready()) {
			String bedLine = bedReader.readLine();
			Gene linc = new Gene(bedLine, false);
			String outName = outFile + "_" + linc.getName();
			FileWriter w = new FileWriter(outName + ".bed");
			fileMap.add(outName);
			w.write(bedLine + "\n");
			w.close();
			cmmd = cmmd + " -g " + outName + ".bed -o " + outName;
			String jobID = Long.valueOf(System.currentTimeMillis()).toString();
			jobs.add(jobID);
			PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, jobID + ".bsub", "week", 32);
			// TODO: need to delete temp bed files
			deleteFiles.add(outName + ".bed");
			cmmd = "java -Xmx30g -Xms15g -Xmn10g -jar " + jar + " -t b";
			
		}
		bedReader.close();
		PipelineUtils.waitForAllJobs(jobs, Runtime.getRuntime());
		// TODO: need to concatenate all matching files into 96 well formats
		// there are 47 different files outputted by my program and I need to concatenate all of them.
		// what is the best/most efficient way to do this
		ArrayList<String> cutLocation = new ArrayList<String>();
		for (String deleteFile : deleteFiles) {
			File delete = new File(deleteFile);
			if (!delete.delete()) {
				System.out.println(deleteFile + " could not be deleted!");
			}
		}
		cutLocation.add("IN");cutLocation.add("TSS");cutLocation.add("PE");
		ArrayList<String> cutTag = new ArrayList<String>();
		cutTag.add("std");cutTag.add("GRNA"); cutTag.add("U6"); cutTag.add("gBlock");
		
		String extension;
		for (Integer i = 1; i<=4; i++) {
			for (String locID : cutLocation) {
				for (String tag : cutTag) {
					if (!tag.equals("std") && !tag.equals("gBlock")) {
						extension = "_" + tag + "_" + locID + "_" + i.toString();
					} else if (tag.equals("std")) {
						extension = "_" + locID + "_" + i.toString();
					} else {
						extension = "_" + tag + "_" + i.toString();
					}
					// concatenate all files with that extension
					concatenateFiles(fileMap, extension, outFile);
				}
			}
		}
		
		for (Integer i = 1; i<=2; i++) {
			for (String locID : cutLocation) {
				extension = "_surveyorPrimers_" + locID + "_" + i.toString();
				
				concatenateFiles(fileMap, extension, outFile);
			}
		}
		extension = "_flaggedLoci";
		concatenateFiles(fileMap,extension, outFile);
	}

}
