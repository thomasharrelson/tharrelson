import java.util.Collection;
import java.util.Iterator;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;
import nextgen.core.feature.Window;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;


public class StaticMethods {
	static Logger logger = Logger.getLogger(StaticMethods.class.getName());
	
	public static Collection<Gene> collapseGenes(Collection<Gene> genes) {
		TreeSet<Gene> rtrn = new TreeSet<Gene>();
		TreeSet<Gene> tmpGenes = new TreeSet<Gene>();
		tmpGenes.addAll(genes);
		IntervalTree<Gene> tree = makeTree(genes);
		Collection<Gene> skip = new TreeSet<Gene>();
		String name;
		
		for (Gene rtrnGene : genes) {
			name = rtrnGene.getName();
			//if (skip.contains(gene)) continue;
			//Gene rtrnGene = new Gene(gene);
			if (tmpGenes.contains(rtrnGene)) {
				skip.add(rtrnGene);
				//Iterator<? extends Window> overlaps = transcriptomeSpace.getFragment(gene).iterator();
				Iterator<Node<Gene>> overItr = tree.overlappers(rtrnGene.getStart(),rtrnGene.getEnd());
				/*
				TreeSet<Annotation> overlapGenes = new TreeSet<Annotation>();
				while (overlaps.hasNext()) {
					overlapGenes.addAll(overlaps.next().getSourceAnnotations());
				}*/
				
				while (overItr.hasNext()) {
					Gene g2 = new Gene(overItr.next().getValue());
					if (rtrnGene.overlapsStranded(g2)) {
						rtrnGene = rtrnGene.takeUnion(g2);
						logger.info("taking union...");
						skip.add(g2);
					}
				}
				/*
				for (Gene checkOverlap : tmpGenes) {
					if (rtrnGene.overlaps(checkOverlap)) {
						rtrnGene = rtrnGene.takeUnion(checkOverlap);
						logger.info("Gene name: " + gene.getName());
						logger.info("collapsing w/ isoform: " + checkOverlap.getName());
						skip.add(checkOverlap);}
				}
				*/
				rtrnGene.setName(name);
				logger.info("finished collapsing gene into: " + rtrnGene);
				rtrn.add(rtrnGene);
				tmpGenes.removeAll(skip);
				Iterator<Gene> skipItr = skip.iterator();
				while (skipItr.hasNext()) {
					Gene rmv = skipItr.next();
					tree.remove(rmv.getStart(),rmv.getEnd(),rmv);
				}
				skip.clear();
			}
		}
		return rtrn;
	}
	
	
	public static IntervalTree<Gene> makeTree(Collection<Gene> genes) {
		IntervalTree<Gene> rtrn=new IntervalTree<Gene>();
		
		for(Gene gene: genes){
			int start=gene.getStart();
			int end=gene.getEnd();
			
			//if not already in then add
			Node<Gene> found=rtrn.find(start, end);
			if(found==null){
				rtrn.put(start, end, gene);
			}
			//else add isoform
			else{
				// not sure what to do here if using annotation and not genes
				Gene ref=found.getValue();
				
				rtrn.put(start, end, ref);
			}
			
		}
		return rtrn;
	}
	
	public static Collection<Window> checkOrientation(Collection<? extends Window> frags, Annotation correctStrand) {
		Collection<Window> rtrn = new TreeSet<Window>();
		Iterator<? extends Window> winItr = frags.iterator();
		while (winItr.hasNext()) {
			Window currWin = winItr.next();
			if (currWin.getStrand().equals(correctStrand.getStrand())) {
				rtrn.add(currWin);
			}
		}
		return rtrn;
	}
	/*
	public static IntervalTree<Gene> makeTree(Collection<? extends Annotation> genes) {
		IntervalTree<Gene> rtrn=new IntervalTree<Gene>();
		
		for(Gene gene: genes){
			int start=gene.getStart();
			int end=gene.getEnd();
			
			//if not already in then add
			Node<Gene> found=rtrn.find(start, end);
			if(found==null){
				rtrn.put(start, end, gene);
			}
			//else add isoform
			else{
				// not sure what to do here if using annotation and not genes
				Gene ref=found.getValue();
				ref.addIsoform(gene);
				rtrn.put(start, end, ref);
			}
			
		}
		return rtrn;
	}*/
}
