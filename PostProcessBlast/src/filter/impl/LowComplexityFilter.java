/**
 * 
 */
package filter.impl;

import java.util.HashSet;
import java.util.Set;

import model.Protein;
import model.PutativeDomain;
import filter.Ifilter;
import global.Global;

/**
 * @author christophe
 *
 */
public class LowComplexityFilter implements Ifilter {

	/* (non-Javadoc)
	 * @see filter.Ifilter#filter(model.PutativeDomain)
	 */
	@Override
	public PutativeDomain filter(PutativeDomain domain) {
		PutativeDomain ret = null;
		
		Protein p = Global.PROTEOME_AIMED.getProteinByName(domain.getQueryName()+domain.getQuerySpecies());
		
		String sequence = "";
		for(int index=domain.getDomainStart(); index<domain.getDomainEnd(); index++) {
			sequence += p.getSequence().charAt(index);
		}
		
		if(!isLowComplexity(sequence)) ret = domain;
		
		return ret;
	}
	
	private boolean isLowComplexity(String sequence) {
		boolean ret;
		Set<Character> chars = new HashSet<Character>();
		for(char c : sequence.toCharArray()) chars.add(c);
		ret = (Math.log(chars.size())/Math.log(2) < 2.2); //Algo SEG
		return ret;
	}

}
