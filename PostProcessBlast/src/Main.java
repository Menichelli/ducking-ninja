import global.Global;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.TimeUnit;

import com.google.common.collect.Sets;

import model.BlastHit;
import model.CouplePfamPutative;
import model.CouplePutativePutative;
import model.PfamFamily;
import model.PutativeDomain;
import model.ValidatedDomain;
import parser.BlastResultsParser;
import parser.PfamParser;
import parser.ProteomeParser;
import parser.RParser;
import tools.CoupleGenerator;
import tools.FDREstimator;
import tools.HitsGathering;
import tools.printer.FastaPrinter;
import tools.printer.ResultsPrinter;
import tools.printer.StatsPrinter;

/**
 * 
 */

/**
 * @author christophe
 *
 */
public class Main {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		long startTime = System.currentTimeMillis();

		if(args.length==1) {
			Global.PATHPROPERTIES = args[0];
		} else Global.PATHPROPERTIES = "ppblast.properties";
		System.out.println("Initialization...");
		Global.init();
		System.out.println("done.");

		try {
			Set<String> proteinsWithAtLeastOnePutativeDomain = new HashSet<String>();
			Set<String> proteinsWithNoCouplePFPT = new HashSet<String>();
			Set<String> proteinsWithNoCouplePTPT = new HashSet<String>();
			Set<String> proteinsWithNoCouple = new HashSet<String>();
			
			//Step 1: Parse les resultats Blast, creation de la collection de Hits par proteine
			if(Global.VERBOSE) System.out.println("Parsing Blast results...");
			Map<String,List<BlastHit>> hitsByProt = BlastResultsParser.getHitsByProt();
			if(Global.VERBOSE) {
				int nbhit = 0;
				for(String s : hitsByProt.keySet()) {
					nbhit += hitsByProt.get(s).size();
				}
				System.out.println("Found "+nbhit+" valid hits on "+hitsByProt.keySet().size()+" different proteins.");
			}

			//Step 2: Parse le proteome de P.falciparum, genere un objet Proteome
			if(Global.VERBOSE) System.out.println("Parsing the species proteome...");
			Global.PROTEOME_AIMED = ProteomeParser.getProteome();
			if(Global.VERBOSE) System.out.println("Found "+Global.PROTEOME_AIMED.getProteins().size()+" proteins in the proteome.");

			//Step 3: Agglomere les hits sous forme de domaines potentiels
			if(Global.VERBOSE) System.out.println("Gathering putative domains...");
			Map<String,Set<PutativeDomain>> putativeDomainsByProt = new HashMap<String, Set<PutativeDomain>>();
			Set<PutativeDomain> tmpset,tmp;
			for(String protName : hitsByProt.keySet()) {
				List<BlastHit> hitsOnThisProt = hitsByProt.get(protName);
				while(true) { //break si le gatherHits ne trouve pas de putative
					tmpset = HitsGathering.gatherHits(hitsOnThisProt,Global.PROTEOME_AIMED.getProteinByName(protName),putativeDomainsByProt.get(protName));
					if(!tmpset.isEmpty()) {
						tmp = putativeDomainsByProt.get(protName);
						if(tmp==null) tmp = new HashSet<PutativeDomain>();
						tmp.addAll(tmpset);
						putativeDomainsByProt.put(protName,tmp);
						for(PutativeDomain pdom : tmpset) { //sur chaque domaine trouve
							for(BlastHit h : pdom.getBlastHits()) { //sur chacun de ses hits
								hitsOnThisProt.remove(h); //pour ne pas reutiliser le meme hit deux fois
							}
						}
					} else break;
				}
			}
			if(Global.VERBOSE) {
				int nbhit = 0;
				for(String s : putativeDomainsByProt.keySet()) {
					nbhit += putativeDomainsByProt.get(s).size();
				}
				System.out.println("Found "+nbhit+" putative domains on "+putativeDomainsByProt.keySet().size()+" different proteins.");
			}
			proteinsWithAtLeastOnePutativeDomain = putativeDomainsByProt.keySet();

			//Step 4: Parse Pfam pour recuprer la liste des proteines couvertes par famille
			if(Global.VERBOSE) System.out.println("Parsing Pfam families...");
			Set<PfamFamily> pfamFamilies = PfamParser.getFamilies();
			if(Global.VERBOSE) System.out.println("Found "+pfamFamilies.size()+" different families.");
			
			//Step 5: Genere tous les couples Pfam-PutativeDomain et envoie les info au StatsPrinter
			if(Global.VERBOSE) System.out.println("Testing all couples Pfam-Putative...");
			proteinsWithNoCouplePFPT.addAll(proteinsWithAtLeastOnePutativeDomain);
			proteinsWithNoCouple.addAll(proteinsWithAtLeastOnePutativeDomain);
			int nbCouplesTotal = 0, nbCouplesRetained = 0;
			int nbProtTreated = 0;
			for(String protName : putativeDomainsByProt.keySet()) {
				Set<CouplePfamPutative> couples = CoupleGenerator.getCouplePfamPutative(putativeDomainsByProt.get(protName), pfamFamilies);
				nbCouplesTotal+=couples.size();
				for(CouplePfamPutative couple : couples) {
					proteinsWithNoCouplePFPT.remove(couple.getPutative().getQueryName()+"_"+couple.getPutative().getQuerySpecies());
					proteinsWithNoCouple.remove(couple.getPutative().getQueryName()+"_"+couple.getPutative().getQuerySpecies());
					String pfamFamilyName = couple.getPfam().getFamilyName();
					String putativeDomainIdentifier = couple.getPutative().getIdentifier();
					int nbProtPfam = couple.getPfam().getAllProteinNames().size();
					Set<String> protsCoveringThePutativeDomain = couple.getPutative().getProteinsCoveringResidue(couple.getPutative().getBestPosition());
					int nbProtPutativeDomain = protsCoveringThePutativeDomain.size();
					int nbProtIntersec = Sets.intersection(couple.getPfam().getAllProteinNames(), protsCoveringThePutativeDomain).size();
					
					if(nbProtIntersec >= Global.NB_SEQ_INTERSECT) {
						nbCouplesRetained++;
						StatsPrinter.getInstance(Global.STATS_PFPT_PATH).addEntry(pfamFamilyName, putativeDomainIdentifier, nbProtPfam, nbProtPutativeDomain, nbProtIntersec);
					}
				}
				nbProtTreated++;
				if(Global.VERBOSE && Global.DYNAMIC_DISPLAY) System.out.print("\r"+nbProtTreated+"/"+putativeDomainsByProt.keySet().size()+" proteins tested.");
			}
			if(Global.VERBOSE && Global.DYNAMIC_DISPLAY) System.out.println();
			StatsPrinter.getInstance(Global.STATS_PFPT_PATH).close();
			if(Global.VERBOSE) System.out.println("Generated "+nbCouplesTotal+" couples Pfam-Putative. "+nbCouplesRetained+" couples with at least "+Global.NB_SEQ_INTERSECT+" proteins in common.");

			//Step 6: Comme precedemment mais avec les couples Putative-Putative
			if(Global.VERBOSE) System.out.println("Testing all couples Putative-Putative...");
			proteinsWithNoCouplePTPT.addAll(proteinsWithAtLeastOnePutativeDomain);
			nbCouplesTotal = 0;
			nbCouplesRetained = 0;
			nbProtTreated = 0;
			for(String protName : putativeDomainsByProt.keySet()) {
				Set<CouplePutativePutative> couples = CoupleGenerator.getCouplePutativePutative(putativeDomainsByProt.get(protName));
				nbCouplesTotal+=couples.size();
				for(CouplePutativePutative couple : couples) {
					proteinsWithNoCouplePTPT.remove(couple.getPutative1().getQueryName()+"_"+couple.getPutative1().getQuerySpecies());
					proteinsWithNoCouple.remove(couple.getPutative1().getQueryName()+"_"+couple.getPutative1().getQuerySpecies());
					String putativeDomainIdentifier1 = couple.getPutative1().getIdentifier();
					String putativeDomainIdentifier2 = couple.getPutative2().getIdentifier();

					Set<String> protsCoveringThePutativeDomain1 = couple.getPutative1().getProteinsCoveringResidue(couple.getPutative1().getBestPosition());
					int nbProtPutativeDomain1 = protsCoveringThePutativeDomain1.size();
					Set<String> protsCoveringThePutativeDomain2 = couple.getPutative2().getProteinsCoveringResidue(couple.getPutative2().getBestPosition());
					int nbProtPutativeDomain2 = protsCoveringThePutativeDomain2.size();

					int nbProtIntersec = Sets.intersection(protsCoveringThePutativeDomain1, protsCoveringThePutativeDomain2).size();

					if(nbProtIntersec >= Global.NB_SEQ_INTERSECT) {
						nbCouplesRetained++;
						StatsPrinter.getInstance(Global.STATS_PTPT_PATH).addEntry(putativeDomainIdentifier1, putativeDomainIdentifier2, nbProtPutativeDomain1, nbProtPutativeDomain2, nbProtIntersec);
					}
				}
				nbProtTreated++;
				if(Global.VERBOSE && Global.DYNAMIC_DISPLAY) System.out.print("\r"+nbProtTreated+"/"+putativeDomainsByProt.keySet().size()+" proteins tested.");
			}
			if(Global.VERBOSE && Global.DYNAMIC_DISPLAY) System.out.println();
			StatsPrinter.getInstance(Global.STATS_PTPT_PATH).close();
			if(Global.VERBOSE) System.out.println("Generated "+nbCouplesTotal+" couples Putative-Putative. "+nbCouplesRetained+" couples with at least "+Global.NB_SEQ_INTERSECT+" proteins in common.");
			
			if(Global.VERBOSE) {
				System.out.println("Number of proteins with at least one putative domain: "+proteinsWithAtLeastOnePutativeDomain.size());
				System.out.println("Number of proteins with at least one putative domain and no couple Pfam-Putative: "+proteinsWithNoCouplePFPT.size());
				System.out.println("Number of proteins with at least one putative domain and no couple Putative-Putative: "+proteinsWithNoCouplePTPT.size());
				System.out.println("Number of proteins with only one putative domain and no Pfam : "+proteinsWithNoCouple.size());
			}
			

			//Step 7: Run the R script to compute coocs
			if(Global.VERBOSE) System.out.print("Starting R script to compute coocurrence scores...");
			//Compute Pfam-Putative
			Process child1 = Runtime.getRuntime().exec("Rscript "+Global.R_SCRIPT_PATH+" "+Global.STATS_PFPT_PATH+" "+Global.R_RESULTS_PFPT_PATH);
			int code1 = child1.waitFor();
			switch (code1) {
			case 0:
				break;
			case 1:
				System.out.println("R script failed for PFPT.");
			}
			//Compute Putative-Putative
			Process child2 = Runtime.getRuntime().exec("Rscript "+Global.R_SCRIPT_PATH+" "+Global.STATS_PTPT_PATH+" "+Global.R_RESULTS_PTPT_PATH);
			int code2 = child2.waitFor();
			switch (code2) {
			case 0:
				break;
			case 1:
				System.out.println("R script failed for PTPT.");
			}
			if(Global.VERBOSE) System.out.println("done.");

			//Step 8: Parse R results to identify validated putative domains
			if(Global.VERBOSE) System.out.println("Parsing R results and store validated domains...");
			//charge les resultats du RParser
			Set<ValidatedDomain> validatedDomainsPFPT = RParser.getValidatedDomains(Global.R_RESULTS_PFPT_PATH);
			Set<ValidatedDomain> validatedDomainsPTPT = RParser.getValidatedDomains(Global.R_RESULTS_PTPT_PATH);
			//init
			Set<String> vDomainsPFPT = new HashSet<String>(); //les id des domaines valides par PFPT
			Set<String> vDomainsPTPT = new HashSet<String>(); //les id des domaines valides par PTPT
			Set<String> vDomainsAll = new HashSet<String>(); //les id de tous les domaines valides
			
			//Garde uniquement la meilleur P-valeur si un domaine est valide plusieurs fois
			Map<String,ValidatedDomain> keepBestPvalue = new HashMap<String, ValidatedDomain>();
			ValidatedDomain tmpVD;
			for(ValidatedDomain vd : validatedDomainsPFPT) {
				tmpVD = keepBestPvalue.get(vd.getIdentifierValidatedDomain());
				if(tmpVD!=null) {
					if(vd.getScore()<tmpVD.getScore()) {
						keepBestPvalue.put(vd.getIdentifierValidatedDomain(), vd);
					}
				} else {
					keepBestPvalue.put(vd.getIdentifierValidatedDomain(), vd);
				}
			}
			validatedDomainsPFPT.clear();
			for(String s : keepBestPvalue.keySet()) {
				validatedDomainsPFPT.add(keepBestPvalue.get(s));
			}
			//--Garde
			
			//Charge les domaines PFPT
			for(ValidatedDomain v : validatedDomainsPFPT) {
				vDomainsPFPT.add(v.getIdentifierValidatedDomain());
			}
			if(Global.VERBOSE) System.out.println("Found "+vDomainsPFPT.size()+" validated domains PFPT.");
			//Charge les domaines PTPT
			for(ValidatedDomain v : validatedDomainsPTPT) {
				vDomainsPTPT.add(v.getIdentifierValidatedDomain());
			}
			if(Global.VERBOSE) System.out.println("Found "+vDomainsPTPT.size()+" validated domains PTPT.");
			//Charge tous les domaines
			vDomainsAll.addAll(vDomainsPFPT);
			vDomainsAll.addAll(vDomainsPTPT);
			if(Global.VERBOSE) System.out.println("-> "+vDomainsAll.size()+" unique validated domains.");
			
			//Step 9: print results of validated putative domains
			Map<PutativeDomain,Set<String>> coocProtsPerDomain = new HashMap<PutativeDomain, Set<String>>();

			if(Global.VERBOSE) System.out.println("Printing results for validated domains PFPT...");
			Map<PutativeDomain,Set<String>> pdomainsPFPT = new HashMap<PutativeDomain, Set<String>>();
			for(ValidatedDomain vDomain : validatedDomainsPFPT) {
				PutativeDomain domain = null;
				//recupere l'objet PutativeDomain correspondant a l'identifiant
				for(PutativeDomain pd : putativeDomainsByProt.get(vDomain.getIdentifierValidatedDomain().split("/")[0])) {
					if(pd.getIdentifier().equals(vDomain.getIdentifierValidatedDomain())) {
						domain = pd;
						break;
					}
				}

				Set<String> validatingDomains = pdomainsPFPT.get(domain); //charge la liste des validants de ce domaine
				if(validatingDomains==null) validatingDomains = new HashSet<String>(); //initialise la list si necessaire
				validatingDomains.add(vDomain.getIdentifierValidatingDomain()); //ajoute l'identifiant du validant
				pdomainsPFPT.put(domain, validatingDomains); //remet la liste dans la map
				ResultsPrinter.getInstance().addEntry(domain,vDomain.getIdentifierValidatingDomain(),vDomain.getScore(),true); //print le resultat
			}

			//reitere la liste des domaines qui ont ete valide par PFPT
			for(PutativeDomain dom : pdomainsPFPT.keySet()) {
				Set<String> allowedProteins = new HashSet<String>();
				for(String fname : pdomainsPFPT.get(dom)) { //itere sur la liste des domaines validants le putative
					for(PfamFamily fam : pfamFamilies) { //itere tous les Pfam
						if(fam.getFamilyName().equals(fname)) { //recherche la famille qui correspond au fname
							allowedProteins.addAll(fam.getAllProteinNames()); //ajoute la liste des proteines ou on retrouve le Pfam
							break;
						}
					}
				}
				coocProtsPerDomain.put(dom, allowedProteins); //stock les proteines de tous les validants du dom
			}

			if(Global.VERBOSE) System.out.println("Printing results for validated domains PTPT...");
			Map<PutativeDomain,Set<String>> pdomainsPTPT = new HashMap<PutativeDomain, Set<String>>();
			for(ValidatedDomain vDomain : validatedDomainsPTPT) {
				PutativeDomain domain = null;
				for(PutativeDomain pd : putativeDomainsByProt.get(vDomain.getIdentifierValidatedDomain().split("/")[0])) {
					if(pd.getIdentifier().equals(vDomain.getIdentifierValidatedDomain())) {
						domain = pd;
						break;
					}
				}

				Set<String> validatingDomains = pdomainsPTPT.get(domain); //charge la liste des validants de ce domaine
				if(validatingDomains==null) validatingDomains = new HashSet<String>(); //init si necessaire
				validatingDomains.add(vDomain.getIdentifierValidatingDomain()); //ajoute l'identifiant
				pdomainsPTPT.put(domain, validatingDomains); //stock dans la map PTPT
				ResultsPrinter.getInstance().addEntry(domain,vDomain.getIdentifierValidatingDomain(),vDomain.getScore(),false); //print..
			}

			//reitere la liste des domaines PTPT
			for(PutativeDomain dom : pdomainsPTPT.keySet()) { //pour chaque putative PTPT
				Set<String> allowedProteins = new HashSet<String>(); //liste des prots des validants
				for(String validantID : pdomainsPTPT.get(dom)) { //pour chaque validant
					for(PutativeDomain validatingDomain : putativeDomainsByProt.get(dom.getQueryName()+"_"+dom.getQuerySpecies())) {
						if(validatingDomain.getIdentifier().equals(validantID) && !validatingDomain.equals(dom)) {
							allowedProteins.addAll(validatingDomain.getProteinsCoveringResidue(validatingDomain.getBestPosition()));
							break;
						}
					}
				}
				Set<String> aProts = coocProtsPerDomain.get(dom); //recupere les resultats de PFPT
				if(aProts==null) aProts = new HashSet<String>(); //s'il n'y a pas de PFPT, faut initialiser
//				aProts.addAll(allowedProteins);
//				coocProtsPerDomain.put(dom, aProts);
			}

			//Step 10: Print Fasta
			Set<String> initSet = new HashSet<String>();
			if(Global.VERBOSE) System.out.print("Initializing the FastaPrinter...");
			for(PutativeDomain dom : coocProtsPerDomain.keySet()) {
				initSet.addAll(coocProtsPerDomain.get(dom));
			}
			FastaPrinter.getInstance().init(initSet);
			if(Global.VERBOSE) System.out.println("Ready.");

			if(Global.VERBOSE) System.out.println("Printing Fasta...");
			int domPrinted = 0;
			for(PutativeDomain dom : coocProtsPerDomain.keySet()) {
				FastaPrinter.getInstance().printFasta(dom,coocProtsPerDomain.get(dom));
				domPrinted++;
				if(Global.VERBOSE && Global.DYNAMIC_DISPLAY) System.out.print("\r"+domPrinted+"/"+coocProtsPerDomain.keySet().size());
			}
			if(Global.VERBOSE && Global.DYNAMIC_DISPLAY) System.out.println();
			if(Global.VERBOSE) System.out.println("Printing done.");
			FastaPrinter.close();

			ResultsPrinter.getInstance().close();


			//Step 11: Estimate FDR
			if(Global.COMPUTE_FDR) {
				if(Global.VERBOSE) System.out.println("Computing the false detection rate...");
				FDREstimator.estimateFDR(putativeDomainsByProt, pfamFamilies, vDomainsPFPT.size());
			}

		} catch (Exception e) {
			e.printStackTrace();
		}

		long stopTime = System.currentTimeMillis();
		long millis = stopTime-startTime;
		String s = String.format("%d min, %d sec",
				TimeUnit.MILLISECONDS.toMinutes(millis),
				TimeUnit.MILLISECONDS.toSeconds(millis) -
				TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(millis))
				);
		System.out.println("\nExecution Time: "+s);
	}

}
