/**
 * 
 */
package tools;

import global.Global;

/**
 * @author christophe
 *
 */
public class CoocScoringModule {
	
	private CoocScoringModule() {}
	
	public static void compute(String stats_path, String results_path) throws Exception {
		Process child = Runtime.getRuntime().exec("Rscript "+Global.R_SCRIPT_PATH+" "+stats_path+" "+results_path);
		int code = child.waitFor();
		switch (code) {
		case 0:
			break;
		case 1:
			System.out.println("R script failed for: "+stats_path+" "+results_path);
		}
	}

}
