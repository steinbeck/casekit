package org.openscience.spectra;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.StringTokenizer;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;

/**
 * SimilarityRanker takes a SpectrumPredictor and parses an SDF file, returning a configurable number of compounds and 
 * their ranked spectrum similarity.
 * 
 * @author steinbeck
 *
 */
public class SimilarityRanker {
	
	public boolean verbose = false;

	public SimilarityRanker() {
		// TODO Auto-generated constructor stub
	}

	public ArrayList<Signal> readSpectrum(File inFile) throws NumberFormatException, IOException
	{
		String line;
		StringTokenizer strtok;
		int linecounter = 0;
		Double shift = null;
		Integer mult = null;
		Signal signal; 
		ArrayList<Signal> spectrum = new ArrayList<Signal>();
		BufferedReader br = new BufferedReader(new FileReader(inFile));
		if (verbose) System.out.println("Start reading spectrum from  " + inFile);
		while((line = br.readLine()) != null)
		{
			strtok = new StringTokenizer(line, ";");
			linecounter++;
			shift = Double.parseDouble(strtok.nextToken());
			mult = Integer.parseInt(strtok.nextToken());
			signal = new Signal(shift, mult);
			spectrum.add(signal);
		}	
		br.close();
		if (verbose) System.out.println("Read " + linecounter + " signals from spectrum in file " +  inFile);
		
		return spectrum;
	}
	

	public void rank(File inFile, ArrayList<Signal> spectrum, File outFile) 
	{
		

		
	}
	
	public static void main(String[] args) {
		SimilarityRanker sr = new SimilarityRanker();
		ArrayList<Signal> spectrum = null;
		try {
			spectrum = sr.readSpectrum(new File(args[0]));
			//sr.rank(new File(args[1]),  spectrum, new File(args[2]));
;			
		} catch (NumberFormatException | IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	

}
