package org.openscience.spectra;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.StringTokenizer;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.interfaces.IAtom;
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
	
	public boolean verbose = true;

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
		String tempString;
		ArrayList<Signal> spectrum = new ArrayList<Signal>();
		BufferedReader br = new BufferedReader(new FileReader(inFile));
		if (verbose) System.out.println("Start reading spectrum from  " + inFile);
		while((line = br.readLine()) != null)
		{
			if (!line.startsWith("#") && line.trim().length() > 0)
			{
				strtok = new StringTokenizer(line, ";");
				if (verbose) System.out.println(line);
				linecounter++;
				
				shift = Double.parseDouble(strtok.nextToken().trim());
				mult = Integer.parseInt(strtok.nextToken().trim());
				signal = new Signal(shift, mult);
				spectrum.add(signal);
			}
		}	
		br.close();
		if (verbose) System.out.println("Read " + linecounter + " signals from spectrum in file " +  inFile);
		
		return spectrum;
	}
	

	public void rank(File inFile, ArrayList<Signal> spectrum, File outFile, File HOSEFile) throws Exception 
	{
		/*
		 * Iterate of SDF file given by input file, predict a spectrum and calculate a similarity with the 
		 * spectrum given in @spectrum. 
		 * Store the 10 most similar spectra in a list and write them to outFile in the end  
		 */
		
		
		HOSECodePredictor predictor = new HOSECodePredictor(HOSEFile);
		IAtomContainer ac = null;
		double similarity = 0.0;
		IAtomContainer bestAc = null;
		double bestSimilarity = 1000000000.0;
		ArrayList<Result> results = new ArrayList<Result>();
		ResultComparator comp = new ResultComparator();
		DepictionGenerator dg = null; 
		IteratingSDFReader iterator = new IteratingSDFReader(
				new FileReader(inFile),
				DefaultChemObjectBuilder.getInstance()
				);
		
		while (iterator.hasNext()) 
		{
			ac = iterator.next();
			predictor.predict(ac);
			similarity = calculateSimilarity(ac, spectrum);
			if (results.size() > 0)
			{
				if (similarity < results.get(results.size()-1).getScore()) 
				{
					bestSimilarity = similarity;
					bestAc = ac;
					ac.setProperty(CDKConstants.TITLE, "Similarity " + similarity);
					results.add(new Result(ac, similarity));
					results.sort(comp);
					//System.out.println(results.size());
					if (results.size() == 100) results.remove(99);
				}
			}
			else results.add(new Result(ac, similarity));
		}
		iterator.close();
		System.out.println("Best similarity = " + bestSimilarity);
		String filename = null;
		for (int f = 0; f < results.size(); f++)
		{
			filename = "~/temp/" + String.format("%03d", f) + "-mol.png";
			dg = new DepictionGenerator().withSize(800, 800).withAtomColors().withAtomValues().withMolTitle().withFillToFit();
			dg.depict(results.get(f).getAc()).writeTo(filename);
		}
		
	}
	
	public double calculateSimilarity(IAtomContainer ac, ArrayList<Signal> spectrum)
	{	
		double similarity = 0.0;
		double lastDiff = 0.0;
	
		DecimalFormat df = new DecimalFormat();
		df.setMaximumFractionDigits(2);
		int counter = 0;
		String shift = null;
		boolean matchFound = false;
		double diff = 0.0;
		double shifts[] = new double[spectrum.size()];
		for (IAtom atom : ac.atoms())
		{
			if (atom.getAtomicNumber() == 6)
			{
				shift = atom.getProperty(CDKConstants.NMRSHIFT_CARBON);
				if (shift != null) shifts[counter] = Double.parseDouble(shift);
				else shifts[counter] = -1.0;
				counter ++;
			}
		}
		//System.out.println(spectrum.size() + " " + counter);
		for (int f = 0; f < spectrum.size(); f++)
		{
			lastDiff = 10000000000.0;
			matchFound = false;
			for (int g = f + 1; g < spectrum.size(); g++)
			{
				if (shifts[f] > spectrum.get(g).getShift().doubleValue()) diff = shifts[f] - spectrum.get(g).getShift().doubleValue();
				else diff = spectrum.get(g).getShift().doubleValue() - shifts[f];
				df.format(diff);
				if (diff < lastDiff) 
				{
					lastDiff = diff;
					matchFound = true;
				}
				//System.out.println(shifts[f] + "-" + spectrum.get(g).getShift().doubleValue() + " " + diff );
			}
			//System.out.println("lastDiff = " + lastDiff);
			if (matchFound) similarity += lastDiff;
		}
		
		return similarity/spectrum.size();
	}
	

	
	class ResultComparator implements Comparator<Result>
	{

		public int compare(Result o1, Result o2) {

			if (o1.getScore() < o2.getScore()) return -1; 
			return 1;
		}
		
	}
	
	
	public static void main(String[] args) {
		SimilarityRanker sr = new SimilarityRanker();
		ArrayList<Signal> spectrum = null;
		try {
			spectrum = sr.readSpectrum(new File(args[1]));
			sr.rank(new File(args[0]),  spectrum, new File(args[2]), new File(args[3]));			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
