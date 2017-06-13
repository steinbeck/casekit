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
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;

/**
 * SimilarityRanker uses a SpectrumPredictor and parses an SDF file, returning a configurable number of compounds and 
 * their ranked spectrum similarity.
 * 
 * @author steinbeck
 *
 */
public class SimilarityRanker {
	
	public boolean verbose = true;
	DecimalFormat df;
	int resultListSize = 100;
	
	public boolean isVerbose() {
		return verbose;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	public int getResultListSize() {
		return resultListSize;
	}

	public void setResultListSize(int resultListSize) {
		this.resultListSize = resultListSize;
	}

	public SimilarityRanker() {
		// TODO Auto-generated constructor stub
		df = new DecimalFormat();
		df.setMaximumFractionDigits(2);
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
	

	public ArrayList<Result> rank(File inFile, ArrayList<Signal> spectrum, File outFile, File HOSEFile) throws Exception 
	{
		/*
		 * Iterate of SDF file given by input file, predict a spectrum and calculate a similarity with the 
		 * spectrum given in @spectrum. 
		 * Store the 10 most similar spectra in a list and write them to outFile in the end  
		 */
		
		HOSECodePredictor predictor = new HOSECodePredictor(HOSEFile);
		IAtomContainer ac = null;
		double similarity = 0.0;
		double bestSimilarity = 1000000000.0;
		ArrayList<Result> results = new ArrayList<Result>();
		ResultComparator comp = new ResultComparator();
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
					ac.setProperty(CDKConstants.TITLE, "Distance " + df.format(similarity));
					results.add(new Result(ac, similarity));
					results.sort(comp);
					//After sorting, we remove the worst entry and thereby trim the results list to resultListSize
					if (results.size() == resultListSize) results.remove(resultListSize - 1);
				}
			}
			else results.add(new Result(ac, similarity));
		}
		iterator.close();
		System.out.println("Best similarity = " + bestSimilarity);
		return results;
	}

	public double calculateSimilarity(IAtomContainer ac, ArrayList<Signal> spectrum)
	{	
		double similarity = 0.0;
		double lastDiff = 0.0;
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
		for (int f = 0; f < spectrum.size(); f++)
		{
			lastDiff = 10000000000.0;
			matchFound = false;
			for (int g = 0; g < spectrum.size(); g++)
			{
				if (shifts[f] > spectrum.get(g).getShift().doubleValue()) diff = shifts[f] - spectrum.get(g).getShift().doubleValue();
				else diff = spectrum.get(g).getShift().doubleValue() - shifts[f];
				df.format(diff);
				if (diff < lastDiff) 
				{
					lastDiff = diff;
					matchFound = true;
				}
			}
			if (matchFound) similarity += lastDiff;
		}
		return similarity/spectrum.size();
	}
	
	public void reportResults(ArrayList<Result> results) throws IOException, CDKException
	{
		String filename = null;
		DepictionGenerator dg = null; 
		for (int f = 0; f < results.size(); f++)
		{
			filename = "~/temp/" + String.format("%03d", f) + "-mol.png";
			dg = new DepictionGenerator().withSize(800, 800).withAtomColors().withAtomValues().withMolTitle().withFillToFit();
			dg.depict(results.get(f).getAc()).writeTo(filename);
		}
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
		sr.setResultListSize(Integer.parseInt(args[4]));
		ArrayList<Signal> spectrum = null;
		ArrayList<Result> results = null;		
		try {
			spectrum = sr.readSpectrum(new File(args[1]));
			results = sr.rank(new File(args[0]),  spectrum, new File(args[2]), new File(args[3]));
			sr.reportResults(results);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
