
/* 
* This Open Source Software is provided to you under the MIT License
 * Refer to doc/mit.license or https://opensource.org/licenses/MIT for more information
 * 
 * Copyright (c) 2017, Christoph Steinbeck 
 */
package casekit;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.StringTokenizer;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

/** 
 * SimilarityRanker uses a SpectrumPredictor and parses an SDF file, returning a configurable number of compounds and 
 * their ranked spectrum similarity.
 *
 * This Open Source Software is provided to you under the MIT License
 * Refer to doc/mit.license or https://opensource.org/licenses/MIT for more information
 * 
 * Copyright (c) 2017, Christoph Steinbeck 
 * 
 * @author steinbeck
 *
 */
public class SimilarityRanker {
	
	public boolean verbose = true;
	DecimalFormat df;
	public int resultListSize = 100;
	public String inFile = null;
	public String outPath = null;
	public String spectrumFile = null;
	public String hoseTSVFile = null;
	ArrayList<Signal> spectrum = null;
	ArrayList<Result> results = null;
	
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

	public void readSpectrum() throws NumberFormatException, IOException
	{
		String line;
		StringTokenizer strtok;
		int linecounter = 0;
		Double shift = null;
		Integer mult = null;
		Signal signal; 
		String tempString;
		ArrayList<Signal> spectrum = new ArrayList<Signal>();
		BufferedReader br = new BufferedReader(new FileReader(spectrumFile));
		if (verbose) System.out.println("Start reading spectrum from  " + spectrumFile);
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
		if (verbose) System.out.println("Read " + linecounter + " signals from spectrum in file " +  spectrumFile);
		
		this.spectrum = spectrum;
	}
	

	public ArrayList<Result> rank() throws Exception 
	{
		/*
		 * Iterate of SDF file given by input file, predict a spectrum and calculate a similarity with the 
		 * spectrum given in @spectrum. 
		 * Store the 10 most similar spectra in a list and write them to outFile in the end  
		 */
		
		HOSECodePredictor predictor = new HOSECodePredictor(hoseTSVFile);
		IAtomContainer ac = null;
		double similarity = 0.0;
		double bestSimilarity = 1000000000.0;
		results = new ArrayList<Result>();
		ResultComparator comp = new ResultComparator();
		IteratingSDFReader iterator = new IteratingSDFReader(
				new FileReader(inFile),
				SilentChemObjectBuilder.getInstance()
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
		if (verbose) System.out.println("Calculation finished. Best similarity = " + bestSimilarity);
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
	
	public void reportResults() throws Exception
	{
		String filename = null;
		DepictionGenerator dg = null; 
		if (!outPath.endsWith(File.separator)) 
			outPath += File.separator;
		for (int f = 0; f < results.size(); f++)
		{
			filename = outPath + String.format("%03d", f) + "-mol.png";
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
	
	private void parseArgs(String[] args) throws ParseException
	{
		Options options = setupOptions(args);	
		CommandLineParser parser = new DefaultParser();
		try {
			CommandLine cmd = parser.parse( options, args);
			this.inFile = cmd.getOptionValue("infile");
			this.hoseTSVFile = cmd.getOptionValue("hosecodes");
			this.outPath = cmd.getOptionValue("outpath");
			this.spectrumFile = cmd.getOptionValue("spectrum");
			if (cmd.hasOption("numbers")) this.resultListSize = Integer.parseInt(cmd.getOptionValue("numbers"));
			if (cmd.hasOption("verbose")) this.verbose = true;
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			HelpFormatter formatter = new HelpFormatter();
			formatter.setOptionComparator(null);
			 String header = "Ranke structures based on given experimental spectrum and similarity to predicted spectrum.\n\n";
			 String footer = "\nPlease report issues at https://github.com/steinbeck/spectra";
			formatter.printHelp( "java -jar spectra.jar casekit.SimilarityRanker", header, options, footer, true );
			throw e;
		}
	}
	
	private Options setupOptions(String[] args)
	{
		Options options = new Options();

		Option infile = Option.builder("i")
			     .required(true)
			     .hasArg()
			     .longOpt("infile")
			     .desc("filename of with SDF/MOL file of structures to be ranked (required)")
			     .build();
		options.addOption(infile);
		Option spectrumfile = Option.builder("p")
			     .required(true)
			     .hasArg()
			     .longOpt("spectrum")
			     .desc("filename of CSV file with spectrum. Format of each line: <shift>;<number of attached H> (required)")
			     .build();
		options.addOption(spectrumfile);	
		Option outpath = Option.builder("o")
			     .required(true)
			     .hasArg()
			     .longOpt("outpath")
			     .desc("path to store pictures of ranked output structures (required)")
			     .build();
		options.addOption(outpath);
		Option hosefile = Option.builder("s")
			     .required(true)
			     .hasArg()
			     .longOpt("hosecodes")
			     .desc("filename of TSV file with HOSE codes (required)")
			     .build();
		options.addOption(hosefile);
		Option outputnumber = Option.builder("n")
			     .hasArg()
			     .longOpt("number")
			     .desc("number of structures in output file. Default is 10, if this option is ommitted")
			     .build();
		options.addOption(outputnumber);

		Option verbose = Option.builder("v")
			     .required(false)
			     .longOpt("verbose")
			     .desc("generate messages about progress of operation")
			     .build();
		options.addOption(verbose);	

		return options;
	}
	
	
	public static void main(String[] args) {
		SimilarityRanker sr = new SimilarityRanker();
		try {
			sr.parseArgs(args);
			sr.readSpectrum();
			sr.rank();
			sr.reportResults();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

}
