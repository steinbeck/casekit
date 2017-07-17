/* 
* This Open Source Software is provided to you under the MIT License
 * Refer to doc/mit.license or https://opensource.org/licenses/MIT for more information
 * 
 * Copyright (c) 2017, Christoph Steinbeck 
 */
package org.openscience.spectra;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.StringTokenizer;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.HOSECodeGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/** 
 * Helper class to parse an NMRShiftDB SDF file with spectra assignments 
 * and convert it to a tab-separated values file with HOSE codes
 * (Bremser, W., HOSE - A Novel Substructure Code, Analytica Chimica Acta, 1978, 103:355-365)
 * and associated shift values.
 * The TSV file can then be used by <code>HOSECodePredictor</code> to predict spectra 
 * 
 * @author Christoph Steinbeck
 */

public class NMRShiftDBSDFParser {
	BufferedWriter bw;
	IMolecularFormula formula = null;
	String comment = null;
	String carbonNMR = null;
	String hydrogenNMR = null;
	String moleculeTitle = null;
	int carbonNMRCount = 0; 
	int hydrogenNMRCount = 0;
	int moleculeCount = 0;
	String report = "";
	String temp = "";
	boolean generatePictures = false;
	String picdir = null;
	int hoseCodeCounter = 0;
	int carbonCounter = 0;
	String inFile = null;
	String outFile = null;
	boolean verbose = false;
	int maxSpheres;
	
	public NMRShiftDBSDFParser(String[] args) throws Exception
	{
		parseArgs(args);
		if (verbose) System.out.println("Starting HOSE code generation with " + maxSpheres + " sphere from " + inFile);
		File fout = new File(outFile);
		FileOutputStream fos = new FileOutputStream(fout);
		bw = new BufferedWriter(new OutputStreamWriter(fos));
		
		IAtomContainer ac = SilentChemObjectBuilder.getInstance().newAtomContainer();
		
		IteratingSDFReader iterator = new IteratingSDFReader(
				new FileReader(inFile),
				SilentChemObjectBuilder.getInstance()
				);
		
		while (iterator.hasNext()) 
		{
			ac = iterator.next();
			carbonNMR = (String)ac.getProperty("Spectrum 13C 0");
			hydrogenNMR = (String)ac.getProperty("Spectrum 1H 0");			
			if (carbonNMR != null) 
			{
				carbonNMRCount++;
				ac = assignCarbonNMR(ac, carbonNMR);
				generateHOSECodes(ac, maxSpheres);
			}		
			if (hydrogenNMR != null) hydrogenNMRCount++;			
			moleculeCount ++;			
			if (generatePictures) generatePicture(ac, picdir);
		}
		iterator.close();
		report = "File contains " + moleculeCount + " molecules with " + carbonNMRCount + " carbon spectra and " + hydrogenNMRCount + " hydrogen spectra.\n";
		report += hoseCodeCounter + " HOSE codes generated for " + carbonCounter + "carbon atoms, and written to file.";
		if (verbose) System.out.println(report);
		bw.close();
	}

	IAtomContainer assignCarbonNMR(IAtomContainer ac, String nmrString) throws IOException, CDKException
	{
		String sigString = null, shiftString = null, multString = null, atomNumString = null;
		StringTokenizer strTok1 = new StringTokenizer(nmrString, "|");
		StringTokenizer strTok2 = null;
		while (strTok1.hasMoreTokens())
		{
			sigString = strTok1.nextToken(); 			//System.out.println(sigString);
			strTok2 = new StringTokenizer(sigString, ";");
			while (strTok2.hasMoreTokens())
			{
				shiftString = strTok2.nextToken();		//System.out.println(shiftString);
				multString = strTok2.nextToken();		//System.out.println(multString);
				atomNumString = strTok2.nextToken();	//System.out.println(atomNumString);
				try
				{
					ac.getAtom(Integer.parseInt(atomNumString)).setProperty(CDKConstants.NMRSHIFT_CARBON, Double.parseDouble(shiftString));
					ac.getAtom(Integer.parseInt(atomNumString)).setProperty(CDKConstants.COMMENT, Double.parseDouble(shiftString));
				}catch(Exception exc)
				{
					System.out.println("Failed to assign shift to atom no. " + Integer.parseInt(atomNumString) + " in molecule no " + moleculeCount + ", title: " + ac.getProperty(CDKConstants.TITLE) + " with " + ac.getAtomCount() + " atoms. ");
				}
			}
		}
		return ac;
	}
	
	public void generateHOSECodes(IAtomContainer ac, int maxSpheres) throws Exception
	{
		String hose = null;
		HOSECodeGenerator hcg = new HOSECodeGenerator();
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(ac);
		Aromaticity.cdkLegacy().apply(ac);
		for (int f = 0; f < ac.getAtomCount(); f++)
		{
			if (ac.getAtom(f).getAtomicNumber() == 6)
			{
				carbonCounter ++;
				for (int g = 0; g < maxSpheres; g++)
				{
					hose = hcg.getHOSECode(ac, ac.getAtom(f), g + 1);
					if (hose != null && ac.getAtom(f).getProperty(CDKConstants.NMRSHIFT_CARBON) != null)
					{	
						bw.write(hose + "\t" + ac.getAtom(f).getProperty(CDKConstants.NMRSHIFT_CARBON));
						bw.newLine();
						hoseCodeCounter++;
					}
				}
			}
		}
	}
	
	public void generatePicture(IAtomContainer ac, String picdir) throws IOException, CDKException
	{
		try
		{
			temp = ac.getProperty(CDKConstants.TITLE);
			if (temp != null && temp.length() > 15) temp = temp.substring(0, 14) + "...";
			ac.setProperty(CDKConstants.TITLE, temp);
		}
		catch(Exception e)
		{
			System.out.println("Problem with title " + temp);
		}
		if (!picdir.endsWith(File.separator)) picdir += File.separator;
		moleculeTitle = picdir + String.format("%03d", moleculeCount) + "-mol.png";
		if (verbose) System.out.println(moleculeTitle);
		DepictionGenerator dg = new DepictionGenerator().withSize(800, 800).withAtomColors().withAtomValues().withMolTitle().withFillToFit();
		dg.depict(ac).writeTo(moleculeTitle);
	}
	
	private void parseArgs(String[] args) throws ParseException
	{
		Options options = setupOptions(args);	
		CommandLineParser parser = new DefaultParser();
		try {
			CommandLine cmd = parser.parse( options, args);
			this.inFile = cmd.getOptionValue("infile");
			this.outFile = cmd.getOptionValue("outfile");
			if (cmd.hasOption("maxspheres")) 
			{
				this.maxSpheres = Integer.parseInt(cmd.getOptionValue("maxspheres"));
			}
			if (cmd.hasOption("verbose")) this.verbose = true;

			if (cmd.hasOption("picdir")) 
			{
				this.generatePictures = true;
				this.picdir = cmd.getOptionValue("picdir");
				
			}
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			HelpFormatter formatter = new HelpFormatter();
			
			formatter.setOptionComparator(null);
			 String header = "Generates a table of HOSE codes and assigned shifts from an NMRShiftDB SDF file from http://nmrshiftdb.nmr.uni-koeln.de/portal/js_pane/P-Help.\n\n";
			 String footer = "\nPlease report issues at https://github.com/steinbeck/spectra";
			formatter.printHelp( "java -jar spectra.jar NMRShiftDBSDFParser", header, options, footer, true );
			throw new ParseException("Problem parsing command line");
		}
	}
	
	private Options setupOptions(String[] args)
	{
		Options options = new Options();
		Option infile = Option.builder("i")
			     .required(true)
			     .hasArg()
			     .longOpt("infile")
			     .desc("filename of NMRShiftDB SDF with spectra (required)")
			     .build();
		options.addOption(infile);
		Option outfile = Option.builder("o")
			     .required(true)
			     .hasArg()
			     .longOpt("outfile")
			     .desc("filename of generated HOSE code table (required)")
			     .build();
		options.addOption(outfile);
		Option verbose = Option.builder("v")
			     .required(false)
			     .longOpt("verbose")
			     .desc("generate messages about progress of operation")
			     .build();
		options.addOption(verbose);	
		Option picdir = Option.builder("d")
			     .required(false)
			     .hasArg()
			     .longOpt("picdir")
			     .desc("store pictures in given directory")
			     .build();
		options.addOption(picdir);
		Option maxspheres = Option.builder("m")
			     .required(false)
			     .hasArg()
			     .longOpt("maxspheres")
			     .desc("maximum sphere size up to which to generate HOSE codes")
			     .build();
		options.addOption(maxspheres);
		return options;
	}
	
	public static void main(String[] args)
	{
		try {
			new NMRShiftDBSDFParser(args);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
