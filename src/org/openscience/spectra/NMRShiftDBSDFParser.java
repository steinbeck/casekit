package org.openscience.spectra;

import java.awt.image.BufferedImageFilter;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.Collection;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Map;
import java.util.StringTokenizer;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.tools.HOSECodeGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

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
	int hoseCodeCounter = 0;
	int carbonCounter = 0;
	
	public NMRShiftDBSDFParser(String[] args) throws Exception
	{
		String inFile = args[0];
		String outFile = args[1];
		File fout = new File(outFile);
		FileOutputStream fos = new FileOutputStream(fout);
		bw = new BufferedWriter(new OutputStreamWriter(fos));
		
		IAtomContainer ac = new AtomContainer();
		
		IteratingSDFReader iterator = new IteratingSDFReader(
				new FileReader(inFile),
				DefaultChemObjectBuilder.getInstance()
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
				generateHOSECodes(ac);
			}		
			if (hydrogenNMR != null) hydrogenNMRCount++;			
			moleculeCount ++;			
			if (generatePictures) generatePicture(ac);
			//System.out.print(".");
			//if (moleculeCount == 100) break;
		}
		iterator.close();
		report = "File contains " + moleculeCount + " molecules with " + carbonNMRCount + " carbon spectra and " + hydrogenNMRCount + " hydrogen spectra.\n";
		report += hoseCodeCounter + " HOSE codes generated for " + carbonCounter + "carbon atoms, and written to file.";
		System.out.println(report);
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
	
	public void generateHOSECodes(IAtomContainer ac) throws Exception
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
				for (int g = 0; g < 4; g++)
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
	
	public void generatePicture(IAtomContainer ac) throws IOException, CDKException
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
		
		moleculeTitle = "~/temp/structures/" + String.format("%03d", moleculeCount) + "-mol.png";
		DepictionGenerator dg = new DepictionGenerator().withSize(800, 800).withAtomColors().withAtomValues().withMolTitle().withFillToFit();
		dg.depict(ac).writeTo(moleculeTitle);
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
