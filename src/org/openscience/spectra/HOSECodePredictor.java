/* 
* This Open Source Software is provided to you under the MIT License
 * Refer to doc/mit.license or https://opensource.org/licenses/MIT for more information
 * 
 * Copyright (c) 2017, Christoph Steinbeck 
 */

package org.openscience.spectra;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.StringTokenizer;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.tools.HOSECodeGenerator;


/**
 * Predicts NMRS spectra by lookup of HOSE codes
 * (Bremser, W., HOSE - A Novel Substructure Code, Analytica Chimica Acta, 1978, 103:355-365)
 * parsed from a tab-separated value file, produced, for example, by the <code>NMRShiftDBSDFParser</code> 
 * tool provided in this package. 
 * The current version is hard coded to predict carbon spectra only.
 * 
 * <p>
 * This Open Source Software is provided to you under the MIT License
 * Refer to doc/mit.license or https://opensource.org/licenses/MIT for more information
 * 
 * Copyright (c) 2017, Christoph Steinbeck 
 * 
 * @author Christoph Steinbeck
 * 
 */
public class HOSECodePredictor {

	Hashtable<String,ArrayList<Double>> hoseLookup;
	boolean verbose = false;
	
	/**
	 * Initializes a HOSECodePredictor by reading HOSECodes and assigned shift values from 
	 * from a tab-separated values (TSV) file, produced, for example, by the <code>NMRShiftDBSDFParser</code> 
	 * tool provided in this package. 
	 * Each line in the TSV file has the format
	 * <p>
	 * <code>String hosecode &#60;&#92;t&#62;double shiftvalue</code>
	 * <p>
	 * @param hoseTSVfile tab-separated values (TSV) file with HOSECodes and assigned shift values
	 * @throws IOException Exception if TSV file cannot be read
	 */
	public HOSECodePredictor(File hoseTSVfile) throws IOException
	{
		String line = null;
		StringTokenizer strtok;
		String hose; 
		Double shift;
		ArrayList<Double> shifts; 
		int linecounter = 0;

		BufferedReader br = new BufferedReader(new FileReader(hoseTSVfile));
		hoseLookup = new Hashtable<String,ArrayList<Double>>();
		if (verbose) System.out.println("Start reading HOSE codes from " + hoseTSVfile);
		while((line = br.readLine()) != null)
		{
			strtok = new StringTokenizer(line, "\t");
			linecounter++;
			hose = strtok.nextToken();
			shift = Double.parseDouble(strtok.nextToken());
			//System.out.println(hose + " ---- " + shift);
			if (hoseLookup.containsKey(hose))
			{
				shifts = hoseLookup.get(hose);
				shifts.add(new Double(shift));
				//System.out.println("HOSE code already in hashtable. Adding to existing list.");
			}
			else
			{
				//System.out.println("HOSE code not in hashtable. Adding to new list.");
				shifts = new ArrayList<Double>();
				shifts.add(new Double(shift));
				hoseLookup.put(hose, shifts);
			}
		}	
		br.close();
		if (verbose) System.out.println("Finished reading " + linecounter + " lines of HOSE codes.");
	}
	/**
	 * Predicts NMR chemical shifts based on a given HOSE code table read by the 
	 * Constructor of this class. 
	 * The predicted chemical shifts are assigned to each atom by
	 * as a property of type <code>CDKConstants.NMRSHIFT_CARBON</code>.
	 * They are also written as <code>CDKConstants.COMMENT</code> to aid the
	 * <code>DepictionGenerator</code> class to generate a picture with shift annotations 
	 * for all carbon atoms.
	 * 
	 * @param ac The IAtomContainer for which to predict the shift values
	 * @throws Exception	Thrown if something goes wrong. 
	 */
	public void predict(IAtomContainer ac) throws Exception
	{
		IAtom atom;
		String hose = null;
		Double shift = null;
		HOSECodeGenerator hcg = new HOSECodeGenerator();
		DecimalFormat df = new DecimalFormat();
		df.setMaximumFractionDigits(2);
		/**
		 * A visual appendix to the <code>CDKConstants.NMRSHIFT_COMMENT</code> annotation to show 
		 * how many HOSE code spheres where used to predict this shift value. 
		 */
		String[] sphereCount = 
		{
				"'",
				"''",
				"'''",
				"''''",
				"'''''",
				"'''''",
		};
		fixExplicitHydrogens(ac);
		if (verbose) System.out.println("Entering prediction module");
		for (int f = 0; f < ac.getAtomCount(); f++)
		{
			atom = ac.getAtom(f);
			if (verbose) System.out.println("Atom no. " + f);
			if (atom.getAtomicNumber() == 6)
			{
				// We descend from 6-sphere HOSE codes to those with lower spheres
				for (int g = 6; g > 0; g--)
				{
					hose = hcg.getHOSECode(ac, atom, g);
					//System.out.println("Look-up for HOSE code " + hose);
					try{
						shift = getShift(hose);
						if (shift != null)
						{
							if (verbose) System.out.println("Shift " + df.format(shift) + " found with " + g + "-sphere HOSE code.");
							ac.getAtom(f).setProperty(CDKConstants.NMRSHIFT_CARBON, df.format(shift));
							ac.getAtom(f).setProperty(CDKConstants.COMMENT, df.format(shift) + sphereCount[g - 1]);
							//If we found a HOSE code of a higher sphere, we take that one and skip the lower ones
							break;
						}
					}
					catch(Exception e)
					{
						e.printStackTrace();
						
					}
				}
			}
		}
	}

	public Double getShift(String hose)
	{
		double shiftvalue = 0.0;
		if (!hoseLookup.containsKey(hose)) return null;
		ArrayList<Double> list = hoseLookup.get(hose);
		for (int f = 0; f < list.size(); f++)
		{
			shiftvalue = shiftvalue + ((Double) list.get(f)).doubleValue();
		}
		shiftvalue = shiftvalue / list.size();
		if (verbose) System.out.println("Predicted HOSE code from " + list.size() + " values");
		
		return new Double(shiftvalue);
	}
	
	public void generatePicture(IAtomContainer ac, String path) throws IOException, CDKException
	{
		String moleculeTitle = "";
		/* Path separators differ in operating systems. Unix uses slash, windows backslash.
		   That's why Java offers this constant File.pathSeparator since it knows what OS it is running on */
		if (path.endsWith(File.pathSeparator)) 
			moleculeTitle = path + "mol.png";
		else
			moleculeTitle = path + File.pathSeparator + "mol.png";
		DepictionGenerator dg = new DepictionGenerator().withSize(800, 800).withAtomColors().withAtomValues().withMolTitle().withFillToFit();
		dg.depict(ac).writeTo(moleculeTitle);
	}
	

	
	/**
	 * This predictor cannot handle explicit hydrogens. Where therefore convert them to implicit first
	 */
	void fixExplicitHydrogens(IAtomContainer ac)
	{
		IAtom atomB;
		for (IAtom atomA : ac.atoms())
		{
			if (atomA.getAtomicNumber() == 1)
			{
				atomB = ac.getConnectedAtomsList(atomA).get(0);
				atomB.setImplicitHydrogenCount(atomB.getImplicitHydrogenCount() +1 );
				ac.removeAtom(atomA);
			}
		}
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		IAtomContainer ac = null;
		try {
			HOSECodePredictor hp = new HOSECodePredictor(new File(args[0]));
			if (args[3].equals("true")) hp.verbose = true;
			IteratingSDFReader iterator = new IteratingSDFReader(
					new FileReader(args[1]),
					DefaultChemObjectBuilder.getInstance()
					);
			
			while (iterator.hasNext()) 
			{
				ac = iterator.next();
				hp.predict(ac);
				hp.generatePicture(ac, args[2]);
			}
			iterator.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			
		}

	}

}
