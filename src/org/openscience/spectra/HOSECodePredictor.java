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

public class HOSECodePredictor {
/**
 * We expect a tab separated text file with HOSE codes and associated shifts as input
 * @throws FileNotFoundException 
 * 
 */
	Hashtable<String,ArrayList<Double>> hoseLookup;
	boolean verbose = false;
	public HOSECodePredictor(String hoseTSVfileName) throws IOException
	{
		String line = null;
		StringTokenizer strtok;
		String hose; 
		Double shift;
		ArrayList<Double> shifts; 
		int linecounter = 0;

		BufferedReader br = new BufferedReader(new FileReader(hoseTSVfileName));
		hoseLookup = new Hashtable<String,ArrayList<Double>>();
		if (verbose) System.out.println("Start reading HOSE codes from " + hoseTSVfileName);
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
	
	public void predict(IAtomContainer ac) throws Exception
	{
		IAtom atom;
		String hose = null;
		Double shift = null;
		HOSECodeGenerator hcg = new HOSECodeGenerator();
		DecimalFormat df = new DecimalFormat();
		df.setMaximumFractionDigits(2);
		String[] sphereCount = 
		{
				"'",
				"''",
				"'''",
				"''''"
		};
		fixExplicitHydrogens(ac);
		if (verbose) System.out.println("Entering prediction module");
		for (int f = 0; f < ac.getAtomCount(); f++)
		{
			atom = ac.getAtom(f);
			if (verbose) System.out.println("Atom no. " + f);
			if (atom.getAtomicNumber() == 6)
			{
				// We descend from 4-sphere HOSE codes to those with lower spheres
				for (int g = 4; g > 0; g--)
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
			HOSECodePredictor hp = new HOSECodePredictor(args[0]);
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
