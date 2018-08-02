/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package casekit.NMR;

import casekit.NMR.model.Spectrum;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.Connection;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.xml.parsers.ParserConfigurationException;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.xml.sax.SAXException;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class test {
    
    public static void main(String[] args) throws ParserConfigurationException, SAXException, CloneNotSupportedException, FileNotFoundException, SQLException, ClassNotFoundException {
        
        final int maxSpheres = 1;
//        final String[] args2 = new String[]{"-i", "/Users/mwenk/Downloads/nmrshiftdb2withsignals.sd", "-o", "/Users/mwenk/Downloads/hose" + maxSpheres + ".tsv", "-m", String.valueOf(maxSpheres), "-v"};
//        try {
//            final NMRShiftDBSDFParser parser = new NMRShiftDBSDFParser(args2);
//        } catch (Exception ex) {
//            Logger.getLogger(test.class.getName()).log(Level.SEVERE, null, ex);
//        }

                
        final String Peaks13C_HJ555 = "/Users/mwenk/work/research/Beemelmanns-HKI-HJ555/13C/50/pdata/1/peaklist.xml";
        final String Peaks13C_HJ777 = "/Users/mwenk/work/research/Beemelmanns-HKI-HJ777/HJ777_13C_NMR.csv";
        final String PeaksH1_HJ555 = "/Users/mwenk/work/research/Beemelmanns-HKI-HJ555/1H/1/pdata/1/peaklist.xml";
        final String Peaks1H_HJ777 = "/Users/mwenk/work/research/Beemelmanns-HKI-HJ777/HJ777_1H_NMR.csv";
        final String PeaksDEPT90_HJ555 = "/Users/mwenk/work/research/Beemelmanns-HKI-HJ555/HJ555_DEPT90_pseudo.xml";
        final String PeaksDEPT90_HJ777 = "/Users/mwenk/work/research/Beemelmanns-HKI-HJ777/HJ777_DEPT90_NMR_pseudo.csv";
        final String PeaksDEPT135_HJ555 = "/Users/mwenk/work/research/Beemelmanns-HKI-HJ555/DEPT135/5/pdata/1/peaklist.xml";
        final String PeaksDEPT135_HJ777 = "/Users/mwenk/work/research/Beemelmanns-HKI-HJ777/HJ777_DEPT135_NMR.csv";
        final String PeaksHSQC_HJ555 = "/Users/mwenk/work/research/Beemelmanns-HKI-HJ555/HSQC/3/pdata/1/peaklist.xml";
        final String PeaksINADEQUATE_HJ555 = "/Users/mwenk/work/research/Beemelmanns-HKI-HJ555/HJ555_INADEQUATE_pseudo.xml";
        final String PeaksHSQC_HJ777 = "/Users/mwenk/work/research/Beemelmanns-HKI-HJ777/HJ777_HSQC_NMR.csv";
        final String PeaksHMBC_HJ555 = "/Users/mwenk/work/research/Beemelmanns-HKI-HJ555/HMBC/4/pdata/1/peaklist.xml";
        final String PeaksHMBC_HJ777 = "/Users/mwenk/work/research/Beemelmanns-HKI-HJ777/HJ777_HMBC_NMR.csv";
        final String PeaksCOSY_HJ555 = "/Users/mwenk/work/research/Beemelmanns-HKI-HJ555/COSY/2/pdata/1/peaklist.xml";
        final String PeaksCOSY_HJ777 = "/Users/mwenk/work/research/Beemelmanns-HKI-HJ777/HJ777_COSY_NMR.csv";

        final String pathToNMRShiftDB = "/Users/mwenk/Downloads/nmrshiftdb2withsignals.sd";
        final String pathToNMRShiftDBTest = "/Users/mwenk/Downloads/test.sdf";
        final String pathToNMRShiftDBHOSE = "/Users/mwenk/Downloads/hose" + maxSpheres + ".tsv";
        
        
        
        final double tolC = 0.5;
        final double tolH = 0.2;
        final String molFormulaString_HJ555 = "C21H19NO8";
        final String molFormulaString_HJ777 = "C28H25NO11S";
        final IMolecularFormula molFormula_HJ555 = MolecularFormulaManipulator.getMolecularFormula(molFormulaString_HJ555, DefaultChemObjectBuilder.getInstance());
        final IMolecularFormula molFormula_HJ777 = MolecularFormulaManipulator.getMolecularFormula(molFormulaString_HJ777, DefaultChemObjectBuilder.getInstance());
        String projectName = "";
        casekit.NMR.Process process = null;
        Spectrum spec;
        

        IAtomContainer ac;
        try {
            // HJ555
            projectName = "HJ555";
            process = new Process(molFormula_HJ555);
            spec = Utils.parseXML(Peaks13C_HJ555, 1, new int[]{1}, new String[]{"C"});
            process.set1DNMRShifts(spec);
//            process.parse1DNMR(Peaks13C_HJ555, "C");

            spec = Utils.parseXML(PeaksDEPT90_HJ555, 1, new int[]{1}, new String[]{"C"});
            Spectrum spec135 = Utils.parseXML(PeaksDEPT135_HJ555, 1, new int[]{1}, new String[]{"C"});
            int assignedHAtoms = process.setDEPT(spec, spec135, tolC);
//            int assignedHAtoms = process.parseDEPT(PeaksDEPT90_HJ555, PeaksDEPT135_HJ555, tolC);          
            System.out.println("assigned protons: " + assignedHAtoms);

            spec = Utils.parseXML(PeaksHSQC_HJ555, 2, new int[]{2, 1}, new String[]{"H", "C"});
            process.setHSQC(spec, tolC);
//            process.parseHSQC(PeaksHSQC_HJ555, "C", tolH);  

            spec = Utils.parseXML(PeaksCOSY_HJ555, 2, new int[]{2, 1}, new String[]{"H", "H"});
            process.setHHCOSY(spec, tolH);
//            process.parseHHCOSY(PeaksCOSY_HJ555, tolH);

            spec = Utils.parseXML(PeaksINADEQUATE_HJ555, 2, new int[]{2, 1}, new String[]{"C", "C"});
            process.setINADEQUATE(spec, tolC);
//            process.parseINADEQUATE(PeaksINADEQUATE_HJ555, tolC);

            spec = Utils.parseXML(PeaksHMBC_HJ555, 2, new int[]{2, 1}, new String[]{"H", "C"});
            process.setHMBC(spec, tolH, tolC);
//            process.parseHMBC(PeaksHMBC_HJ555, "C", tolH, tolC);

            process.setEquivalentProperties();
            process.setBonds(new String[]{CDKConstants.NMRSPECTYPE_2D_HHCOSY, CDKConstants.NMRSPECTYPE_2D_INADEQUATE, CDKConstants.NMRSPECTYPE_2D_HMBC}); // without hybridizations
            process.createLSDFile(projectName, "/Users/mwenk/Downloads/testLSD", new String[]{"/Users/mwenk/work/software/LSD-3.4.9/Filters/", "/Users/mwenk/work/software/LSD-3.4.9/Filters/MOLGEN/badlist1/"});
                        
            
//            // definition of all possible bond combinations up to 6 valences
//            final String[] bondsSet = {"-", "--", "---", "----", "=", "==", "=-", "=--", "%", "%-"}; // up to 4 valences (carbon)
//            //"-----", "------", "=---", "=----", "==-", "==--", "===", "%%", "%--", "%---", "%=", "%=-"}; // up to 6 valences (e.g. sulfur)
//            final String[] neighborElems = new String[]{"C", "O", "N", "S", "P", "Br", "Cl"};
////            final IAtomContainerSet acSet = NMR.DB.getStructuresFromNMRShiftDBFile(pathToNMRShiftDB, 10); // ring size of 10 in aromaticity search (pubchem txt file)
//            final Connection DBConnection = NMR.DB.getDBConnection("jdbc:mysql://localhost/nmrshiftdb", "useUnicode=true&useJDBCCompliantTimezoneShift=true&useLegacyDatetimeCode=false&serverTimezone=UTC&useSSL=false", "root", "jmd2017a");
////            NMR.Utils.getSpectraIDsFromNMRShiftDB(DBConnection, 155.0, 156.0, "C");
////            final HashMap<String, ArrayList<Double>> lookup = NMR.Utils.getLookupTableFromNMRShiftDB(DBConnection, "C");
////            NMR.Utils.getRMS(lookup);
//            final int minShift = 0, maxShift = 220, stepSize = 10;
//            final String elem = "C";
//            NMR.DB.getRMS(DBConnection, minShift, maxShift, elem);
//            final int[][] neighborhoodCountsMatrix = NMR.DB.countNeighborhoodBonds(DBConnection, bondsSet, elem, neighborElems, minShift, maxShift, stepSize);
//            NMR.Utils.writeNeighborhoodBondsCountMatrix("/Users/mwenk/Downloads/countMatrix_" + elem + "_SQL.csv", neighborhoodCountsMatrix, bondsSet, elem, neighborElems, minShift, maxShift, stepSize);
            
            // create 1D spectrum
            // coffein: 27.8;0.0Q;9|29.6;0.0Q;10|33.5;0.0Q;11|107.8;0.0S;5|144.3;0.0D;7|147.5;0.0S;4|151.6;0.0S;2|155.3;0.0S;0|
//            final ArrayList<NMR.Signal> spectrum = new ArrayList<>();
//            spectrum.add(new Signal(elem, 27.8, "Q", null));
//            spectrum.add(new Signal(elem, 29.6, "Q", null));
//            spectrum.add(new Signal(elem, 33.5, "Q", null));
//            spectrum.add(new Signal(elem, 107.8, "S", null));
//            spectrum.add(new Signal(elem, 144.3, "D", null));
//            spectrum.add(new Signal(elem, 147.5, "S", null));
//            spectrum.add(new Signal(elem, 151.6, "S", null));
//            spectrum.add(new Signal(elem, 155.3, "S", null));
//            NMR.DB.matchSpectrumAgainstDB(DBConnection, spectrum, 0.1, null, stepSize);
//            proc.countNeighborhoodBonds(acSet, bondsSet, elem, neighborElems, minShift, maxShift, stepSize);
//            NMR.Utils.writeNeighborhoodBondsCountMatrix("/Users/mwenk/Downloads/countMatrix_" + elem + ".csv", proc.getNeighborhoodBondsCountMatrix(), bondsSet, elem, neighborElems, minShift, maxShift, stepSize);
            
            

        } catch (IOException ex) {
            Logger.getLogger(test.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        ac = process.getAtomContainer();
//        final HashMap<String, ArrayList<Integer>> atomTypeIndices = proc.getAtomTypeIndices();
        System.out.println("\n");
        System.out.println(process.getAtomTypeIndices());
        for (int i = 0; i< ac.getAtomCount(); i++) {
                System.out.println("i: " + i + " -> atom: " + ac.getAtom(i).getSymbol() + ", shift: " + ac.getAtom(i).getProperty(casekit.NMR.Utils.getNMRShiftConstant("C")) + ", #H: " + ac.getAtom(i).getImplicitHydrogenCount() + 
                        ", H shifts: " + ac.getAtom(i).getProperty(CDKConstants.NMRSPECTYPE_2D_HSQC) + ", Hybrid.: " + ac.getAtom(i).getHybridization() + ", HHCOSY: " + ac.getAtom(i).getProperty(CDKConstants.NMRSPECTYPE_2D_HHCOSY) + 
                        ", INADEQUATE: " + ac.getAtom(i).getProperty(CDKConstants.NMRSPECTYPE_2D_INADEQUATE) + ", HMBC: " + ac.getAtom(i).getProperty(CDKConstants.NMRSPECTYPE_2D_HMBC) + ", EQUAL: " + ac.getAtom(i).getProperty(casekit.NMR.ParseRawData.PROP_EQUIVALENCE));
        }
        System.out.println("\nbond count: " + ac.getBondCount() + ":");
        for (IBond bond : ac.bonds()) {
            System.out.println("bond: " + bond);
        }
        

//        System.out.println("\n\nOpen Bonds:\n");    
//        for (int i = 0; i < ac.getAtomCount(); i++) {
//            Utils.getOpenBonds(ac, i);
//        }
            

//        try {
//          Utils.convertSDFtoLSD("/Users/mwenk/work/software/molgen5.02/badlist2.sdf", "/Users/mwenk/Downloads/", "/Users/mwenk/work/software/LSD-3.4.9/Mol2abSrc");
//        } catch (CDKException ex) {
//            Logger.getLogger(test.class.getName()).log(Level.SEVERE, null, ex);
//        } catch (IOException ex) {
//            Logger.getLogger(test.class.getName()).log(Level.SEVERE, null, ex);
//        }
        
    }
}
