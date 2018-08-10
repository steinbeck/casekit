/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package casekit.NMR;

import casekit.NMR.model.Assignment;
import casekit.NMR.model.Spectrum;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.xml.parsers.ParserConfigurationException;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.xml.sax.SAXException;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class test {
    
    public static void main(String[] args) throws ParserConfigurationException, SAXException, CloneNotSupportedException, FileNotFoundException, SQLException, ClassNotFoundException {
        
        final int maxSpheres = 1;
                
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
        Spectrum spec = null;
        

        IAtomContainer ac;
        try {
            // HJ555
            projectName = "HJ555";
            process = new Process(molFormula_HJ555);
            spec = Process.parse1DNMRviaXML(Peaks13C_HJ555, "C");
            spec.setSpecType(CDKConstants.NMRSPECTYPE_1D_ARBITRARY);
            final ArrayList<Integer> indicesInAtomContainer1D_13C = process.set1DNMR(spec);
            final Assignment assignments1D_13C = new Assignment(spec);
            assignments1D_13C.setAssignments(0, indicesInAtomContainer1D_13C);
            System.out.println("assignments spectrum 13C: " + Arrays.toString(assignments1D_13C.getAssignments(0)));
            System.out.println("equivalences: " + spec.getEquivalences());
            
            spec = Process.parseDEPTviaXML(PeaksDEPT90_HJ555);
            spec.setSpecType(CDKConstants.NMRSPECTYPE_1D_DEPT90);
            Spectrum spec135 = Process.parseDEPTviaXML(PeaksDEPT135_HJ555);
            spec135.setSpecType(CDKConstants.NMRSPECTYPE_1D_DEPT135);
            final HashMap<String, ArrayList<Integer>> matches1D_DEPT = process.setDEPT(spec, spec135, tolC);
            final Assignment assignments1D_DEPT90 = new Assignment(spec);
            final Assignment assignments1D_DEPT135 = new Assignment(spec135);
            assignments1D_DEPT90.setAssignments(0, matches1D_DEPT.get("DEPT90"));
            assignments1D_DEPT135.setAssignments(0, matches1D_DEPT.get("DEPT135"));
            System.out.println("assignments spectrum DEPT90: " + Arrays.toString(assignments1D_DEPT90.getAssignments(0)));
            System.out.println("assignments spectrum DEPT135: " + Arrays.toString(assignments1D_DEPT135.getAssignments(0)));

            spec = Process.parseHSQCviaXML(PeaksHSQC_HJ555, "C");
            spec.setSpecType(CDKConstants.NMRSPECTYPE_2D_HSQC);
            final HashMap<String, ArrayList<Integer>> matches2D_HSQC = process.setHSQC(spec, tolC);
            final Assignment assignments2D_HSQC = new Assignment(spec);
            assignments2D_HSQC.setAssignments(0, matches2D_HSQC.get(Utils.getAtomTypeFromSpectrum(spec, 0)));
            assignments2D_HSQC.setAssignments(1, matches2D_HSQC.get(Utils.getAtomTypeFromSpectrum(spec, 1)));
            System.out.println("assignments spectrum HSQC dim1: " + Arrays.toString(assignments2D_HSQC.getAssignments(0)));
            System.out.println("assignments spectrum HSQC dim2: " + Arrays.toString(assignments2D_HSQC.getAssignments(1)));
            
            spec = Process.parseHHCOSYviaXML(PeaksCOSY_HJ555);
            spec.setSpecType(CDKConstants.NMRSPECTYPE_2D_HHCOSY);
            final HashMap<String, ArrayList<Integer>> matches2D_HHCOSY = process.setHHCOSY(spec, tolH);
            final Assignment assignments2D_HHCOSY = new Assignment(spec);
            assignments2D_HHCOSY.setAssignments(0, matches2D_HHCOSY.get(Utils.getAtomTypeFromSpectrum(spec, 0)));
            assignments2D_HHCOSY.setAssignments(1, matches2D_HHCOSY.get(Utils.getAtomTypeFromSpectrum(spec, 1)));
            System.out.println("assignments spectrum HHCOSY dim1: " + Arrays.toString(assignments2D_HHCOSY.getAssignments(0)));
            System.out.println("assignments spectrum HHCOSY dim2: " + Arrays.toString(assignments2D_HHCOSY.getAssignments(1)));

            spec = Process.parseINADEQUATEviaXML(PeaksINADEQUATE_HJ555);
            spec.setSpecType(CDKConstants.NMRSPECTYPE_2D_INADEQUATE);            
            final HashMap<String, ArrayList<Integer>> matches2D_INADEQUATE = process.setINADEQUATE(spec, tolC);
            final Assignment assignments2D_INADEQUATE = new Assignment(spec);
            assignments2D_INADEQUATE.setAssignments(0, matches2D_INADEQUATE.get(Utils.getAtomTypeFromSpectrum(spec, 0)));
            assignments2D_INADEQUATE.setAssignments(1, matches2D_INADEQUATE.get(Utils.getAtomTypeFromSpectrum(spec, 1)));
            System.out.println("assignments spectrum INADEQUATE dim1: " + Arrays.toString(assignments2D_INADEQUATE.getAssignments(0)));
            System.out.println("assignments spectrum INADEQUATE dim2: " + Arrays.toString(assignments2D_INADEQUATE.getAssignments(1)));

            spec = Process.parseHMBCviaXML(PeaksHMBC_HJ555, "C");
            spec.setSpecType(CDKConstants.NMRSPECTYPE_2D_HMBC);
            final HashMap<String, ArrayList<Integer>> matches2D_HMBC = process.setHMBC(spec, tolH, tolC);
            final Assignment assignments2D_HMBC = new Assignment(spec);
            assignments2D_HMBC.setAssignments(0, matches2D_HMBC.get(Utils.getAtomTypeFromSpectrum(spec, 0)));
            assignments2D_HMBC.setAssignments(1, matches2D_HMBC.get(Utils.getAtomTypeFromSpectrum(spec, 1)));
            System.out.println("assignments spectrum HMBC dim1: " + Arrays.toString(assignments2D_HMBC.getAssignments(0)));
            System.out.println("assignments spectrum HMBC dim2: " + Arrays.toString(assignments2D_HMBC.getAssignments(1)));
            
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
                        ", INADEQUATE: " + ac.getAtom(i).getProperty(CDKConstants.NMRSPECTYPE_2D_INADEQUATE) + ", HMBC: " + ac.getAtom(i).getProperty(CDKConstants.NMRSPECTYPE_2D_HMBC));
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
