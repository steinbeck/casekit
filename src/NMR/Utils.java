/*
 * The MIT License
 *
 * Copyright 2018 Michael Wenk [https://github.com/michaelwenk].
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package NMR;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.w3c.dom.Document;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.lang3.StringUtils;
import org.openscience.cdk.Atom;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.KlekotaRothFingerprinter;
import org.openscience.cdk.fingerprint.SubstructureFingerprinter;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.qsar.descriptors.atomic.AtomValenceDescriptor;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.smarts.parser.SMARTSParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Utils {
    
       
    public Utils() throws Exception {

    }
    
    /**
     * Splits an SDF into single molecular files and converts each of them into the LSD substructure format. 
     * Therefore, the mol2ab executable provided by LSD is required.
     * 
     * @param pathSDF path to SDF to split
     * @param pathOut path to directory which should contain the splitted and converted structure files
     * @param pathMol2ab path to mol2ab executable provided by LSD
     * @throws FileNotFoundException
     * @throws CDKException
     * @throws IOException
     */
    public static void convertSDFtoLSD(final String pathSDF, final String pathOut, final String pathMol2ab) throws FileNotFoundException, CDKException, IOException{
        
        
        System.out.println("Conversion from SDF format to LSD format... ");
        

        IAtomContainer ac;
        
        IteratingSDFReader iterator = new IteratingSDFReader(
                new FileReader(pathSDF),
                SilentChemObjectBuilder.getInstance()
        );

        
        File fout;
        FileOutputStream fos;
        BufferedWriter bw;
        File foutPilot = new File(pathOut + "/pilot");
        FileOutputStream fosPilot = new FileOutputStream(foutPilot);
        BufferedWriter bwPilot = new BufferedWriter(new OutputStreamWriter(fosPilot));
        
        
        int i = 0;
        while (iterator.hasNext()) {
            i++;
            ac = iterator.next();
            String molID = String.valueOf(i);//(String) ac.getProperties().get("cdk:Remark");
//            molID = molID.replace(" ", "_");
            fout = new File(pathOut + "/" + molID + ".sdf");
            fos = new FileOutputStream(fout);
            bw = new BufferedWriter(new OutputStreamWriter(fos));
            
            SDFWriter wtr = new SDFWriter(bw);
//            Properties sdfWriterProps = new Properties();
//            sdfWriterProps.put("WriteAromaticBondTypes", "true");
//            wtr.addChemObjectIOListener(new PropertiesListener(sdfWriterProps));
//            wtr.customizeJob();
            
            wtr.write(ac);
            wtr.close();
            bw.close();
            
            bwPilot.write(molID + " " + fout.getPath());
            bwPilot.newLine();
            
        }
       
        iterator.close();
        bwPilot.close();
        System.out.println("Input file contained " + i + " molecules!\nSingle files created!");
         
        
        // should be replaced by e.g. the process command because:
        // 1. for very long files the program ends long before the conversion process (command) ends
        // 2. no control or output here
        Runtime.getRuntime().exec(pathMol2ab + "/mol2ab " + pathOut + " " + foutPilot.getPath());
        
        
        System.out.println("Conversion from SDF format to LSD format... DONE!");
        
    }
    
 
    /**
     * Creates an IAtomContainer object containing atoms without any bond
     * information, given by a molecular formula.
     *
     * @param molFormula Molecular Formula
     * @return
     * @deprecated
     */
    @Deprecated
    public static IAtomContainer createAtomContainer(final String molFormula) {

        HashMap<String, Integer> hash = NMR.Utils.getAtomCountsInMolecularFormula(molFormula);
        IAtomContainer ac = SilentChemObjectBuilder.getInstance().newAtomContainer();

        for (String elem : hash.keySet()) {
            // add atoms of current element
            ac = NMR.Utils.addAtoms(ac, elem, hash.get(elem));
        }

        return ac;
    }

    /**
     * Creates a HashMap with the number of atoms for each occurring atom type.
     *
     * @deprecated
     * @param molFormula
     * @return
     */
    @Deprecated
    public static HashMap<String, Integer> getAtomCountsInMolecularFormula(final String molFormula) {

        HashMap<String, Integer> hash = new HashMap<>();
        String[] molFormSplit = molFormula.split("[A-Z]");
        Matcher m = Pattern.compile("[A-Z]").matcher(molFormula);
        String elem;
        int noAtoms;
        int k = 1;

        while (m.find()) {
            // name of current element
            elem = molFormula.substring(m.start(), m.end());
            if (k >= molFormSplit.length || molFormSplit[k].isEmpty()) {
                // if no atom number is given then assume only one atom
                noAtoms = 1;
            } else if (Character.isLowerCase(molFormSplit[k].charAt(0))) {
                // if element's name contains two letters then extend it
                elem += molFormSplit[k].charAt(0);
                // if more than one atoms of that element with two letters exist
                if (molFormSplit[k].length() > 1) {
                    // check given atom number
                    noAtoms = Integer.parseInt(molFormSplit[k].substring(1));
                } else {
                    noAtoms = 1;
                }
            } else {
                // if atom number is given
                noAtoms = Integer.parseInt(molFormSplit[k].substring(0));
            }
            try {
                // add atom type and frequency to class hashmap
                hash.put(elem, noAtoms);
            } catch (Exception e) {
                System.err.println("Illegal element \"" + elem + "\" will be ignored!!!");
            }

            k++;
        }

        return hash;
    }

    /**
     *
     * @param ac
     * @param atomType
     * @return
     * @deprecated
     */
    @Deprecated
    public static int getAtomTypeCount(final IAtomContainer ac, final String atomType) {

        int noAtoms = 0;
        for (int i = 0; i < ac.getAtomCount(); i++) {
            if (ac.getAtom(i).getSymbol().equals(atomType)) {
                noAtoms++;
            }
        }

        return noAtoms;
    }

    /**
     * Creates atoms of the same atom type and store it into an atom container.
     *
     * @param ac Atom container
     * @param noAtoms Number of atoms to create
     * @param atomType Atom type (element's name, e.g. C or Br)
     * @return
     */
    public static IAtomContainer addAtoms(final IAtomContainer ac, final String atomType, final int noAtoms) throws IllegalArgumentException {

        for (int i = 0; i < noAtoms; i++) {
            ac.addAtom(new Atom(atomType));
        }

        return ac;
    }
    
    
    /**
     * Removes atoms from a given atom type from an atom container.
     *
     * @param ac IAtomContainer object where to remove the atoms
     * @param atomType Atom type (element's name, e.g. C or Br)
     * @return IAtomContainer where the atoms were removed
     */
    public static IAtomContainer removeAtoms(final IAtomContainer ac, final String atomType){
        
        List<IAtom> toRemoveList = new ArrayList<>();
        for (IAtom atomA : ac.atoms()) {
            if (atomA.getSymbol().equals(atomType)){// detect wether the current atom A is a from the given atom type 
                toRemoveList.add(atomA);
            }
        }

        for (IAtom iAtom : toRemoveList) {
            ac.removeAtom(iAtom);
        }
        
        return ac;
    }
    
    
    /**
     * Returns a list of atom indices in an atom container for a given atom 
     * type (e.g. C or Br)
     *
     * @param ac IAtomContainer to use for search
     * @param atomType Atom type to find in atom container 
     * @return
     */
    public static ArrayList<Integer> getAtomTypeIndicesInAtomContainer(final IAtomContainer ac, final String atomType){
        
        final ArrayList<Integer> indices = new ArrayList<>();
        for (int i = 0; i < ac.getAtomCount(); i++) {
            if(ac.getAtom(i).getSymbol().equals(atomType)){
                indices.add(i);
            }
        }
        
        return indices;
    }
    

    /**
     * Reads a specific column of a NMR peak table and stores it into an
     * ArrayList object.
     *
     * @param pathToPeakList path to NMR peak table
     * @param column column to select in peak table
     * @return ArrayList of Double shift values
     * @throws IOException
     */
    public static ArrayList<Double> parsePeakTable(final String pathToPeakList, final int column) throws IOException {

        final ArrayList<Double> shifts = new ArrayList<>();
        String line;
        String[] tokens;
        BufferedReader fileReader = new BufferedReader(new FileReader(pathToPeakList));
        while ((line = fileReader.readLine()) != null) {
            tokens = line.split(",");
            // get shift value 
            if (tokens[column].trim().matches("^[+|-]{0,1}\\d+\\.{0,1}\\d*")) {
                shifts.add(Double.parseDouble(tokens[column].trim()));
            }
        }
        fileReader.close();

        return shifts;
    }
    
    
    /**
     * Reads a NMR peak XML file and stores it into an
     * ArrayList object.
     * The XML file must be in Bruker's TopSpin format.
     *
     * @param pathToXML Path to XML file
     * @param ndim number of dimensions of given data 1 (1D) or 2 (2D)
     * @param attribute which attribute index in XML peak nodes should be used: 
     * 1 (shift of 1st dimension), 2 (shift of 2nd dimension if 2D data, 
     * intensity if 1D data) or 3 (intensity if 2D data)
     * 
     * @return ArrayList of Double shift values
     * @throws IOException
     * @throws javax.xml.parsers.ParserConfigurationException
     * @throws org.xml.sax.SAXException
     */
    public static ArrayList<Double> parseXML(final String pathToXML, final int ndim, final int attribute) throws IOException, ParserConfigurationException, SAXException {

        final ArrayList<Double> shifts = new ArrayList<>();
        final DocumentBuilderFactory docBuilderFactory = DocumentBuilderFactory.newInstance();
        final DocumentBuilder docBuilder = docBuilderFactory.newDocumentBuilder();
        final Document doc = docBuilder.parse(new File(pathToXML));

        // normalize text representation
        doc.getDocumentElement().normalize();
        final NodeList peakLists = doc.getElementsByTagName("Peak" + ndim + "D");
        for (int i = 0; i < peakLists.getLength(); i++) {
            shifts.add(Double.parseDouble(peakLists.item(i).getAttributes().item(attribute - 1).getNodeValue()));
        }

        return shifts;
    }
    
    
    /**
     * Corrects a match list regarding a given shift list and an atom container.
     * This is useful when two ore more shift values (e.g. DEPT shifts) match
     * with the same atom in the atom container. So the purpose here is to
     * enable more unambiguous matches. This method first looks for unambiguous
     * matches and calculates the median of the difference values between the
     * shift list values and the shifts of atom container. Then, all shift list
     * values are adjusted (+/-) with this median value.
     *
     * @param ac IAtomContainer to search
     * @param shifts Shift value list to match
     * @param matches Match list to correct
     * @param tol Tolerance value
     * @param atomType Element name (e.g. "C") which also occurrs in 
     * {@link testkit.Utils#getNMRShiftConstant(java.lang.String)}
     * @return
     */
    public static ArrayList<Integer> correctShiftMatches(final IAtomContainer ac, final ArrayList<Double> shifts, final ArrayList<Integer> matches, final double tol, final String atomType) {

        int matchIndex;
        // get differences of unique matches between query shift and ac shifts
        ArrayList<Double> diffs = new ArrayList<>();
        final HashSet<Integer> uniqueMatchIndicesSet = new HashSet<>(matches);
        for (Integer matchIndexAtomContainer : uniqueMatchIndicesSet) {
            if (Collections.frequency(matches, matchIndexAtomContainer) == 1) {
                matchIndex = matches.indexOf(matchIndexAtomContainer);
                if (matches.get(matchIndex) >= 0) {
                    diffs.add(shifts.get(matchIndex) - Double.parseDouble(ac.getAtom(matches.get(matchIndex)).getProperty(NMR.Utils.getNMRShiftConstant(atomType)).toString()));
                }
            }
        }
        // calculate the median of found unique match differences
        if (diffs.size() > 0) {
            final double median = NMR.Utils.getMedian(diffs);
            // add or subtract the median of the differences to all shift list values (input) and match again then
            for (int i = 0; i < shifts.size(); i++) {
                shifts.set(i, shifts.get(i) - median);
            }
            // rematch
            return NMR.Utils.findShiftMatches(ac, shifts, tol, atomType);
        } 
        
        return matches;
    }
    

    /**
     * Finds the matches with the lowest deviations between a given shift value
     * set and the atoms of an atom container. A tolerance value and NMRSHIFT
     * constant must be set.
     *
     * @param ac IAtomContainer to search
     * @param shiftList shift value list to match
     * @param tol Tolerance value [ppm]
     * @param atomType Element name (e.g. "C") which also occurrs in 
     * {@link testkit.Utils#getNMRShiftConstant(java.lang.String)}
     * @return List of match indices for every query shift within the IAtomContainer
     */
    public static ArrayList<Integer> findShiftMatches(final IAtomContainer ac, final ArrayList<Double> shiftList, final double tol, final String atomType) {

        final ArrayList<Integer> matches = new ArrayList<>();
        for (int i = 0; i < shiftList.size(); i++) {
            matches.add(NMR.Utils.findSingleShiftMatch(ac, shiftList.get(i), tol, atomType));
        }

        return matches;
    }

    /**
     * Finds the match with the lowest deviation between a given shift value and
     * the atoms of an atom container. A tolerance value and NMRSHIFT constant
     * must be set.
     *
     * @param ac IAtomContainer to search
     * @param shift Shift value to match [ppm]
     * @param tol Tolerance value [ppm]
     * @param atomType Element name (e.g. "C") which also occurrs in 
     * {@link testkit.Utils#getNMRShiftConstant(java.lang.String)}
     * @return Match index of a query shift within the IAtomContainer
     */
    public static int findSingleShiftMatch(final IAtomContainer ac, final double shift, final double tol, final String atomType) {

        int matchIndex = -1;
        double minDiff = tol, acShift;
        for (int k = 0; k < ac.getAtomCount(); k++) {
            // skip other atom types than given
            if (ac.getAtom(k).getProperty(NMR.Utils.getNMRShiftConstant(atomType)) == null) {
                continue;
            }
            // figure out the atom with lowest shift deviation 
            acShift = Double.parseDouble(ac.getAtom(k).getProperty(NMR.Utils.getNMRShiftConstant(atomType)).toString());
            if ((shift - tol <= acShift) && (acShift <= shift + tol) && (Math.abs(shift - acShift) < minDiff)) {
                minDiff = Math.abs(shift - acShift);
                matchIndex = k;
            }
        }

        return matchIndex;
    }

    /**
     * Finds match indices between a given shift list from a peak table and an atom container.
     * Wrapper function for {@link #parsePeakTable(String, int)},
     * {@link #findShiftMatches(IAtomContainer, ArrayList, double, String)}
     * and
     * {@link #correctShiftMatches(IAtomContainer, ArrayList, ArrayList, double, String)}.
     *
     * @param ac IAtomContainer to search for matches
     * @param pathToPeakList Path to peak table
     * @param atomType Element name (e.g. "C") which also occurrs in 
     * {@link testkit.Utils#getNMRShiftConstant(java.lang.String)}
     * @param tol Tolerance value [ppm]
     * @param column Column number of shift values in peak table
     * @return Indices of matches for each shift within the IAtomContainer 
     * @throws IOException
     */
    public static ArrayList<Integer> matchShiftsFromPeakTable(final IAtomContainer ac, final String pathToPeakList, final String atomType, final double tol, final int column) throws IOException {

        final ArrayList<Double> shiftsAtomType = NMR.Utils.parsePeakTable(pathToPeakList, column);
        ArrayList<Integer> matchesAtomType = NMR.Utils.findShiftMatches(ac, shiftsAtomType, tol, atomType);
        matchesAtomType = NMR.Utils.correctShiftMatches(ac, shiftsAtomType, matchesAtomType, tol, atomType);

        return matchesAtomType;
    }
    
    
    /**
     * Finds match indices between a given shift list from a XML file and an
     * atom container. Wrapper function for {@link #parseXML(String, int)},
     * {@link #findShiftMatches(IAtomContainer, ArrayList, double, String)} and
     * {@link #correctShiftMatches(IAtomContainer, ArrayList, ArrayList, double, String)}.
     *
     * @param ac IAtomContainer to search for matches
     * @param pathToXML
     * @param atomType Element name (e.g. "C") which also occurrs in 
     * {@link testkit.Utils#getNMRShiftConstant(java.lang.String)}
     * @param tol Tolerance value [ppm]
     * @param ndim number of dimensions of given data 1 (1D) or 2 (2D)
     * @param attribute which attribute index in XML peak nodes should be used: 
     * 1 (shift of 1st dimension), 2 (shift of 2nd dimension if 2D data, 
     * intensity if 1D data) or 3 (intensity if 2D data)
     * @return Indices of matches for each shift within the IAtomContainer
     * @throws IOException
     * @throws javax.xml.parsers.ParserConfigurationException
     * @throws org.xml.sax.SAXException
     */
    public static ArrayList<Integer> matchShiftsFromXML(final IAtomContainer ac, final String pathToXML, final String atomType, final double tol, final int ndim, final int attribute) throws IOException, ParserConfigurationException, SAXException {

        final ArrayList<Double> shiftsAtomType = NMR.Utils.parseXML(pathToXML, ndim, attribute);
        ArrayList<Integer> matchesAtomType = NMR.Utils.findShiftMatches(ac, shiftsAtomType, tol, atomType);
        matchesAtomType = NMR.Utils.correctShiftMatches(ac, shiftsAtomType, matchesAtomType, tol, atomType);

        return matchesAtomType;
    }
    
    
    /**
     * Creates a two dimensional array of a given NMRShiftDB NMR entry 
     * with all shift values and atom indices.
     *
     * @param shiftsString
     * @return two dimensional array: 
     * 1. dimension: shift entry (row); 
     * 2. dimension: shift value (column 1), atom index in atom container (column 2)
     */
    public static String[][] parseShiftsNMRShiftDB(final String shiftsString){
        
        if(shiftsString.trim().length() == 0){
            return new String[][]{};
        }
        
        String[] signalSplit;
        final String[] shiftsSplit = shiftsString.split("\\|");
        final String[][] values = new String[shiftsSplit.length][3];
        for (int i = 0; i < shiftsSplit.length; i++) {
            signalSplit = shiftsSplit[i].split(";");
            values[i][0] = signalSplit[0];
            values[i][1] = signalSplit[1];
            values[i][2] = signalSplit[2];
        }
        
        return values;
    }
    
    /**
     * Returns the NMR shift constant value for a given element. As far as 
     * it is defined, the value from CDKConstants.NMRSHIFT_* (e.g.
     * {@link org.openscience.cdk.CDKConstants#NMRSHIFT_CARBON}) will be used.
     * Otherwise the same format is used for other atom types.
     * Elements defined so far: C, H, N, P, F, D, O, S, Si, B, Pt.
     * @param element element's symbol (e.g. "C")
     * @return
     */
    public static String getNMRShiftConstant(final String element){
        switch(element){
            case "C": return CDKConstants.NMRSHIFT_CARBON;
            case "H": return CDKConstants.NMRSHIFT_HYDROGEN;
            case "N": return CDKConstants.NMRSHIFT_NITROGEN;
            case "P": return CDKConstants.NMRSHIFT_PHOSPORUS;
            case "F": return CDKConstants.NMRSHIFT_FLUORINE;
            case "D": return CDKConstants.NMRSHIFT_DEUTERIUM;
            case "O": return "oxygen nmr shift";
            case "S": return "sulfur nmr shift";
            case "Si": return "silicon nmr shift";
            case "B": return "boron nmr shift";
            case "Pt": return "platinum nmr shift";
            default:
                return null;
        }
    }
    
    
    /**
     * Returns a bond type for two bond atoms from its hybridization.
     * CURRENTLY ONLY SINGLE BOND DETECTION POSSIBLE!!!
     * This function detects single, double and triple bonds and returns a 
     * bond order from {@link org.openscience.cdk.interfaces.IBond.Order}.
     * If no bond type could be detected then 
     * {@link org.openscience.cdk.interfaces.IBond.Order#UNSET} will be returned.
     * For single and double bond detection, the following elements are defined so far: C, O, N, S. 
     * For triple bond detection, the following elements are defined so far: C, N.
     *
     * 
     * @param atom1
     * @param atom2
     * @return
     */
    public static IBond.Order getBondTypeFromHybridizations(final IAtom atom1, final IAtom atom2){
        
        final String atomType1 = atom1.getSymbol();
        final IAtomType.Hybridization hybridization1 = atom1.getHybridization();
        final String atomType2 = atom2.getSymbol();
        final IAtomType.Hybridization hybridization2 = atom2.getHybridization();
        
        if(hybridization1 == null || hybridization2 == null){
            return IBond.Order.UNSET;
        }
        IBond.Order bondOrder1 = IBond.Order.UNSET;
        IBond.Order bondOrder2 = IBond.Order.UNSET;
        // single bond detection, the "3" means all SP3 hybrdidizations like SP3, SP3D2 or PLANAR3
        if ((atomType1.equals("C") || atomType1.equals("O") || atomType1.equals("N") || atomType1.equals("S"))
                && hybridization1.toString().contains("3")) {
            return IBond.Order.SINGLE;
        }
        if ((atomType2.equals("C") || atomType2.equals("O") || atomType2.equals("N") || atomType2.equals("S"))
                && hybridization2.toString().contains("3")) {
            return IBond.Order.SINGLE;
        }
//        // double bond detection
//        if ((atomType1.equals("C") && (hybridization1.equals(IAtomType.Hybridization.SP1) || hybridization1.equals(IAtomType.Hybridization.SP2)))
//                || ((atomType1.equals("O") || atomType1.equals("N") || atomType1.equals("S")) && (hybridization1.equals(IAtomType.Hybridization.SP2)))) {
//            bondOrder1 = IBond.Order.DOUBLE;
//        }
//        if ((atomType2.equals("C") && (hybridization2.equals(IAtomType.Hybridization.SP1) || hybridization2.equals(IAtomType.Hybridization.SP2)))
//                || ((atomType2.equals("O") || atomType2.equals("N") || atomType2.equals("S")) && hybridization2.equals(IAtomType.Hybridization.SP2))) {
//            bondOrder2 = IBond.Order.DOUBLE;
//        }
//        // triple bond detection
//        if ((atomType1.equals("C") && (hybridization1.equals(IAtomType.Hybridization.SP1)))
//                && (atomType2.equals("N") && hybridization2.equals(IAtomType.Hybridization.SP1))) {
//            bondOrder1 = IBond.Order.TRIPLE;
//        }
//        if ((atomType2.equals("N") && (hybridization2.equals(IAtomType.Hybridization.SP1)))
//                && (atomType1.equals("C") && hybridization1.equals(IAtomType.Hybridization.SP1))) {
//            bondOrder2 = IBond.Order.TRIPLE;
//        }
        
        if (bondOrder1.equals(bondOrder2)) {
            return bondOrder1;
        }
        
        return IBond.Order.UNSET;
    }
    
    /**
     * Returns a list of open bonds of an atom.
     *
     * @param ac atom container
     * @param atomIndex index of the atom to test
     * @return 
     */
    public static ArrayList<IBond.Order> getOpenBonds(final IAtomContainer ac, final int atomIndex){
        
        final IAtom atom = ac.getAtom(atomIndex);
        if(atom.getHybridization() == null){
            return null;
        }
        final ArrayList<IBond.Order> bondOrderList = new ArrayList<>();
        final AtomValenceDescriptor valenceDesc = new AtomValenceDescriptor();
        final int valence = Integer.valueOf(valenceDesc.calculate(atom, ac).getValue().toString());
        int electronsLeft = (8 - (valence + atom.getImplicitHydrogenCount()));
        
        if (electronsLeft == 0) {
//            System.out.println(atom.getSymbol() + ": " + atomIndex + " (" + atom.getHybridization() + "): " + bondOrderList);
            return bondOrderList;
        }
        // only one single bond left; possible at SP1, SP2 and SP3
        if (electronsLeft == 1) {
            bondOrderList.add(IBond.Order.SINGLE);
//            System.out.println(atom.getSymbol() + ": " + atomIndex + " (" + atom.getHybridization() + "): " + bondOrderList);
            return bondOrderList;
        }
        // with SP3 are only single bonds possible
        if (atom.getHybridization().equals(IAtomType.Hybridization.SP3)) {
            // subtract the single bonded neighbor number
            electronsLeft -= ac.getConnectedAtomsList(atom).size();
            for (int k = 0; k < electronsLeft; k++) {
                bondOrderList.add(IBond.Order.SINGLE);
            }
//            System.out.println(atom.getSymbol() + ": " + atomIndex + " (" + atom.getHybridization() + "): " + bondOrderList);
            return bondOrderList;
        }
        
        if (atom.getHybridization().equals(IAtomType.Hybridization.SP2)) {
            switch (atom.getSymbol()) {
                case "O":
                case "S":
                    bondOrderList.add(IBond.Order.DOUBLE);
                    return bondOrderList;
                case "C":
                    bondOrderList.add(IBond.Order.SINGLE);
                    bondOrderList.add(IBond.Order.SINGLE);
                    bondOrderList.add(IBond.Order.DOUBLE);
                    break;
                case "N":
                    bondOrderList.add(IBond.Order.SINGLE);
                    bondOrderList.add(IBond.Order.DOUBLE);
                    break;
                default:
                    break;
            }
        } else if (atom.getHybridization().equals(IAtomType.Hybridization.SP1)) {
            switch (atom.getSymbol()) {
                case "C":
                    bondOrderList.add(IBond.Order.DOUBLE);
                    bondOrderList.add(IBond.Order.DOUBLE);
                    // or
                    bondOrderList.add(IBond.Order.SINGLE);
                    bondOrderList.add(IBond.Order.TRIPLE);
                    break;
                case "N":
                    bondOrderList.add(IBond.Order.TRIPLE);
                    break;
                default:
                    break;
            }
        }
        for (IAtom neighbor : ac.getConnectedAtomsList(atom)) {
            bondOrderList.remove(ac.getBond(atom, neighbor).getOrder());
            electronsLeft -= NMR.Utils.getElectronNumberByBondOrder(ac.getBond(atom, neighbor).getOrder());
        }
        
        int theoCounter = 0;
        for (IBond.Order order : bondOrderList) {
            theoCounter += NMR.Utils.getElectronNumberByBondOrder(order);
        }
        
        switch(Math.abs(theoCounter - electronsLeft)){
            case 1:
                bondOrderList.remove(IBond.Order.SINGLE);
                theoCounter -= 1;
                break;
            case 2:
                
                break;
            case 3:
                
                break;
        }
        
        
//        System.out.println(atom.getSymbol() + ": " + atomIndex + " (" + atom.getHybridization() + "): " + bondOrderList + " -> e: " + theoCounter + " (theo) vs. " + electronsLeft + " (real), bond counter: " + ac.getConnectedAtomsList(atom).size() + " (+" + atom.getImplicitHydrogenCount() + "H)");
        
    
        return bondOrderList;
    }
    
    
    public static int getElectronNumberByBondOrder(final IBond.Order order) {
        switch (order) {
            case SINGLE:
                return 1;
            case DOUBLE:
                return 2;
            case TRIPLE:
                return 3;
            case QUADRUPLE:
                return 4;    
            case QUINTUPLE:
                return 5;
            case SEXTUPLE:
                return 6;    
            default:
                return 0;
        }
    }
    
    
    /**
     * Returns the NMR isotope identifier for a given element, e.g. C -> 13C. 
     * Elements defined so far: C, H, N, P, F, D, O, S, Si, B, Pt.
     *
     * @param element element's symbol (e.g. "C")
     * @return
     */
    public static String getNMRIsotopeIdentifier(final String element) {
        switch(element){
            case "C": return "13C";
            case "H": return "1H";
            case "N": return "15N";
            case "P": return "31P";
            case "F": return "19F";
            case "O": return "17O";
            case "S": return "33S";
            case "Si": return "29Si";
            case "B": return "11B";
            case "Pt": return "195Pt";
            default:
                return null;
        }
    }
    
//    /**
//     * Returns the hybridization level of each heavy atom in given molecule which has
//     * its own shift value.
//     * First it compares the number of attached (implicit) hydrogens and sets 
//     * the hybridization level from it directly. This is only possible for 
//     * carbons with three or four attached hydrogens (sp3). [CURRENTLY DISABLED]
//     * 
//     * If less than three hydrogens are attached or in case of other heavy 
//     * atoms then a NMRShiftDB file will be used to obtain the
//     * frequencies of the different hybridization levels from the database. 
//     * This happens for directly bonded neighbors too.
//     * 
//     *
//     * @param ac
//     * @param pathToNMRShiftDB
//     * @param tol
//     * @param molFormula
//     * @return
//     * @throws FileNotFoundException
//     */
//    public static HashMap<String, HashMap<Integer, HashMap<String, ArrayList<Integer>>>> getHybridizationsFromNMRShiftDB(final IAtomContainer ac, final String pathToNMRShiftDB, final double tol, final IMolecularFormula molFormula) throws FileNotFoundException{
//        
//        final HashMap<Integer, HashMap<String, ArrayList<Integer>>> elementsHybridCounter = new HashMap<>();
//        final HashMap<Integer, HashMap<String, ArrayList<Integer>>> elementsBondTypeCounter = new HashMap<>();
//        final HashMap<Integer, HashMap<String, Integer>> expactedNeighbors = new HashMap<>();
//        String NMRSHIFT_ATOMTYPE;
//        // initializations only
//        for (int i = 0; i < ac.getAtomCount(); i++) {
//            // sure case for carbon: 3 or 4 hydrogens -> sp3
////            if (ac.getAtom(i).getSymbol().equals("C") && ac.getAtom(i).getImplicitHydrogenCount() >= 3) {
////                ac.getAtom(i).setHybridization(IAtomType.Hybridization.SP3);
////                continue;
////            }
//            NMRSHIFT_ATOMTYPE = testkit.Utils.getNMRShiftConstant(ac.getAtom(i).getSymbol());
//            // is the NMR shift constant defined and does the nmr shift property entry in an atom exist?
//            if ((NMRSHIFT_ATOMTYPE == null) || (ac.getAtom(i).getProperty(NMRSHIFT_ATOMTYPE) == null)) {
//                continue;
//            }
//            elementsHybridCounter.put(i, new HashMap<>());
//            elementsBondTypeCounter.put(i, new HashMap<>());
//            elementsHybridCounter.get(i).put("query", new ArrayList<>());
//            elementsHybridCounter.get(i).put("queryH", new ArrayList<>());
//            // create an array list for each atom type in given molecular formula
//            for (IElement elem : MolecularFormulaManipulator.getHeavyElements(molFormula)) {
//                elementsHybridCounter.get(i).put(elem.getSymbol(), new ArrayList<>());
//                elementsBondTypeCounter.get(i).put(elem.getSymbol(), new ArrayList<>());
//            }
//
//            expactedNeighbors.put(i, new HashMap<>());
//            for (IAtom expNeighbor : ac.getConnectedAtomsList(ac.getAtom(i))) {
//                if (!expactedNeighbors.get(i).keySet().contains(expNeighbor.getSymbol())) {
//                    expactedNeighbors.get(i).put(expNeighbor.getSymbol(), 0);
//                }
//                expactedNeighbors.get(i).put(expNeighbor.getSymbol(), expactedNeighbors.get(i).get(expNeighbor.getSymbol()) + 1);
//            }
//        }
//        // beginning of DB search
//        String shiftsDB;
//        double shiftDB, shiftQ;
//        int atomIndexDB;
//        boolean add, toContinue;
//        final AtomHybridizationDescriptor hybridDesc = new AtomHybridizationDescriptor();
//        IAtom qAtom;
//        IAtomContainer acDB;
//        final IteratingSDFReader iterator = new IteratingSDFReader(
//                new FileReader(pathToNMRShiftDB),
//                SilentChemObjectBuilder.getInstance()
//        );
//        while (iterator.hasNext()) {
//            acDB = iterator.next();
//            ArrayList<String> props = (ArrayList<String>) (ArrayList<?>) (new ArrayList<>(acDB.getProperties().keySet()));
//            Collections.sort(props);
//            // the DB entry should at least contain one carbon spectrum 
//            toContinue = false;
//            for (String prop : props) {
//                if (prop.contains("Spectrum " + testkit.Utils.getNMRIsotopeIdentifier("C"))) {
//                    toContinue = true;
//                    break;
//                }
//            }
//            if (!toContinue) {
//                continue;
//            }
//
//            for (int i : elementsHybridCounter.keySet()) {
//                qAtom = ac.getAtom(i);
//                // check wether the DB entry contains a spectrum for the current query atom type
//                shiftsDB = null;
//                for (String prop : props) {
//                    if (prop.contains("Spectrum " + testkit.Utils.getNMRIsotopeIdentifier(qAtom.getSymbol()))) {
//                        shiftsDB = acDB.getProperty(prop);
//                        break;
//                    }
//                }
//                if(shiftsDB == null){
//                    continue;
//                }
//                // ignore the already set sp3 hybridizations at carbon atoms with at least 3 implicit hydrogens
////                    if (qAtom.getSymbol().equals("C") && qAtom.getImplicitHydrogenCount() >= 3) {
////                        continue;
////                    }
//                shiftQ = qAtom.getProperty(testkit.Utils.getNMRShiftConstant(ac.getAtom(i).getSymbol()));
//                
//                // check wether the DB entry contains a proton spectrum
//                String shiftsDBHydrogen = null;
//                for (String prop : props) {
//                    if (prop.contains("Spectrum " + testkit.Utils.getNMRIsotopeIdentifier("H"))) {
//                        shiftsDBHydrogen = acDB.getProperty(prop);
//                        break;
//                    }
//                }
//                
//                String[][] shiftsDBvalues = testkit.Utils.parseShiftsNMRShiftDB(shiftsDB);
//                for (String[] shiftsDBvalue : shiftsDBvalues) {
//                    shiftDB = Double.parseDouble(shiftsDBvalue[0]);
//                    atomIndexDB = Integer.parseInt(shiftsDBvalue[2]);
//                    add = true;
//                    // shift match within a shift tolerance range
//                    if ((shiftQ - tol <= shiftDB) && (shiftDB <= shiftQ + tol)) {
//                        // matched atom should have the same number of attached (implicit) hydrogens
//                        if (acDB.getAtom(atomIndexDB).getImplicitHydrogenCount().intValue() == qAtom.getImplicitHydrogenCount().intValue()) {
//                            // count next neighbors
//                            HashMap<String, Integer> foundNeighbors = new HashMap<>();
//                            for (IAtom neighborAtomDB : acDB.getConnectedAtomsList(acDB.getAtom(atomIndexDB))) {
//                                if (!foundNeighbors.keySet().contains(neighborAtomDB.getSymbol())) {
//                                    foundNeighbors.put(neighborAtomDB.getSymbol(), 0);
//                                }
//                                foundNeighbors.put(neighborAtomDB.getSymbol(), foundNeighbors.get(neighborAtomDB.getSymbol()) + 1);
//                            }
//                            // check whether the number of expacted next neighbors is higher than the number of found next neighbor, if yes then skip this DB atom match
//                            for (String elemExpNeighbor : expactedNeighbors.get(i).keySet()) {
//                                if (foundNeighbors.get(elemExpNeighbor) == null || (expactedNeighbors.get(i).get(elemExpNeighbor) > foundNeighbors.get(elemExpNeighbor))) {
//                                    add = false;
//                                }
//                            }
//                            if(!add){
//                                continue;
//                            }
//                            // only elements which occur in molecular formula of the unknown are allowed, otherwise skip this matched DB atom
//                            for (IAtom neighborAtomDB : acDB.getConnectedAtomsList(acDB.getAtom(atomIndexDB))) {
//                                if (MolecularFormulaManipulator.getElementCount(molFormula, neighborAtomDB.getSymbol()) == 0) {
//                                    add = false;
//                                    break;
//                                }
//                                // ignore explicit protons; ignore query atoms here, add them as below -> otherwise multiple counting 
//                                if (!neighborAtomDB.getSymbol().equals("H")){// && !neighborAtomDB.getSymbol().equals(qAtom.getSymbol())) {
//                                    elementsHybridCounter.get(i).get(neighborAtomDB.getSymbol()).add(Integer.parseInt(hybridDesc.calculate(neighborAtomDB, acDB).getValue().toString()));
//                                    elementsBondTypeCounter.get(i).get(neighborAtomDB.getSymbol()).add(acDB.getBond(acDB.getAtom(atomIndexDB), neighborAtomDB).getOrder().numeric());
//                                }
//                            }
//                            if(!add){
//                                continue;
//                            }
//                            // likely allowed to add hybridization for query atom
//                            // check whether the shifts of attached hydrogens are equal to hydrogen shifts of query atom -> higher priority at hybridization assignment step later
//                            boolean added = false;
//                            if(shiftsDBHydrogen != null){
//                                String[][] shiftsDBvaluesHydrogen = testkit.Utils.parseShiftsNMRShiftDB(shiftsDBHydrogen);
//                                if(qAtom.getProperty("HydrogenShifts") != null){
//                                    ArrayList<Double> shiftsQAtomvaluesHydrogen = qAtom.getProperty("HydrogenShifts");
//                                    for (int j = 0; j < shiftsQAtomvaluesHydrogen.size(); j++) {
//                                        for (String[] shiftsDBvalueHydrogen : shiftsDBvaluesHydrogen) {
//                                            shiftDB = Double.parseDouble(shiftsDBvalueHydrogen[0]);
//                                            if((shiftsQAtomvaluesHydrogen.get(j) - 0.1 <= shiftDB) && (shiftDB <= shiftsQAtomvaluesHydrogen.get(j) + 0.1)){
//                                                elementsHybridCounter.get(i).get("queryH").add(Integer.parseInt(hybridDesc.calculate(acDB.getAtom(atomIndexDB), acDB).getValue().toString()));
//                                                added = true;
//                                                break;
//                                            }
//                                        }
//                                        if(added){
//                                            break;
//                                        }
//                                    }
//                                }
//
//                            } else {
//                                elementsHybridCounter.get(i).get("query").add(Integer.parseInt(hybridDesc.calculate(acDB.getAtom(atomIndexDB), acDB).getValue().toString()));
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        final HashMap<String, HashMap<Integer, HashMap<String, ArrayList<Integer>>>> toReturn = new HashMap<>();
//        toReturn.put("hybridCounter", elementsHybridCounter);
//        toReturn.put("bondTypeCounter", elementsBondTypeCounter);
//        
//        return toReturn;
//    }
    
    
    
    
    
    public static int[] getNeighborhoodBondsCount(final IAtomContainer ac, final int indexAC, final String[] bondsSet, final String[] neighborElems){
        final int[] counts = new int[neighborElems.length * bondsSet.length];
        String foundBonds;
        // for all given neighbor element types
        for (int n = 0; n < neighborElems.length; n++) {
            foundBonds = "";
            // for all next neighbors of a specific element
            for (IAtom neighborAtom : ac.getConnectedAtomsList(ac.getAtom(indexAC))) {
                // skip if not the right neighborhood element or bond type is unknown/unset
                if ((!neighborAtom.getSymbol().equals(neighborElems[n])) || (NMR.Utils.getStringFromBondOrder(ac.getBond(ac.getAtom(indexAC), neighborAtom).getOrder()) == null)) {
                    continue;
                }
                foundBonds += NMR.Utils.getStringFromBondOrder(ac.getBond(ac.getAtom(indexAC), neighborAtom).getOrder());
            }
            for (int k = 0; k < bondsSet.length; k++) {
                counts[n * bondsSet.length + k] = 0;
                if (NMR.Utils.sortString(foundBonds).equals(NMR.Utils.sortString(bondsSet[k]))) {
                    counts[n * bondsSet.length + k] = 1;
                    break;
                }
            }
        }
        
        return counts;
    }
    
    
    public static void writeNeighborhoodBondsCountMatrix(final String pathToOutput, final int[][] m, final String[] bondsSet, final String elem, String[] neighborElems, final int min, final int max, final int stepSize) throws IOException{
        
        if(stepSize < 1){
            System.err.println("stepSize < 1 not allowed!!!");
            return;
        }
        
        final StringBuilder sb = new StringBuilder();
        sb.append("shift [" + elem + "] (" + stepSize + "),nTotal,inRing,isArom,q" + elem + "," + elem + "H," + elem + "H2," + elem + "H3,");
        for (int i = 0; i < neighborElems.length; i++) {
            for (int j = 0; j < bondsSet.length; j++) {
                sb.append(bondsSet[j] + "[" + neighborElems[i] + "]");
                if (j < bondsSet.length - 1) {
                    sb.append(",");
                }
            }
            if (i < neighborElems.length - 1) {
                sb.append(",");
            }
        }
        sb.append("\n");
        for (int i = 0; i < stepSize * (max - min) + 1; i++) {
            sb.append((i + min) + ",");
            for (int j = 0; j < 3 + 4 + neighborElems.length * bondsSet.length; j++) {
                sb.append(m[i][j]);
                if (j < 3 + 4 + neighborElems.length * bondsSet.length - 1) {
                    sb.append(",");
                }
            }
            sb.append("\n");
        }
        
        final FileWriter writer = new FileWriter(pathToOutput);
        writer.append(sb.toString());
        writer.flush();
        writer.close();
    }
    
    
    public static String sortString(final String s) {
        final char[] c = s.toCharArray();
        Arrays.sort(c);
        return new String(c);
    }
    
    
    public static ArrayList<ArrayList<IBond.Order>> getBondOrderSets(final String[] valences) {
        
        final ArrayList<ArrayList<IBond.Order>> bondOrderSets = new ArrayList<>();
        for (int i = 0; i < valences.length; i++) {
            bondOrderSets.add(new ArrayList<>());
            for (int k = 0; k < StringUtils.countMatches(valences[i], "-"); k++) {
                bondOrderSets.get(i).add(IBond.Order.SINGLE);
            }
            for (int k = 0; k < StringUtils.countMatches(valences[i], "="); k++) {
                bondOrderSets.get(i).add(IBond.Order.DOUBLE);
            }
            for (int k = 0; k < StringUtils.countMatches(valences[i], "%"); k++) {
                bondOrderSets.get(i).add(IBond.Order.TRIPLE);
            }
        }
        
        return bondOrderSets;
    }
    
    
    public static String getStringFromBondOrder(final IBond.Order order) {
        switch (order) {
            case SINGLE:
                return "-";
            case DOUBLE:
                return "=";
            case TRIPLE:
                return "%";
            default:
                return null;
        }
    }
    
    
    public static IBond.Order getBondOrderFromString(final String order){
        switch(order){
            case "-": return IBond.Order.SINGLE;
            case "=": return IBond.Order.DOUBLE;
            case "%": return IBond.Order.TRIPLE;
            default:  return null;
        }
    }
    
    
    public static void writeCSV(final String pathToOutput, final String table) throws IOException {
        FileWriter fr = new FileWriter(new File(pathToOutput));
        BufferedWriter br = new BufferedWriter(fr);
        br.write(table);
        br.close();
    }
    
    
    /**
     *
     * @param input
     * @return
     */
    public static ArrayList<Integer> getOutliers(ArrayList<Integer> input) {
        final ArrayList<Integer> outliers = new ArrayList<>();
        if(input.size() <= 1){
            return outliers;
        }
        Collections.sort(input);
        final List<Integer> data1 = input.subList(0, input.size() / 2);
        final List<Integer> data2;
        if (input.size() % 2 == 0) {
            data2 = input.subList(input.size() / 2, input.size());
        } else {
            data2 = input.subList(input.size() / 2 + 1, input.size());
        }
        final double q1 = getMedian(data1);
        final double q3 = getMedian(data2);
        final double iqr = q3 - q1;
        final double lowerFence = q1 - 1.5 * iqr;
        final double upperFence = q3 + 1.5 * iqr;
        for (int i = 0; i < input.size(); i++) {
            if ((input.get(i) < lowerFence) || (input.get(i) > upperFence)) {
                outliers.add(input.get(i));
            }
        }
//        System.out.println("input size: " + input.size());
//        System.out.println("output size: " + outliers.size());
        return outliers;
    }

    /**
     *
     * @param data
     * @return
     */
    public static double getMedian(final List<Integer> data) {
        if(data.size() == 1){
            return data.get(0);
        }
        Collections.sort(data);
        if (data.size() % 2 == 1) {
            return data.get(data.size() / 2);
        } else {
            return (data.get(data.size() / 2 - 1) + data.get(data.size() / 2)) / 2.0;
        }
    }
    
    
    /**
     *
     * @param data
     * @return
     */
    public static double getMedian(final ArrayList<Double> data) {
        if(data.size() == 1){
            return data.get(0);
        }
        Collections.sort(data);
        if (data.size() % 2 == 1) {
            return data.get(data.size() / 2);
        } else {
            return (data.get(data.size() / 2 - 1) + data.get(data.size() / 2)) / 2.0;
        }
    }
    
    
    /**
     *
     * @param data
     * @return
     */
    public static double getRMS(final ArrayList<Double> data) {
        if(data.size() == 1){
            return data.get(0);
        }
        double qSum = 0;
        for (final Double d : data) {
            qSum += d*d;
        }
        
        return Math.sqrt(qSum/data.size());
    }
    
    
    /**
     *
     * @param data
     * @return
     */
    public static double getMean(final ArrayList<Double> data) {
        double sum = 0;
        for (Double d : data) {
            sum += d;
        }
        return sum/data.size();
    }
    
    
    /**
     * Tests whether two array lists of integers are equal which also means
     * bidirectional values to each other.
     *
     * @param shiftMatches1
     * @param shiftMatches2
     * @return
     */
    public static boolean isBidirectional(final ArrayList<Integer> shiftMatches1, final ArrayList<Integer> shiftMatches2) {
        final ArrayList<Integer> temp1 = new ArrayList<>(shiftMatches1);
        final ArrayList<Integer> temp2 = new ArrayList<>(shiftMatches2);
        Collections.sort(temp1);
        Collections.sort(temp2);

        return temp1.equals(temp2);
    }

    /**
     *
     * @param ac
     * @param shiftMatches1
     * @param shiftMatches2
     * @param prop
     */
    public static void setBidirectionalLinks(final IAtomContainer ac, final ArrayList<Integer> shiftMatches1, final ArrayList<Integer> shiftMatches2, final String prop) {

        ArrayList<Integer> propList1, propList2;
        for (int i = 0; i < shiftMatches1.size(); i++) {
            if (shiftMatches1.get(i) >= 0 && shiftMatches2.get(i) >= 0) {
                if (ac.getAtom(shiftMatches1.get(i)).getProperty(prop) == null) {
                    ac.getAtom(shiftMatches1.get(i)).setProperty(prop, new ArrayList<>());
                }
                if (ac.getAtom(shiftMatches2.get(i)).getProperty(prop) == null) {
                    ac.getAtom(shiftMatches2.get(i)).setProperty(prop, new ArrayList<>());
                }
                propList1 = ac.getAtom(shiftMatches1.get(i)).getProperty(prop);
                propList2 = ac.getAtom(shiftMatches2.get(i)).getProperty(prop);
                if (!propList1.contains(shiftMatches2.get(i))) {
                    propList1.add(shiftMatches2.get(i));
                }
                if (!propList2.contains(shiftMatches1.get(i))) {
                    propList2.add(shiftMatches1.get(i));
                }
            }
        }
    }
    
    
    public static ArrayList<Double> countSetShiftInAtomContainer(final IAtomContainer ac, final ArrayList<Integer> indices){
        
        final ArrayList<Double> shifts = new ArrayList<>();
        for (final Integer index : indices) {
                shifts.add(ac.getAtom(index).getProperty(Utils.getNMRShiftConstant(ac.getAtom(index).getSymbol())));
        }
        return shifts;
    }
    
    
    
    public static String getFileFormat(final String pathToFile) {

        if(pathToFile == null || pathToFile.trim().isEmpty()){
            return "";
        }        
        final String[] split = pathToFile.split("\\.");

        return split[split.length - 1];
    }
//    
//    public static IAtomContainer getStructureFromINCHICode(final String inchi) throws CDKException {
//        final InChIToStructure intostruct = InChIGeneratorFactory.getInstance().getInChIToStructure(
//                inchi, SilentChemObjectBuilder.getInstance()
//        );
//
//        INCHI_RET ret = intostruct.getReturnStatus();
//        if (ret == INCHI_RET.WARNING) {
//        // Structure generated, but with warning message
//            System.out.println("InChI warning: " + intostruct.getMessage());
//        } else if (ret != INCHI_RET.OKAY) {
//        // Structure generation failed
//            throw new CDKException("Structure generation failed: " + ret.toString()
//                    + " [" + intostruct.getMessage() + "]");
//        }
//
//        System.out.println("inchi ac: " + intostruct.getAtomContainer().getAtomCount());
//    
//        return SilentChemObjectBuilder.getInstance().newAtomContainer();//intostruct.getAtomContainer();
//    
//    }
    
   
    
    
    
    /**
     *
     * @param lookup
     * @return
     */
    public static HashMap<String, Double> getRMS(final HashMap<String, ArrayList<Double>> lookup){
        
        final HashMap<String, Double> rms = new HashMap<>();
        for (final String key : lookup.keySet()) {
            rms.put(key, NMR.Utils.getRMS(lookup.get(key)));
//            System.out.println("count: " + lookup.get(key).size() + ", mean: " + NMR.Utils.getMean(lookup.get(key)) + ", rms: " +  rms.get(key) + ", median: " + NMR.Utils.getMedian(lookup.get(key)));
        }
      
        return rms;
    }
    
    
    
    
    
    
    
    
    public static IAtomContainer setAromaticitiesInAtomContainer(final IAtomContainer ac, final int maxCycleSize) throws CDKException {
        
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(ac);
        final ElectronDonation model = ElectronDonation.cdkAllowingExocyclic();
        final CycleFinder cycles = Cycles.or(Cycles.all(), Cycles.all(maxCycleSize));
        final Aromaticity aromaticity = new Aromaticity(model, cycles);
        aromaticity.apply(ac);
        
        return ac;
    }
    
    
    public static double getTanimotoCoefficient(final IAtomContainer a, final IAtomContainer b) throws CDKException, IOException, CloneNotSupportedException{
        
        // pubchem fingerprinter expects 
        // 1. explicit hydrogens
//        AtomContainerManipulator.convertImplicitToExplicitHydrogens(a);
        // 2. set atom type names -> done during setting of aromaticities
        // 3. set aromaticity -> done during DB scanning
        
        SubstructureFingerprinter substructfp  = new SubstructureFingerprinter();
        IBitFingerprint fingerprint = substructfp.getBitFingerprint(a);
        System.out.println("\n\ndefault substructure bitstring: " + fingerprint.asBitSet());
        for (int setbit : fingerprint.getSetbits()) {
            System.out.println("default substructure of index " + setbit + ": " + substructfp.getSubstructure(setbit));// + " -> " + SMARTSParser.parse(substructfp.getSubstructure(setbit), SilentChemObjectBuilder.getInstance()));
        }
        
//        DepictionGenerator dg = new DepictionGenerator().withSize(800, 800).withAtomColors().withAtomValues().withMolTitle().withFillToFit();
        IAtomContainer ac = SilentChemObjectBuilder.getInstance().newAtomContainer();
        QueryAtomContainer qac = SMARTSParser.parse(substructfp.getSubstructure(fingerprint.getSetbits()[1]), SilentChemObjectBuilder.getInstance());
        
        System.out.println("qac: " + qac.getAtomCount() + ", " + qac.getBondCount() + " -> " + qac.getProperties());
        for (IAtom atom : qac.atoms()) {
            ac.addAtom(atom);
            System.out.println("qac atom: ");
        }
        for (IBond bond : qac.bonds()) {
            ac.addBond(bond);
            System.out.println("qac bond: " + bond);
        }
        System.out.println("ac: " + ac.getAtomCount() + ", " + ac.getBondCount() + " -> " + ac.getProperties());
                
        
//        dg.depict(ac).writeTo("/Users/mwenk/Downloads/test.png");
        
        System.out.println("\n\n");
        SubstructureFingerprinter klekotasubstructfp = new KlekotaRothFingerprinter();
        fingerprint = klekotasubstructfp.getBitFingerprint(a);
        System.out.println("Klekota substructure bitstring: " + fingerprint.asBitSet());
        for (int setbit : fingerprint.getSetbits()) {
            System.out.println("Klekota substructure of index " + setbit + ": " + klekotasubstructfp.getSubstructure(setbit));
        }
        
        
        return 0.0;//Tanimoto.calculate(pubchemfp.getBitFingerprint(a), pubchemfp.getBitFingerprint(b));
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
     /**
     * Returns 
     *
     * @param values
     * @return
     */
    public static HashMap<Integer, Double> getValueFrequencies(final ArrayList<Integer> values) {
        
        final HashMap<Integer, Double> freqs = new HashMap<>();
        final HashSet<Integer> valueLevels = new HashSet<>(values);
        int sum = 0;
        for (int value : valueLevels) {
            sum += Collections.frequency(values, value);
        }
        for (int value : valueLevels) {            
            freqs.put(value, (Collections.frequency(values, value) / (double) sum));
        }
        
        return freqs;
    }
    
    
    
    
  
    
}
