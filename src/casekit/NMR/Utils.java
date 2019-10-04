/*
 * The MIT License
 *
 * Copyright 2019 Michael Wenk [https://github.com/michaelwenk]
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
package casekit.NMR;


import casekit.NMR.model.Spectrum;
import org.apache.commons.lang3.StringUtils;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Utils {
        
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
    public static void SDFtoLSD(final String pathSDF, final String pathOut, final String pathMol2ab) throws FileNotFoundException, CDKException, IOException{
        
        
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
     * Returns a hashmap constisting of lists of atom indices in an atom container.
     * This is done for all atom types (e.g. C or Br) in given atom container.
     *
     * @param ac IAtomContainer to look in 
     * @return
     * @see #getAtomTypeIndicesByElement(org.openscience.cdk.interfaces.IAtomContainer, java.lang.String) 
     */
    public static HashMap<String, ArrayList<Integer>> getAtomTypeIndices(final IAtomContainer ac) {
        
        final HashMap<String, ArrayList<Integer>> atomTypeIndices = new HashMap<>();
        final HashSet<String> atomTypes = new HashSet<>();
        for (final IAtom heavyAtom : AtomContainerManipulator.getHeavyAtoms(ac)) {
            atomTypes.add(heavyAtom.getSymbol());
        }
        for (final String atomType : atomTypes) {
            atomTypeIndices.put(atomType, Utils.getAtomTypeIndicesByElement(ac, atomType));
        }
        
        return atomTypeIndices;
    }
    
    
    /**
     * Returns a list of atom indices in an atom container for a given atom 
     * type (e.g. C or Br)
     *
     * @param ac IAtomContainer to use for search
     * @param atomType Atom type to find in atom container 
     * @return
     */
    public static ArrayList<Integer> getAtomTypeIndicesByElement(final IAtomContainer ac, final String atomType){
        
        final ArrayList<Integer> indices = new ArrayList<>();
        for (int i = 0; i < ac.getAtomCount(); i++) {
            if(ac.getAtom(i).getSymbol().equals(atomType)){
                indices.add(i);
            }
        }
        
        return indices;
    }


    public static String getAtomTypeFromSpectrum(final Spectrum spectrum, final int dim){
        if(spectrum.containsDim(dim)){
            return Utils.getAtomTypeFromNucleus(spectrum.getNuclei()[dim]);
        }
        
        return null;
    }
    
    public static String getAtomTypeFromNucleus(final String nucleus){
        final String[] nucleusSplit = nucleus.split("\\d");
        return nucleusSplit[nucleusSplit.length - 1];
    }
    
    public static int getDifferenceSpectrumSizeAndMolecularFormulaCount(final Spectrum spectrum, final IMolecularFormula molFormula, final int dim) throws CDKException{
        if(!spectrum.containsDim(dim)){
            throw new CDKException(Thread.currentThread().getStackTrace()[2].getClassName() + "." + Thread.currentThread().getStackTrace()[2].getMethodName() + ": invalid dimension in spectrum given");
        }
        final String atomType = Utils.getAtomTypeFromSpectrum(spectrum, dim);
        int atomsInMolFormula = 0;
        if(molFormula != null){
            atomsInMolFormula = MolecularFormulaManipulator.getElementCount(molFormula, atomType);
        }
        return atomsInMolFormula - spectrum.getSignalCount();
    }
    
    public static void editSignalsInSpectrum(final Spectrum spectrum, final IMolecularFormula molFormula, final int dim) throws Exception {
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in)); int n;
        final ArrayList<Integer> validIndices = new ArrayList<>();
        int diff = Utils.getDifferenceSpectrumSizeAndMolecularFormulaCount(spectrum, molFormula, dim);
        // walk through all signals in spectrum add missing or to remove signals
        while (diff != 0) {
            // display all selectable signal indices in spectrum
            if(diff > 0){
                System.out.println("\n" + diff + " " + spectrum.getNuclei()[0] + " signals are missing!\nWhich signal is not unique?");
            } else {
                System.out.println("\n" + (-1 * diff) + " " + spectrum.getNuclei()[0] + " signals are to be removed!\nWhich signal is to remove?");
            }
            for (int s = 0; s < spectrum.getSignalCount(); s++) {
                System.out.print("index: " + s);
                for (int d = 0; d < spectrum.getNDim(); d++) {
                    System.out.print(", shift dim " + (d+1) + ": " + spectrum.getShift(s, d));
                }
                System.out.println();
                validIndices.add(s);
            }
            // get selected index by user input
            n = -1;
            while(!validIndices.contains(n)){
                System.out.println("Enter the index: ");
                n = Integer.parseInt(br.readLine());
            }
            // add/remove signals in spectrum
            if(diff > 0){
                spectrum.addSignal(spectrum.getSignal(validIndices.indexOf(n)).getClone());
                spectrum.setEquivalence(spectrum.getSignalCount() - 1, validIndices.indexOf(n));
            } else {
                spectrum.removeSignal(validIndices.indexOf(n));
            }
            diff = Utils.getDifferenceSpectrumSizeAndMolecularFormulaCount(spectrum, molFormula, dim);
        }
    }
    
    /**
     * Specified for carbons only -> not generic!!!
     *
     * @param mult
     * @return
     */
    public static Integer getHydrogenCountFromMultiplicity(final String mult){
        
        if(mult == null){
            System.out.println("null!!!");
            return null;
        }
        switch(mult){
            case "Q": 
                return 3;
            case "T": 
                return 2;
            case "D": 
                return 1;
            case "S": 
                return 0;
            default: 
                System.out.println("unknown symbol!!");
                return null;
        }
    }
    
    /**
     * Specified for carbons only -> not generic!!!
     *
     * @param hCount
     * @return
     */
    public static String getMultiplicityFromHydrogenCount(final int hCount) {
        switch (hCount) {
            case 0:
                return "S";
            case 1:
                return "D";
            case 2:
                return "T";
            case 3:
                return "Q";
            default:
                return null;
        }
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
//            case "S": return CDKConstants.NMRSHIFT_SULFUR;
            default:
                return null;
        }
    }
    
    /**
     * Returns the NMR isotope identifier for a given element, e.g. C -> 13C. 
     * Elements defined so far: C, H, N, P, F, D, O, S, Si, B, Pt.
     *
     * @param element element's symbol (e.g. "C")
     * @return
     */
    public static String getIsotopeIdentifier(final String element) {
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
                return element;
        }
    }    
    
    
    public static HashSet<String> getAtomTypesInAtomContainer(final IAtomContainer ac) {
        final HashSet<String> atomTypes = new HashSet<>();
        for (IAtom atom : ac.atoms()) {
            atomTypes.add(atom.getSymbol());
        }
        
        return atomTypes;
    }
    
    
    public static boolean checkMinMaxValue(final double min, final double max, final double value){
        
        return (value >= min && value <= max);
    }
    
    /**
     *
     * @param ac
     * @param indexAC
     * @param bondsSet
     * @param neighborElems
     * @return
     * 
     * @deprecated 
     */
    public static int[] getNeighborhoodBondsCount(final IAtomContainer ac, final int indexAC, final String[] bondsSet, final ArrayList<String> neighborElems){
        final int[] counts = new int[neighborElems.size() * bondsSet.length];
        String foundBonds;
        // for all given neighbor element types
        for (int n = 0; n < neighborElems.size(); n++) {
            foundBonds = "";
            // for all next neighbors of a specific element
            for (IAtom neighborAtom : ac.getConnectedAtomsList(ac.getAtom(indexAC))) {
                // skip if not the right neighborhood element or bond type is unknown/unset
                if ((!neighborAtom.getSymbol().equals(neighborElems.get(n))) || (casekit.NMR.Utils.getStringFromBondOrder(ac.getBond(ac.getAtom(indexAC), neighborAtom).getOrder()) == null)) {
                    continue;
                }
                foundBonds += casekit.NMR.Utils.getStringFromBondOrder(ac.getBond(ac.getAtom(indexAC), neighborAtom).getOrder());
            }
            for (int k = 0; k < bondsSet.length; k++) {
                counts[n * bondsSet.length + k] = 0;
                if (casekit.NMR.Utils.sortString(foundBonds).equals(casekit.NMR.Utils.sortString(bondsSet[k]))) {
                    counts[n * bondsSet.length + k] = 1;
                    break;
                }
            }
        }
        
        return counts;
    }
        
    /**
     *
     * @param pathToOutput
     * @param m
     * @param bondsSet
     * @param elem
     * @param neighborElems
     * @param min
     * @param max
     * @param stepSize
     * 
     * @throws IOException
     * 
     * @deprecated 
     */
    public static void writeNeighborhoodBondsCountMatrix(final String pathToOutput, final int[][] m, final String[] bondsSet, final String elem, final ArrayList<String> neighborElems, final int min, final int max, final int stepSize) throws IOException{
        
        if(stepSize < 1){
            System.err.println("stepSize < 1 not allowed!!!");
            return;
        }
        final StringBuilder sb = new StringBuilder();
        sb.append("shift [" + elem + "] (" + stepSize + "),nTotal,inRing,isArom,q" + elem + "," + elem + "H," + elem + "H2," + elem + "H3,");
        for (int i = 0; i < neighborElems.size(); i++) {
            for (int j = 0; j < bondsSet.length; j++) {
                sb.append(bondsSet[j] + "[" + neighborElems.get(i) + "]");
                if (j < bondsSet.length - 1) {
                    sb.append(",");
                }
            }
            if (i < neighborElems.size() - 1) {
                sb.append(",");
            }
        }
        sb.append("\n");
        for (int i = 0; i < stepSize * (max - min) + 1; i++) {
            sb.append((i + min) + ",");
            for (int j = 0; j < 3 + 4 + neighborElems.size() * bondsSet.length; j++) {
                sb.append(m[i][j]);
                if (j < 3 + 4 + neighborElems.size() * bondsSet.length - 1) {
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
    
    /**
     *
     * @param s
     * @return
     * 
     * @deprecated 
     */
    public static String sortString(final String s) {
        final char[] c = s.toCharArray();
        Arrays.sort(c);
        return new String(c);
    }
    
    /**
     *
     * @param valences
     * @return
     * 
     * @deprecated 
     */
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
    
    /**
     *
     * @param order
     * @return
     * 
     * @deprecated 
     */
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
    
    
    public static void writeTextFile(final String pathToOutputFile, final String content) throws IOException {
        FileWriter fr = new FileWriter(new File(pathToOutputFile));
        BufferedWriter br = new BufferedWriter(fr);
        br.write(content);
        br.close();
    }
    
    /**
     * Simple function without any settings to generate a picture from a structure 
     * given as IAtomcontainer.
     *
     * @param ac Atom container
     * @param path Path to file for storing
     * @throws IOException
     * @throws CDKException
     */
    public static void generatePicture(final IAtomContainer ac, final String path) throws IOException, CDKException {        
        final DepictionGenerator dg = new DepictionGenerator().withSize(1200, 1200).withAtomColors().withFillToFit().withAtomNumbers();
        dg.depict(ac).writeTo(path);
    }


    /**
     * Detects outliers in given array list of input values and removes them. <br>
     * Here, outliers are those which are outside of a calculated lower and upper bound (whisker).
     * The interquartile range (IQR) of the input values is therefore multiplied with a given value
     * for whisker creation.
     *
     * @param input list of values to process
     * @param multiplierIQR multiplier for IQR to use for lower and upper bound creation
     * @return new array list without values outside the generated boundaries
     */
    public static ArrayList<Double> removeOutliers(final ArrayList<Double> input, final double multiplierIQR){
        final ArrayList<Double> inputWithoutOutliers = new ArrayList<>(input);
        inputWithoutOutliers.removeAll(Utils.getOutliers(inputWithoutOutliers, multiplierIQR));

        return inputWithoutOutliers;
    }

    /**
     *
     * @param input
     * @return
     */
    public static ArrayList<Double> getOutliers(final ArrayList<Double> input, final double multiplierIQR) {
        final ArrayList<Double> outliers = new ArrayList<>();
        if(input.size() <= 1){
            return outliers;
        }
        Collections.sort(input);
        final ArrayList<Double> data1 = new ArrayList<>(input.subList(0, input.size() / 2));
        final ArrayList<Double> data2;
        if (input.size() % 2 == 0) {
            data2 = new ArrayList<>(input.subList(input.size() / 2, input.size()));
        } else {
            data2 = new ArrayList<>(input.subList(input.size() / 2 + 1, input.size()));
        }
        final double q1 = getMedian(new ArrayList<>(data1));
        final double q3 = getMedian(new ArrayList<>(data2));
        final double iqr = q3 - q1;
        final double lowerBound = q1 - multiplierIQR * iqr;
        final double upperBound = q3 + multiplierIQR * iqr;
        for (int i = 0; i < input.size(); i++) {
            if ((input.get(i) < lowerBound) || (input.get(i) > upperBound)) {
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
    public static Double getMedian(final ArrayList<Double> data) {
        if((data == null) || data.isEmpty()) {
            return null;
        }
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
    public static Double getMean(final Collection<Double> data) {
        if((data == null) || data.isEmpty()){
            return null;
        }
        double sum = 0;
        int nullCounter = 0;
        for (final Double d : data) {
            if(d != null){
                sum += d;
            } else {
                nullCounter++;
            }
        }
        return ((data.size() - nullCounter) != 0) ? (sum/(data.size() - nullCounter)) : null;
    }
    
    /**
     *
     * @param data
     * @return
     */
    public static Double getStandardDeviation(final ArrayList<Double> data) {
        if ((data == null) || data.isEmpty()) {
            return null;
        }        
        final Double variance = Utils.getVariance(data);
        
        return (variance != null) ? Math.sqrt(variance) : null;
    }
    
    public static Double getVariance(final Collection<Double> data) {
        if ((data == null) || data.isEmpty()) {
            return null;
        }
        final int nullCounter = Collections.frequency(data, null);
        double quadrSum = 0.0;        
        final Double mean = Utils.getMean(data);
        if(mean == null){
            return null;
        }        
        for (final Double d : data) {
            if (d != null) {
                quadrSum += Math.pow(d - mean, 2);
            }
        }
        
        return ((data.size() - nullCounter) != 0) ? (quadrSum / (data.size() - nullCounter)) : null;        
    }
    
    
    /**
     *
     * @param data
     * @return
     */
    public static Double getMean(final Double[] data) {
        if((data == null) || (data.length == 0)){
            return null;
        }
        double sum = 0;
        int nullCounter = 0;
        for (final Double d : data) {
            if(d != null){
                sum += d;
            } else {
                nullCounter++;
            }
        }
        return ((data.length - nullCounter) != 0) ? (sum / (data.length - nullCounter)) : null;
    }

    public static HashMap<String, Double> getMean(final HashMap<String, ArrayList<Double>> lookup) {

        final HashMap<String, Double> means = new HashMap<>();
        Double meanInList;
        for (final String key : lookup.keySet()) {
            meanInList = Utils.getMean(lookup.get(key));
            if (meanInList != null) {
                means.put(key, meanInList);
            }
        }

        return means;
    }

    public static boolean isValidBondAddition(final IAtomContainer ac, final int atomIndex, final IBond bondToAdd){

        if(ac.getAtom(atomIndex).isAromatic()){
            System.out.println(atomIndex + " --> (" + Utils.getBondOrderSum(ac, atomIndex, true) + " + " + Utils.getBondOrderAsNumeric(bondToAdd) + " - 1) = " + (Utils.getBondOrderSum(ac, atomIndex, true) + Utils.getBondOrderAsNumeric(bondToAdd) - 1) + " <= " + ac.getAtom(atomIndex).getValency() + " ? -> " + ((Utils.getBondOrderSum(ac, atomIndex, true) + Utils.getBondOrderAsNumeric(bondToAdd) - 1) <= ac.getAtom(atomIndex).getValency()));
            return ((Utils.getBondOrderSum(ac, atomIndex, true) - 1) + Utils.getBondOrderAsNumeric(bondToAdd)) <= ac.getAtom(atomIndex).getValency();
        }

        System.out.println(atomIndex + " --> " + Utils.getBondOrderSum(ac, atomIndex, true) + " + " + Utils.getBondOrderAsNumeric(bondToAdd) + " = " + (Utils.getBondOrderSum(ac, atomIndex, true) + Utils.getBondOrderAsNumeric(bondToAdd)) + " <= " + ac.getAtom(atomIndex).getValency() + " ? -> " + ((Utils.getBondOrderSum(ac, atomIndex, true) + Utils.getBondOrderAsNumeric(bondToAdd)) <= ac.getAtom(atomIndex).getValency()));
        return (Utils.getBondOrderSum(ac, atomIndex, true) + Utils.getBondOrderAsNumeric(bondToAdd)) <= ac.getAtom(atomIndex).getValency();
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
     * @deprecated 
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
    
    
    /**
     *
     * @param pathToFile
     * @return
     */
    public static String getFileFormat(final String pathToFile) {

        if(pathToFile == null || pathToFile.trim().isEmpty()){
            return "";
        }        
        final String[] split = pathToFile.split("\\.");

        return split[split.length - 1];
    }
   

    /**
     *
     * @param data
     * @return
     */
    public static Double getRMS(final ArrayList<Double> data) {
        if((data == null) || data.isEmpty()){
            return null;
        }
        if (data.size() == 1) {
            return data.get(0);
        }
        int nullCounter = 0;
        double qSum = 0;
        for (final Double d : data) {
            if(d != null){
                qSum += d * d;
            } else {
                nullCounter++;
            }
        }

        return ((data.size() - nullCounter) != 0) ? Math.sqrt(qSum / (data.size() - nullCounter)) : null;
    }
    
    
    /**
     *
     * @param lookup
     * @return
     */
    public static HashMap<String, Double> getRMS(final HashMap<String, ArrayList<Double>> lookup){
        final HashMap<String, Double> rms = new HashMap<>();
        Double rmsInList;
        for (final String key : lookup.keySet()) {
            rmsInList = Utils.getRMS(lookup.get(key));
            if(rmsInList != null) {
                rms.put(key, rmsInList);
            }
        }
      
        return rms;
    }
    
    public static boolean isSaturated(final IAtomContainer ac, final int atomIndex)  {
        return Utils.getBondOrderSum(ac, atomIndex, true).intValue() >= ac.getAtom(atomIndex).getValency();
    }
    
    public static void addImplicitHydrogens(final IAtomContainer ac) throws CDKException{
        final CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(ac.getBuilder());
        IAtomType type;
        for (IAtom atom : ac.atoms()) {
            type = matcher.findMatchingAtomType(ac, atom);
            AtomTypeManipulator.configure(atom, type);
        }
        final CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(ac.getBuilder());
        adder.addImplicitHydrogens(ac);
    }
    
//    public static int countElements(final String input){
//        int counter = 0;
//        for (int k = 0; k < input.length(); k++) {
//            // Check for uppercase letters
//            if (Character.isLetter(input.charAt(k)) && Character.isUpperCase(input.charAt(k))) {
//                counter++;
//            }
//        }
//
//        return counter;
//    }
    
//    public static ArrayList<String> getComponents(final String symbols){
//        final ArrayList<String> components = new ArrayList<>();
//        for (int i = 0; i < symbols.length(); i++) {
//            if ((i + 1 < symbols.length())
//                    && Character.isLowerCase(symbols.charAt(i + 1))) {
//                components.add(symbols.substring(i, i + 2));
//                i++;
//            } else {
//                components.add(symbols.substring(i, i + 1));
//            }
//        }
//        
//        return components;
//    }
    
    /**
     *
     * @param lookup
     * @return
     * 
     */
    public static HashMap<String, Double> getMedian(final HashMap<String, ArrayList<Double>> lookup) {

        final HashMap<String, Double> medians = new HashMap<>();
        Double medianInList;
        for (final String key : lookup.keySet()) {
            medianInList = Utils.getMedian(lookup.get(key));
            if (medianInList != null) {
                medians.put(key, medianInList);
            }
        }

        return medians;
    }

    /**
     *
     *
     * @param hoseLookupToExtend
     * @param hoseLookup
     *
     * @deprecated
     */
    public static void combineHashMaps(final HashMap<String, ArrayList<Double>> hoseLookupToExtend, final HashMap<String, ArrayList<Double>> hoseLookup){
        for (final String hose : hoseLookup.keySet()) {
            if(!hoseLookupToExtend.containsKey(hose)){
                hoseLookupToExtend.put(hose, new ArrayList<>());
            }
            hoseLookupToExtend.get(hose).addAll(hoseLookup.get(hose));
        }        
    } 
    
    public static Double roundDouble(final Double value, final int decimalPlaces){
        if(value == null){
            return null;
        }
        final int decimalFactor = (int) (Math.pow(10, decimalPlaces));

        return (Math.round(value * decimalFactor) / (double) decimalFactor);
    }
    
    /**
     * Checks whether a structure contains explicit hydrogen atoms or not.
     *
     * @param ac structure to check
     * @return
     */
    public static boolean containsExplicitHydrogens(final IAtomContainer ac){
        for (final IAtom atomA : ac.atoms()) {
            // check each atom whether it is an hydrogen
            if (atomA.getSymbol().equals("H")) {
                return true;
            }
        }
        
        return false;
    }
    
    /**
     * Stores all explicit hydrogens as implicit counter for the bonded heavy 
     * atoms and removes those from the atom container. Also, a HashMap 
     * containing non-hydrogen atoms and its indices 
     * before the removals will be returned which one can use for atom index 
     * comparison (before and after the removals) later. 
     *
     * @param ac the structure to convert
     * @return 
     * 
     * @see #containsExplicitHydrogens(org.openscience.cdk.interfaces.IAtomContainer) 
     */
    public static HashMap<IAtom, Integer> convertExplicitToImplicitHydrogens(final IAtomContainer ac){
        // create a list of atom indices which one can use for index comparison (before vs. after) after removing the explict hydrogens
        final HashMap<IAtom, Integer> atomIndices = new HashMap<>();
        final List<IAtom> toRemoveList = new ArrayList<>();
        IAtom atomB;
        for (final IAtom atomA : ac.atoms()) {
            // check each atom whether it is an hydrogen;
            // if yes then store (increase) the number of implicit hydrogens 
            // for its bonded heavy atom
            if (atomA.getSymbol().equals("H")) {
                atomB = ac.getConnectedAtomsList(atomA).get(0);
                if(atomB.getImplicitHydrogenCount() == null){
                    atomB.setImplicitHydrogenCount(0);
                }
                atomB.setImplicitHydrogenCount(atomB.getImplicitHydrogenCount() + 1);
                toRemoveList.add(atomA);                                                
            } else {
                // store all non-hydrogen atoms and their indices
                atomIndices.put(atomA, atomA.getIndex());
            }
            
        }
        // remove all explicit hydrogen atoms
        for (final IAtom iAtom : toRemoveList) {
            ac.removeAtom(iAtom);
        }
        
        return atomIndices;
    }
    
    /**
     *
     * @param ac
     * @return 
     */
    public static int getExplicitHydrogenCount(final IAtomContainer ac){
        final List<IAtom> toRemoveList = new ArrayList<>();
        IAtom atomB;
        for (final IAtom atomA : ac.atoms()) {
            if (atomA.getAtomicNumber() == 1) {
                atomB = ac.getConnectedAtomsList(atomA).get(0);
                if(atomB.getImplicitHydrogenCount() == null){
                    atomB.setImplicitHydrogenCount(0);
                }
                atomB.setImplicitHydrogenCount(atomB.getImplicitHydrogenCount() + 1);
                toRemoveList.add(atomA);
            }
        }
        
        return toRemoveList.size();
    }
    
    
    public static void setAromaticity(final IAtomContainer ac) throws CDKException {
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(ac);
        final ElectronDonation model = ElectronDonation.cdkAllowingExocyclic();
        final CycleFinder cycles = Cycles.all(ac.getAtomCount());
        final Aromaticity aromaticity = new Aromaticity(model, cycles);
        aromaticity.apply(ac);
    }

    public static void setAromaticityAndKekulize(final IAtomContainer ac) throws CDKException {
        Utils.setAromaticity(ac);
        Kekulization.kekulize(ac);
    }
    
    
    /**
     * Removes atoms from a given atom type from an atom container.
     *
     * @param ac IAtomContainer object where to remove the atoms
     * @param atomType Atom type (element's name, e.g. C or Br)
     * @return IAtomContainer where the atoms were removed
     */
    public static IAtomContainer removeAtoms(final IAtomContainer ac, final String atomType) {

        final ArrayList<IAtom> toRemoveList = new ArrayList<>();
        for (IAtom atomA : ac.atoms()) {
            if (atomA.getSymbol().equals(atomType)) {// detect whether the current atom A is a from the given atom type
                toRemoveList.add(atomA);
            }
        }
        for (IAtom iAtom : toRemoveList) {
            ac.removeAtom(iAtom);
        }

        return ac;
    }
    
    /**
     *
     * @param array
     * @return
     * 
     */
    public static ArrayList<Integer> ArrayToArrayList(final int[] array){
        
        final ArrayList<Integer> list = new ArrayList<>();
        for (int i = 0; i < array.length; i++) {
            list.add(array[i]);
        }
        
        return list;
    }
    
    
    public static String getSpectrumNucleiAsString(final Spectrum spectrum){
        String specID = "";
        for (int i = 0; i < spectrum.getNDim(); i++) {
            specID += spectrum.getNuclei()[i];
            if(i < spectrum.getNDim() - 1){
                specID += "-";
            }
        }
        
        return specID;
    }
    
    public static boolean checkIndexInAtomContainer(final IAtomContainer ac, final int atomIndex){
        return ((atomIndex >= 0) && atomIndex < ac.getAtomCount());
    }    
    
    public static ExecutorService initExecuter(final int nThreads) {
        return Executors.newFixedThreadPool(nThreads);
    }

    public static void stopExecuter(final ExecutorService executor, final long seconds) {      
        executor.shutdown();
        try {
            if (!executor.awaitTermination(seconds, TimeUnit.SECONDS)) {
                System.err.println("killing non-finished tasks!");
                executor.shutdownNow();
            }
        } catch (InterruptedException e) {
            System.err.println("killing non-finished tasks!");
            executor.shutdownNow();
        }
    }        
    
     /**
     * Returns the bond order for a numeric order value.
     *
     * @param orderAsNumeric 
     * @return
     */
    public static IBond.Order getBondOrder(final int orderAsNumeric) {
        for (IBond.Order order : IBond.Order.values()){
            if(order.numeric() == orderAsNumeric){
                return order;
            }
        }            
        
        return null;
    }

    public static Float getBondOrderAsNumeric(final IBond bond) {
        if(bond == null){
            return null;
        }
        float bondOrderAsNumeric;
        if (bond.isAromatic()) {
            bondOrderAsNumeric = (float) 1.5;
        } else {
            bondOrderAsNumeric = bond.getOrder().numeric();
        }
        
        return bondOrderAsNumeric;
    }
    
    public static Float getBondOrderSum(final IAtomContainer ac, final int atomIndex, final boolean includeImplicitHydrogenCount) {        
        if(!Utils.checkIndexInAtomContainer(ac, atomIndex)){
            return null;
        }
        float bondsOrderSum = 0;
        final IAtom atom = ac.getAtom(atomIndex);
        for (final IBond bond : ac.getConnectedBondsList(atom)) {
            bondsOrderSum += Utils.getBondOrderAsNumeric(bond);
        }
        if(includeImplicitHydrogenCount && (atom.getImplicitHydrogenCount() != null)){
            bondsOrderSum += atom.getImplicitHydrogenCount();
        }
        
        return bondsOrderSum;
    }
    
}
