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
package casekit.NMR;


import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import org.apache.commons.lang3.ArrayUtils;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.Fingerprinter;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.qsar.descriptors.atomic.AtomHybridizationDescriptor;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.similarity.Tanimoto;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class ShiftMatcher {
    
    public ShiftMatcher(){
        this.HOSESymbolHash = new Hashtable();
        this.verbose = true;
        
//        this.HOSESymbolHash.put("H", "H");
//        this.HOSESymbolHash.put("C", "C");
        this.HOSESymbolHash.put("O", "O");
        this.HOSESymbolHash.put("N", "N");
        this.HOSESymbolHash.put("S", "S");
        this.HOSESymbolHash.put("P", "P");
        this.HOSESymbolHash.put("Si", "Q");
        this.HOSESymbolHash.put("B", "B");
        this.HOSESymbolHash.put("F", "F");
        this.HOSESymbolHash.put("Cl", "X");
        this.HOSESymbolHash.put("Br", "Y");
        this.HOSESymbolHash.put("I", "I");
        
//        this.HOSESymbolHash.put("Se", "Se");
//        this.HOSESymbolHash.put("Sn", "Sn");
//        this.HOSESymbolHash.put("Ge", "Ge");
//        this.HOSESymbolHash.put("te", "Te");
//        this.HOSESymbolHash.put("Zn", "Zn");
//        this.HOSESymbolHash.put("As", "As");
//        this.HOSESymbolHash.put("Li", "Li");
//        this.HOSESymbolHash.put("Ti", "Ti");
//        this.HOSESymbolHash.put("Pb", "Pb");
//        this.HOSESymbolHash.put("Hg", "Hg");
//        this.HOSESymbolHash.put("Mg", "Mg");
        
    }

    private final boolean verbose;
    private Hashtable<String,ArrayList<Double>> hoseLookup;
    private int carbonNMRCount;
    private int hydrogenNMRCount;
    private int molListDBCount;
    private final Hashtable<String, String> HOSESymbolHash;
    private ArrayList<IAtomContainer> molListDB;
    
    
    
    
    
    private ArrayList<String> convertToHOSECodeTerm(final String molFormula){
        
        final ArrayList<String> convertedList = new ArrayList<>();
        Matcher m = Pattern.compile("[A-Z]").matcher(molFormula);
        String s;
        
        while (m.find()) {  
            s = molFormula.substring(m.start(), m.end());
            if(m.end() + 1 <= molFormula.length()){
                Character c = molFormula.subSequence(m.end(), m.end() + 1).charAt(0);
                if(Character.isLowerCase(c)){
                    s = molFormula.substring(m.start(), m.end() + 1);
                } 
            } 
            // ignore C and H atoms
            if(s.equals("C") || s.equals("H")) continue;
            
            if(this.HOSESymbolHash.containsKey(s)) {
                convertedList.add(this.HOSESymbolHash.get(s));
            } else {
                convertedList.add(s);
            }
        }

        return convertedList;
    }
    
    
    public void matchHOSEShifts(final double[] nmrShiftValues, final String molFormula) throws Exception {
        if (verbose) {
            System.out.println("Starting shift matching (2) against hose code table");
        }
        
        Arrays.parallelSort(nmrShiftValues);
        
        final Set<String> keys = this.hoseLookup.keySet();
        ArrayList<Double> shifts;
        
        final Hashtable<String, Integer> hoseHeteroCount = new Hashtable<>();
        final Hashtable<String, Integer> hoseHCount = new Hashtable<>();
        
        StringTokenizer strtok;
        StringTokenizer strtok2;

        final ArrayList<String> atomTypesHetero = convertToHOSECodeTerm(molFormula);
        System.err.println("found heavy atoms in molecular formula: " + molFormula + " -> " + atomTypesHetero.toString());
        
        
        Matcher m;
        int heteroAtomCount;
        int carbonAtomCount;
        int HAtomCount;
        int totalBindPartnCount;
        int allAtomCount;

        int hitCounter;
//        int middle;
//        double mean;
//        double median;
//        double d1;
//        double d2;
//        double sum;
//        double sd;
        
        String neighbors;

        
        for (int i = 0; i < nmrShiftValues.length; i++) {
            final ArrayList<String> hoseMatch = new ArrayList<>();
        
            hitCounter = 0;

            for (String hose : keys) {
                
                strtok = new StringTokenizer(hose, ";");
                totalBindPartnCount = Integer.parseInt(strtok.nextToken().substring(2, 3));
                strtok2 = new StringTokenizer(strtok.nextToken(), "(");

                String directNeighbors = strtok2.nextToken();

                heteroAtomCount = 0;
                for (String s : atomTypesHetero) {
                    if (s.length() == 1) {
                        m = Pattern.compile(s +"(?![a-z])").matcher(directNeighbors);
                    } else {
                        m = Pattern.compile(s).matcher(directNeighbors);
                    }
                    while (m.find()) {
                        heteroAtomCount++;
                    }
                }
                
                m = Pattern.compile("C(?![a-z])").matcher(directNeighbors);
                carbonAtomCount = 0;
                while (m.find()) {
                    carbonAtomCount++;
                }
                
                //check (count) whether other atom types than C or needed hetero atoms are present or not
                m = Pattern.compile("[A-Z]").matcher(directNeighbors);
                allAtomCount = 0;
                while (m.find()) {
                    allAtomCount++;
                }                
                if(allAtomCount - (carbonAtomCount + heteroAtomCount) != 0){
//                    System.out.println(hose + " -> " + allAtomCount + " - (" + carbonAtomCount + " + " + heteroAtomCount + ")");
//                    System.out.println(" -> REFUSED!!!");
                    continue;
                }
                
                HAtomCount = totalBindPartnCount - (carbonAtomCount + heteroAtomCount);
                
                hoseHeteroCount.put(hose, heteroAtomCount);
                hoseHCount.put(hose, HAtomCount);

                
                shifts = hoseLookup.get(hose);
                if ((nmrShiftValues[i] >= shifts.get(shifts.indexOf(Collections.min(shifts)))) && nmrShiftValues[i] <= shifts.get(shifts.indexOf(Collections.max(shifts)))) {

                    hitCounter++;
                    hoseMatch.add(hose);
                }
            }
            System.out.println("\nThere are " + hitCounter + " hits for shift value " + nmrShiftValues[i]);
            
            for (int j = 4; j >= 0; j--) {
            
                System.out.println(j + ":");
                final ArrayList<String> directNeighbors = new ArrayList<>();
                
                for (String hose : hoseMatch) {
                    if (hoseHeteroCount.get(hose) == j) {
//                        System.out.println(hose + ": " + hoseHeteroAtoms.get(hose));
//                        shifts = hoseLookup.get(hose);
//                        
//                        middle = shifts.size() / 2;
//                        if (shifts.size() % 2 == 1) {
//                            median = shifts.get(middle);
//                        } else {
//                            median = (shifts.get(middle - 1) + shifts.get(middle)) / 2.0;
//                        }
//
//                        sum = 0;
//                        for (int k = 0; k < shifts.size(); k++) {
//                            sum += shifts.get(k);
//                        }
//                        mean = (sum / shifts.size());
//
//                        d1 = 0;
//                        d2 = 0;
//                        sum = 0;
//                        for (int k = 0; k < shifts.size(); k++) {
//                            d2 = (mean - shifts.get(k)) * (mean - shifts.get(k));
//                            d1 = d2 + d1;
//                        }
//                        sd = 0;
//                        if (shifts.size() > 1) {
//                            sd = Math.sqrt((d1 / (shifts.size() - 1)));
//                        }
//
//                        System.out.println(hose + " [" + shifts.get(shifts.indexOf(Collections.min(shifts))) + ", " + shifts.get(shifts.indexOf(Collections.max(shifts))) + "] (" + shifts.size() + ") -> mean: " + mean + ", median: " + median + ", sd: " + sd + "\n");

                        strtok = new StringTokenizer(hose, ";");
                        strtok.nextToken();
                        strtok2 = new StringTokenizer(strtok.nextToken(), "(");
                        neighbors = strtok2.nextToken();
                        
                        neighbors = neighbors + "\t\t -> het: " + j + ", H: " + hoseHCount.get(hose);
                        
                        directNeighbors.add(neighbors);
                    }
                }
                
                //unique HOSE codes
                HashSet hs = new HashSet();
                hs.addAll(directNeighbors);
                directNeighbors.clear();
                directNeighbors.addAll(hs);

                for (String n : directNeighbors) {
                    System.out.println("\t" + n);
                }
                
            }
        }
//        System.out.println("\n");

        
    }
 
    
    public void matchDBBits(final IAtomContainer acQ, final String elementSymbol) throws CDKException{
        
        final double[] nmrShiftValuesQ = new double[acQ.getAtomCount()];
        for (int i = 0; i < acQ.getAtomCount(); i++) {
            nmrShiftValuesQ[i] = acQ.getAtom(i).getProperty(CDKConstants.NMRSHIFT_CARBON);
        }
        System.out.println("Q shifts: " + Arrays.toString(nmrShiftValuesQ));
        
        Arrays.parallelSort(nmrShiftValuesQ);
        
        ArrayList<Double>  nmrShiftValuesDBList;
        double[] nmrShiftValuesDB;
        Double[] temp;
        for(int i = 0; i< this.molListDBCount; i++){
            IAtomContainer acDB = this.molListDB.get(i);
            System.out.println("acDB atom count: " + acDB.getAtomCount());
            nmrShiftValuesDBList = new ArrayList<>();
            
            Fingerprinter fp = new Fingerprinter();
            IBitFingerprint ifpDB = fp.getBitFingerprint(acDB);
            IBitFingerprint ifpQ = fp.getBitFingerprint(acQ);
            System.out.println(Arrays.toString(ifpQ.getSetbits()));
            System.out.println(Arrays.toString(ifpDB.getSetbits()) + "\n");

            for(int j = 0; j < acDB.getAtomCount(); j++) {
                if(acDB.getAtom(j).getSymbol().equals(elementSymbol)){
                    nmrShiftValuesDBList.add(acDB.getAtom(j).getProperty(CDKConstants.NMRSHIFT_CARBON));
                }
            }
            temp = nmrShiftValuesDBList.toArray(new Double[nmrShiftValuesDBList.size()]);
            nmrShiftValuesDB = ArrayUtils.toPrimitive(temp);
            Arrays.parallelSort(nmrShiftValuesDB);

            System.out.println("DB shifts: " + Arrays.toString(nmrShiftValuesDB));
            
            double tanimoto_coefficient = Tanimoto.calculate(nmrShiftValuesQ, nmrShiftValuesDB);
            System.out.println(i + ": Tanimo result: " + tanimoto_coefficient);
        }
        
    }
    
    
    
    
    

}
