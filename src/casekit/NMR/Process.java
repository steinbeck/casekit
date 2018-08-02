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

import casekit.NMR.model.Spectrum;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;
import org.openscience.cdk.Atom;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.qsar.descriptors.atomic.AtomHybridizationDescriptor;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Process extends ParseRawData {
    
    final private IAtomContainer mol;
    final private IMolecularFormula molFormula;
    private HashMap<String, ArrayList<Integer>> atomTypeIndices = new HashMap<>();    
    private int[][] neighborhoodCountsMatrix;
    final private HashMap<Integer, ArrayList<Integer[]>> shiftIndicesInACSet = new HashMap<>(); // holding of all indices of each ac set (DB) entry [first value] and it's atom indices [second value] too
    
    
    public Process(){
        super();
        this.molFormula = super.getMolecularFormula();
        this.mol = super.getAtomContainer();
    }
    
    public Process(final IMolecularFormula molFormula){
        super(molFormula);
        this.molFormula = molFormula;
        this.mol = Utils.removeAtoms(MolecularFormulaManipulator.getAtomContainer(this.molFormula), "H");
        this.setAtomTypeIndices();
    }

    
    
    
    
//    /**
//     * Sets the hybridization level of each heavy atom in the molecule which has
//     * its own shift value (property), only if a frequency threshold value for one 
//     * hybridization level is reached.
//     * For further details see {@link testkit.Utils#getHybridizationsFromNMRShiftDB(IAtomContainer, String, double, IMolecularFormula)}
//     * Two threshold value are used to accept a hybridization level and a found 
//     * neighbor as real neighbor (see thrs parameter descriptions).
//     *
//     * @param pathToNMRShiftDB path to NMRShiftDB sdf file
//     * @param tol tolerance value [ppm] for atom shift matching in DB
//     * @param thrsHybrid threshold for accepting a hybridization frequency rate, e.g.
//     * the value 0.9 means that 90% of all found hybridizations for given carbon
//     * shift must be from the same hybridization level
//     * @param thrsNeighbor threshold for accepting a found neighbor frequency rate 
//     * (atom type) in database as real neighbor for the unknown structure. E.g.
//     * the value 0.9 means that 90% of all found neighbors for given carbon
//     * shift must be from the same atom type, like nitrogen or oxygen.
//     *   
//     * @throws IOException
//     */
//    public void setHybridizationsFromNMRShiftDB(final String pathToNMRShiftDB, final double tol, final double thrsHybrid, final double thrsNeighbor) throws IOException {
//
//        final HashMap<String, HashMap<Integer, HashMap<String, ArrayList<Integer>>>> elementsHybridAndBondTypeCounter = testkit.Utils.getHybridizationsFromNMRShiftDB(this.mol, pathToNMRShiftDB, tol, this.molFormula);
//        final HashMap<Integer, HashMap<String, ArrayList<Integer>>> elementsHybridCounter = elementsHybridAndBondTypeCounter.get("hybridCounter");
//        final HashMap<Integer, HashMap<String, ArrayList<Integer>>> elementsBondTypeCounter = elementsHybridAndBondTypeCounter.get("bondTypeCounter");
//
//        ArrayList<Integer> hybrids;
//        final HashMap<String, Double[]> probsNeighbors = new HashMap<>();
//        int idx = 0;
//        Double[] temp;
//        final HashMap<String, Integer> elementsHybridCounterSum = new HashMap<>();
//        for (int qAtomIndex : elementsHybridCounter.keySet()) {
//            for (String keyValue : elementsHybridCounter.get(qAtomIndex).keySet()) {
//                if(keyValue.equals("query") || keyValue.equals("queryH")){
//                    continue;
//                }
//                if(!probsNeighbors.containsKey(keyValue)){
//                    probsNeighbors.put(keyValue, new Double[elementsHybridCounter.keySet().size()]);
//                }
//                hybrids = elementsHybridCounter.get(qAtomIndex).get(keyValue);
//                temp = probsNeighbors.get(keyValue);
//                temp[idx] = (double) hybrids.size();
//                probsNeighbors.put(keyValue, temp);
//                
//                if(!elementsHybridCounterSum.containsKey(keyValue)){
//                    elementsHybridCounterSum.put(keyValue, 0);
//                }
//                elementsHybridCounterSum.put(keyValue, elementsHybridCounterSum.get(keyValue) + hybrids.size());
//            }
//            idx++;
//        }
//        int sumQueryAtom;
//        for (int i = 0; i < idx; i++) {
//            sumQueryAtom = 0;
//            for (String keyValue : probsNeighbors.keySet()) {
//                sumQueryAtom += probsNeighbors.get(keyValue)[i];
//            }
//            for (String keyValue : probsNeighbors.keySet()) {
//                temp = probsNeighbors.get(keyValue);
//                temp[i] = temp[i]/sumQueryAtom;//0.5 * (temp[i]/sumQueryAtom + temp[i]/elementsHybridCounterSum.get(keyValue));
//                probsNeighbors.put(keyValue, temp);
//            }
//        }
//        
//      
//        HashMap<Integer, Double> hybridFreqs;
//        int maxFreqHybridValue;
//        double maxFreq;
//        IAtom qAtom;
//        // for all query atoms which have their own NMR shift value
//        idx = 0;
//        for (int qAtomIndex : elementsHybridCounter.keySet()) {
//            qAtom = this.mol.getAtom(qAtomIndex);
//            System.out.println("\nmain key: " + qAtomIndex +  " -> H: " + qAtom.getImplicitHydrogenCount() + ", nmr shift: " + qAtom.getProperty(testkit.Utils.getNMRShiftConstant(qAtom.getSymbol())));
//            // for all possible neighbors
//            for (String keyValue : elementsHybridCounter.get(qAtomIndex).keySet()) {
//                hybrids = elementsHybridCounter.get(qAtomIndex).get(keyValue);
//                if(hybrids.isEmpty()){
//                    continue;
//                }
//                hybridFreqs = testkit.Utils.getValueFrequencies(hybrids);
//                maxFreqHybridValue = -1;
//                maxFreq = Collections.max(hybridFreqs.values());
//                for (int hybridValue : hybridFreqs.keySet()) {
//                    if(hybridFreqs.get(hybridValue) == maxFreq){
//                        maxFreqHybridValue = hybridValue;
//                        break;
//                    }
//                }
//                
//                // set hybridization for a query atom which has at least one match with an attached hydrogen shift
//                // value for a matched heavy atom in DB; this method is preferred
//                switch (keyValue) {
//                    case "queryH":
//                        System.out.println("queryH -> " + hybrids.size() + " -> " + IAtomType.Hybridization.values()[maxFreqHybridValue] + " (" + hybridFreqs.get(maxFreqHybridValue) + ")");
//                        if(hybridFreqs.get(maxFreqHybridValue) >= thrsHybrid){
//                            qAtom.setHybridization(IAtomType.Hybridization.values()[maxFreqHybridValue]);
////                        System.out.println("queryH -> " + hybrids.size() + " -> " + qAtom.getHybridization() + " (" + hybridFreqs.get(maxFreqHybridValue) + ")");
//                        }   break;
//                    case "query":
//                        System.out.println("query -> " + hybrids.size() + " -> " + IAtomType.Hybridization.values()[maxFreqHybridValue] + " (" + hybridFreqs.get(maxFreqHybridValue) + ")");
//                        if(qAtom.getHybridization() == null && (hybridFreqs.get(maxFreqHybridValue) >= thrsHybrid)){
//                            // set hybridization from DB entries without attached hydrogen shift matches for an heavy atom
//                            qAtom.setHybridization(IAtomType.Hybridization.values()[maxFreqHybridValue]);
////                        System.out.println("query -> " + hybrids.size() + " -> " + qAtom.getHybridization() + " (" + hybridFreqs.get(maxFreqHybridValue) + ")");
//                        }   break;
//                    default:
//                        System.out.println(idx + ": " + keyValue + ": " + probsNeighbors.get(keyValue)[idx] + " (" + IAtomType.Hybridization.values()[maxFreqHybridValue] + ", " + hybridFreqs.get(maxFreqHybridValue) + ")");
//                        HashMap<Integer, Double> freqs = testkit.Utils.getValueFrequencies(elementsBondTypeCounter.get(qAtomIndex).get(keyValue));
//                        for (Integer bondType : freqs.keySet()) {
//                            System.out.println(" -> " + IBond.Order.values()[bondType - 1] + " (" + freqs.get(bondType) + ")");
//                        }   break;
//                }
//            }
//            idx++;
//        }
//
//        System.out.println("\n");
//        idx = 0;
//        for (int qAtomIndex : elementsHybridCounter.keySet()) {
//            qAtom = this.mol.getAtom(qAtomIndex);
//            String output = qAtomIndex + "\t(" + String.format( "%.3f", (double) qAtom.getProperty(testkit.Utils.getNMRShiftConstant(qAtom.getSymbol()))) + ",\t" + qAtom.getHybridization() + ",\tH:" + qAtom.getImplicitHydrogenCount() + "):\t";
//            for (String keyValue : elementsHybridCounter.get(qAtomIndex).keySet()) {
//                if (keyValue.equals("queryH") || keyValue.equals("query")) {
//                    continue;
//                }
//                if(probsNeighbors.get(keyValue)[idx] >= 0.1){
//                    output += keyValue + ": " + String.format( "%.3f", probsNeighbors.get(keyValue)[idx]) + " ";
//                } else {
//                    output += keyValue + ": ----- ";
//                }   
//            }  
//            
////            for (IAtom neighbor : this.mol.getConnectedAtomsList(qAtom)) {
////                if(neighbor.getProperty(testkit.Utils.getNMRShiftConstant(neighbor.getSymbol())) == null){
////                    output += " -> " + neighbor.getSymbol();
////                }
////            }
//            
//            
//            System.out.println(output);
//            idx++;
//            
//            testkit.Utils.getOpenBonds(this.mol, qAtomIndex);
//        }
//    }
   
    
    /**
     * Sets bonds from already set experiment information (H,H-COSY, INADEQUATE and HMBC).
     * Additionally, this function is build for bond type recognition, 
     * for details see {@link testkit.Utils#getBondTypeFromHybridizations(java.lang.String, org.openscience.cdk.interfaces.IAtomType.Hybridization, java.lang.String, org.openscience.cdk.interfaces.IAtomType.Hybridization)}.
     * 
     */
    public void setBonds(final String[] experiments){
        
        String NMRSHIFT_ATOMTYPE;
        ArrayList<Integer> signalList;
        for (int e = 0; e < experiments.length; e++) {
            for (int i = 0; i < this.mol.getAtomCount(); i++) {
                NMRSHIFT_ATOMTYPE = casekit.NMR.Utils.getNMRShiftConstant(this.mol.getAtom(i).getSymbol());
                // is the NMR shift constant defined and does the nmr shift property entry in an atom exist?
                if (NMRSHIFT_ATOMTYPE != null && this.mol.getAtom(i).getProperty(NMRSHIFT_ATOMTYPE) != null) {
                    if (this.mol.getAtom(i).getProperties().containsKey(experiments[e])) {
                        signalList = this.mol.getAtom(i).getProperty(experiments[e]);
                        for (int bondPartnerIndex : signalList) {
                            // no bonds on one and the same atom; ignore already set bonds if no override wanted
                            if ((i == bondPartnerIndex)) {// || (this.mol.getBond(this.mol.getAtom(i), this.mol.getAtom(bondPartnerIndex)) != null)) {
                                continue;
                            }
                            if(experiments[e].equals(CDKConstants.NMRSPECTYPE_2D_HMBC)){
                                System.out.println("HMBC bond setting: still to come!!!");
                            } else {
                                this.setBond(i, bondPartnerIndex);
                            }
                        }
                    }
                }
            }
        }
    }
    
    
    private void setBond(final int index1, final int index2){
        
        if(this.mol.getBond(this.mol.getAtom(index1), this.mol.getAtom(index2)) != null){
                this.mol.removeBond(this.mol.getAtom(index1), this.mol.getAtom(index2));
        }
        this.mol.addBond(index1, index2, casekit.NMR.Utils.getBondTypeFromHybridizations(this.mol.getAtom(index1), this.mol.getAtom(index2)));
    }
    
    
    /**
     * Adds a bond manually after reading the experimental data and
     * setting bonds from that automatically.
     *
     * @param atomType1 Element name (e.g. "C") for the first heavy atom
     * which also occurrs in
     * {@link testkit.Utils#getNMRShiftConstant(java.lang.String)}
     * @param shift1 shift value [ppm] of the first heavy atom for matching
     * @param tol1 tolerance value for first heavy atom matching
     * @param atomType2 Element name (e.g. "C") for the second heavy atom
     * which also occurrs in
     * {@link testkit.Utils#getNMRShiftConstant(java.lang.String)}
     * @param shift2 shift value [ppm] of the second heavy atom for matching
     * @param tol2 tolerance value for second heavy atom matching
     * @return returns false if no matches were found and no bond could
     * be stored or the matched atom indices are the same, otherwise true
     * @deprecated
     */
    public boolean addBond(final String atomType1, final double shift1, final double tol1, final String atomType2, final double shift2, final double tol2) {

        final String NMRSHIFT_ATOMTYPE1 = casekit.NMR.Utils.getNMRShiftConstant(atomType1);
        final String NMRSHIFT_ATOMTYPE2 = casekit.NMR.Utils.getNMRShiftConstant(atomType2);
        if ((NMRSHIFT_ATOMTYPE1 == null) || (NMRSHIFT_ATOMTYPE2 == null)) {
            return false;
        }
        int atomIndex1 = casekit.NMR.Utils.findSingleShiftMatch(this.mol, shift1, tol1, atomType1);
        int atomIndex2 = casekit.NMR.Utils.findSingleShiftMatch(this.mol, shift2, tol2, atomType2);
        if ((atomIndex1 < 0) || (atomIndex2 < 0) || (atomIndex1 == atomIndex2)) {
            return false;
        }
        this.setBond(atomIndex1, atomIndex2);

        
        return true;
    }
    
    
    /**
     * Adds a H,H-COSY signal and bond between two heavy atoms. To add such a signal, 
     * at least the first heavy atom has to have a shift value match within 
     * the atom container of the unknown.
     * For matching the second heavy atom and creating a (pseudo) HSQC signal,
     * the function {@link #addHSQC(java.lang.String, double, double, double) } 
     * is used.
     *
     * @param atomType1 Element name (e.g. "C") for the first heavy atom 
     * which also occurrs in {@link testkit.Utils#getNMRShiftConstant(java.lang.String)}
     * @param shift1 shift value [ppm] of the first heavy atom for matching 
     * @param tol1 tolerance value for first heavy atom matching 
     * @param atomType2 Element name (e.g. "C") for the second heavy atom
     * which also occurrs in {@link testkit.Utils#getNMRShiftConstant(java.lang.String)}
     * @param shift2 shift value [ppm] of the second heavy atom for matching 
     * @param tol2  tolerance value for second heavy atom matching 
     * @param shiftH proton shift value [ppm] to store
     * @return returns false if no matches were found and no COSY signal could 
     * be stored or the matched atom indices are the same, otherwise true
     * @deprecated
     */
    public boolean addCOSY(final String atomType1, final double shift1, final double tol1, final String atomType2, final Double shift2, final double tol2, final Double shiftH){
        
        final String NMRSHIFT_ATOMTYPE1 = casekit.NMR.Utils.getNMRShiftConstant(atomType1);
        final String NMRSHIFT_ATOMTYPE2 = casekit.NMR.Utils.getNMRShiftConstant(atomType2);
        if ((NMRSHIFT_ATOMTYPE1 == null) || (NMRSHIFT_ATOMTYPE2 == null)) {
            return false;
        }
        int atomIndex1 = casekit.NMR.Utils.findSingleShiftMatch(this.mol, shift1, tol1, atomType1);
        int atomIndex2 = this.addHSQC(atomType2, shift2, tol2, shiftH);
        if ((atomIndex1 < 0) || (atomIndex2 < 0) || (atomIndex1 == atomIndex2)) {
            return false;
        }
        if(this.mol.getAtom(atomIndex1).getProperty(CDKConstants.NMRSPECTYPE_2D_HHCOSY) == null){
            this.mol.getAtom(atomIndex1).setProperty(CDKConstants.NMRSPECTYPE_2D_HHCOSY, new ArrayList<>());
        }
        if(this.mol.getAtom(atomIndex2).getProperty(CDKConstants.NMRSPECTYPE_2D_HHCOSY) == null){
            this.mol.getAtom(atomIndex2).setProperty(CDKConstants.NMRSPECTYPE_2D_HHCOSY, new ArrayList<>());
        }
        
        final ArrayList<Integer> COSYList = this.mol.getAtom(atomIndex1).getProperty(CDKConstants.NMRSPECTYPE_2D_HHCOSY);
        final ArrayList<Integer> COSYListX = this.mol.getAtom(atomIndex2).getProperty(CDKConstants.NMRSPECTYPE_2D_HHCOSY);
        COSYList.add(atomIndex2);
        COSYListX.add(atomIndex1);
        
        this.setBond(atomIndex1, atomIndex2);
        
        // set new hybridization of the COSY partner 
        final AtomHybridizationDescriptor desc = new AtomHybridizationDescriptor();
        this.mol.getAtom(atomIndex1).setHybridization(IAtomType.Hybridization.values()[Integer.parseInt(desc.calculate(this.mol.getAtom(atomIndex1), this.mol).getValue().toString())]);
        this.mol.getAtom(atomIndex2).setHybridization(IAtomType.Hybridization.values()[Integer.parseInt(desc.calculate(this.mol.getAtom(atomIndex2), this.mol).getValue().toString())]);
        
        return true;
    }
    
    
    /**
     * Adds a HSQC signal manually after reading the experimental data and setting bonds from that automatically.
     * If a shift value for a heavy atom is >0.0 then this shift value will be used
     * to find a heavy atom match between this given shift value and atoms 
     * of the atom container of the unknown. Otherwise the first heavy atom without stored
     * NMR shift entry and without stored proton shifts in the atom container 
     * is used for attaching a proton. Additionally, a given proton shift value >0.0 
     * is used to store it into the matched heavy atom's proton shift list.
     * 
     * @param atomType atom type used for matching
     * @param shift shift valuen [ppm] of the heavy atom 
     * @param tol tolerance value [ppm] for matching
     * @param shiftH proton shift value [ppm] to store
     * @return index of matched heavy atom within the atom container; returns -1 if no heavy atom match was found
     * @deprecated 
     */
    public int addHSQC(final String atomType, final Double shift, final double tol, final Double shiftH ){
        
        int atomIndex = -1;
        final String NMRSHIFT_ATOMTYPE = casekit.NMR.Utils.getNMRShiftConstant(atomType);
        if ((NMRSHIFT_ATOMTYPE == null) || (this.atomTypeIndices.get(atomType) == null)) {
            return -1;
        }
        // set additional HSQC for an atom with already set shift value
        if(shift != null){
            atomIndex = casekit.NMR.Utils.findSingleShiftMatch(this.mol, shift, tol, atomType);
        } else {
            // set HSQC for the first atom of given atom type without a already set shift value and without attached proton shifts
            for (Integer i : this.atomTypeIndices.get(atomType)) {
                if ((this.mol.getAtom(i).getProperty(NMRSHIFT_ATOMTYPE) == null) && (this.mol.getAtom(i).getProperty(CDKConstants.NMRSPECTYPE_2D_HSQC) == null)) {
                    atomIndex = i;
                    break;
                }
            }
        }
        // if no atom found to attach a proton
        if (atomIndex < 0) {
            return -1;
        }
        // add the proton shift value if it is higher than 0
        if(shiftH != null){
            if (this.mol.getAtom(atomIndex).getProperty(CDKConstants.NMRSPECTYPE_2D_HSQC) == null) {
                this.mol.getAtom(atomIndex).setProperty(CDKConstants.NMRSPECTYPE_2D_HSQC, new ArrayList<>());
            }
            final ArrayList<Double> protonShifts = this.mol.getAtom(atomIndex).getProperty(CDKConstants.NMRSPECTYPE_2D_HSQC);
            protonShifts.add(shiftH);
        }
        // increase the implicit proton number
        if(this.mol.getAtom(atomIndex).getImplicitHydrogenCount() == null){
            this.mol.getAtom(atomIndex).setImplicitHydrogenCount(0);
        }
        this.mol.getAtom(atomIndex).setImplicitHydrogenCount(this.mol.getAtom(atomIndex).getImplicitHydrogenCount() + 1);
        // set the (new) hybridization
        final AtomHybridizationDescriptor desc = new AtomHybridizationDescriptor();
        this.mol.getAtom(atomIndex).setHybridization(IAtomType.Hybridization.values()[Integer.parseInt(desc.calculate(this.mol.getAtom(atomIndex), this.mol).getValue().toString())]);

        
        return atomIndex;
    }
    
    /**
     *
     * @param atomType
     * @param shift
     * @deprecated 
     */
    public void addAtom(final String atomType, final Double shift){
       
        this.mol.addAtom(new Atom(atomType));
        if(shift != null){
            this.mol.getAtom(this.mol.getAtomCount() - 1).setProperty(casekit.NMR.Utils.getNMRShiftConstant(atomType), shift);
        }
        this.setAtomTypeIndices();
    }
    
    
    /**
     *
     * @param projectName
     * @param pathToOutputFile
     * @param pathsToFilters
     * @throws FileNotFoundException
     * @throws UnsupportedEncodingException
     */
    public void createLSDFile(final String projectName, final String pathToOutputFile, final String[] pathsToFilters) throws FileNotFoundException, UnsupportedEncodingException{

        PrintWriter writer = new PrintWriter(pathToOutputFile, "UTF-8");
        ArrayList<Integer> idxs;
        String hybrid, protons, MULT = "", HSQC = "", COSY = "", BOND = "", HMBC = "";
        final int[][] bondTable = new int[this.mol.getAtomCount()][this.mol.getAtomCount()];
        for (int i = 0; i < this.mol.getAtomCount(); i++) {
            for (int j = 0; j < this.mol.getAtomCount(); j++) {
                bondTable[i][j] = 0;
            }
        }
        writer.println("; " + projectName);
        if(this.molFormula != null){
            writer.println("; " + MolecularFormulaManipulator.getString(this.molFormula) + "\n\n");
        } else {
            writer.println("; unknown molecular formula");
        }
        for (int i = 0; i < this.mol.getAtomCount(); i++) {
            // set MULT section in LSD input file
            // set hybridization level
            if(this.mol.getAtom(i).getHybridization() == null){
                hybrid = "X";
            } else {
                switch (this.mol.getAtom(i).getHybridization()) {
                    case SP1:
                    case S:
                        hybrid = "1";
                        break;
                    case SP2:
                        hybrid = "2";
                        break;
                    default:
                        hybrid = "3";
                }
            }
            // set implicit proton number
            if(this.mol.getAtom(i).getImplicitHydrogenCount() == null){
                protons = "X";
            } else {
                protons = String.valueOf(this.mol.getAtom(i).getImplicitHydrogenCount());
            }
            MULT += "MULT " + (i+1) + " " + this.mol.getAtom(i).getSymbol() + " " + hybrid + " " + protons;
            if(this.mol.getAtom(i).getProperty(casekit.NMR.Utils.getNMRShiftConstant(this.mol.getAtom(i).getSymbol())) != null){
                MULT += ";\t" + this.mol.getAtom(i).getProperty(casekit.NMR.Utils.getNMRShiftConstant(this.mol.getAtom(i).getSymbol()));
            }
            MULT += "\n";
            // set HSQC section in LSD input file
            if((this.mol.getAtom(i).getImplicitHydrogenCount() != null) && (this.mol.getAtom(i).getImplicitHydrogenCount() > 0)){
                HSQC += "HSQC " + (i+1) + " " + (i+1) + ";\t" + this.mol.getAtom(i).getSymbol() + "H" + this.mol.getAtom(i).getImplicitHydrogenCount() + "\n";
            }
            // set BOND section in LSD input file from INADEQUATE
            if (this.mol.getAtom(i).getProperty(CDKConstants.NMRSPECTYPE_2D_INADEQUATE) != null) {
                idxs = this.mol.getAtom(i).getProperty(CDKConstants.NMRSPECTYPE_2D_INADEQUATE);
                for (Integer idx : idxs) {
                    if (bondTable[i][idx] == 0 && bondTable[idx][i] == 0) {
                        bondTable[i][idx] = 1;
                        BOND += "BOND " + (i+1) + " " + (idx+1) + ";\t" + this.mol.getAtom(i).getSymbol() + "H" + this.mol.getAtom(i).getImplicitHydrogenCount() + " - " + this.mol.getAtom(idx).getSymbol() + "H" + this.mol.getAtom(idx).getImplicitHydrogenCount() + "\n";
                    }
                }
            }
            // set BOND section in LSD input file from COSY
            if(this.mol.getAtom(i).getProperty(CDKConstants.NMRSPECTYPE_2D_HHCOSY) != null){
                idxs = this.mol.getAtom(i).getProperty(CDKConstants.NMRSPECTYPE_2D_HHCOSY);
                for (Integer idx : idxs) {
                    if(bondTable[i][idx] == 0 && bondTable[idx][i] == 0){
                        bondTable[i][idx] = 1;
                        COSY += "COSY " + (i+1) + " " + (idx+1) + ";\t" + this.mol.getAtom(i).getSymbol() + "H" + this.mol.getAtom(i).getImplicitHydrogenCount() + " - " + this.mol.getAtom(idx).getSymbol() + "H" + this.mol.getAtom(idx).getImplicitHydrogenCount() + "\n";
                    } else {
                        COSY += ";COSY " + (i+1) + " " + (idx+1) + ";\t" + this.mol.getAtom(i).getSymbol() + "H" + this.mol.getAtom(i).getImplicitHydrogenCount() + " - " + this.mol.getAtom(idx).getSymbol() + "H" + this.mol.getAtom(idx).getImplicitHydrogenCount() + "\n";
                    }
                }
            }
            // set HMBC section in LSD input file
            // sets only HMBC signals which are not represented by a bond
            boolean test3JviaNextNeighborBond;
            if (this.mol.getAtom(i).getProperty(CDKConstants.NMRSPECTYPE_2D_HMBC) != null) {
                idxs = this.mol.getAtom(i).getProperty(CDKConstants.NMRSPECTYPE_2D_HMBC);
                for (Integer idx : idxs) {
                    if (bondTable[i][idx] == 0 && bondTable[idx][i] == 0) {
                        test3JviaNextNeighborBond = false;
                        for (IAtom neighbor : this.mol.getConnectedAtomsList(this.mol.getAtom(i))) {
                            if(this.mol.getBond(neighbor, this.mol.getAtom(idx)) != null){
                               test3JviaNextNeighborBond = true;
                               break;
                            }
                        }
                        if(test3JviaNextNeighborBond){
                            HMBC += ";HMBC " + (idx+1) + " " + (i+1) + "; 3J\t\n";
                        } else {
                            HMBC += "HMBC " + (idx+1) + " " + (i+1) + ";\n";
                        }
                    } else {
                        HMBC += ";HMBC " + (idx+1) + " " + (i+1) + "; 2J\t\n";
                    }
                }
            }
        }
        writer.println(MULT);
        writer.println(HSQC);
        writer.println(BOND);
        writer.println(COSY);
        writer.println(HMBC);
        
        String DEFF = "";
        String FEXP = "";
        if(pathsToFilters.length > 0){
            int fragmentCounter = 1;
            for (String pathToFilter : pathsToFilters) {
                File folder = new File(pathToFilter);
                File[] listOfFiles = folder.listFiles();
                for (File file : listOfFiles) {
                    if (file.isFile() && !file.getName().toLowerCase().contains(".")) {
                        DEFF += "DEFF F" + fragmentCounter + " \"" + file.getAbsolutePath() + "\"\n";
                        fragmentCounter++;
                    }
                }
            }
            FEXP = "FEXP \"NOT F1";
            for (int i = 2; i < fragmentCounter; i++) {
                FEXP += " and NOT F" + i;
            }
            FEXP += "\"";
        }
        
        writer.println(DEFF);
        writer.println(FEXP);
        writer.close();
        
    }
    
    
    public int[][] getNeighborhoodBondsCountMatrix(){
        
        return this.neighborhoodCountsMatrix;
    }
    
    
    
    
    
    public void countNeighborhoodBonds(final IAtomContainerSet acSet, final String[] bondsSet, final String elem, String[] neighborElems, final int minShift, final int maxShift, final int stepSize) throws FileNotFoundException, IOException{
        
        if (stepSize < 1) {
            System.err.println("stepSize < 1 not allowed!!!");
            return;
        }
        // creation of frequency counting matrix and shift indices holder
        this.neighborhoodCountsMatrix = new int[stepSize * (maxShift - minShift + 1)][3 + 4 + neighborElems.length * bondsSet.length];
        this.shiftIndicesInACSet.clear();
        for (int i = 0; i < stepSize * maxShift; i++) {
            for (int j = 0; j < 3 + 4 + neighborElems.length * bondsSet.length; j++) {
                neighborhoodCountsMatrix[i][j] = 0;
            }
            this.shiftIndicesInACSet.put(i, new ArrayList<>());
        }
        int atomIndexDB, shiftDBInt; double shiftDBDouble; IAtomContainer acDB;
        // go through all molecules in DB
        for (int k = 0; k < acSet.getAtomContainerCount(); k++) {
            acDB = acSet.getAtomContainer(k);
            // for all DB entries containing a spectrum for the current query atom type
            for (final String shiftsDB : casekit.NMR.DB.getSpectraFromNMRShiftDBEntry(acDB, elem)) {
                if (shiftsDB == null) {
                    continue;
                }
                String[][] shiftsDBvalues = casekit.NMR.Utils.parseShiftsNMRShiftDB(shiftsDB);
                for (String[] shiftsDBvalue : shiftsDBvalues) {
                    atomIndexDB = Integer.parseInt(shiftsDBvalue[2]);
                    // sometimes the DB atom index is wrong and out of array range 
                    if (atomIndexDB > acDB.getAtomCount() - 1) {
                        continue;
                    }
                    shiftDBDouble = Math.round(Double.parseDouble(shiftsDBvalue[0]) * stepSize) / (double) stepSize;;
                    // if DB shift value out of min-max-range then skip this shift 
                    if(shiftDBDouble < minShift || shiftDBDouble > maxShift - 1){
                        continue;
                    }
                    shiftDBInt = (int) (shiftDBDouble * stepSize);
                    this.neighborhoodCountsMatrix[shiftDBInt - minShift][0] += 1; // increase number of this shift occurence
                    this.neighborhoodCountsMatrix[shiftDBInt - minShift][1] += (acDB.getAtom(atomIndexDB).isInRing()) ? 1 : 0; // increase if atom is a ring member
                    this.neighborhoodCountsMatrix[shiftDBInt - minShift][2] += (acDB.getAtom(atomIndexDB).isAromatic()) ? 1 : 0; // increase if atom is aromatic
                    this.neighborhoodCountsMatrix[shiftDBInt - minShift][3] += ((acDB.getAtom(atomIndexDB).getImplicitHydrogenCount() != null) && (acDB.getAtom(atomIndexDB).getImplicitHydrogenCount() == 0)) ? 1 : 0; // qC count or equivalents, e.g. qN
                    this.neighborhoodCountsMatrix[shiftDBInt - minShift][4] += ((acDB.getAtom(atomIndexDB).getImplicitHydrogenCount() != null) && (acDB.getAtom(atomIndexDB).getImplicitHydrogenCount() == 1)) ? 1 : 0; // CH count or equivalents, e.g. NH
                    this.neighborhoodCountsMatrix[shiftDBInt - minShift][5] += ((acDB.getAtom(atomIndexDB).getImplicitHydrogenCount() != null) && (acDB.getAtom(atomIndexDB).getImplicitHydrogenCount() == 2)) ? 1 : 0; // CH2 count or equivalents, e.g. NH2
                    this.neighborhoodCountsMatrix[shiftDBInt - minShift][6] += ((acDB.getAtom(atomIndexDB).getImplicitHydrogenCount() != null) && (acDB.getAtom(atomIndexDB).getImplicitHydrogenCount() == 3)) ? 1 : 0; // CH3 count or equivalents, e.g. NH3
                    // add counts for a specific atom to matrix m
                    int[] counts = casekit.NMR.Utils.getNeighborhoodBondsCount(acDB, atomIndexDB, bondsSet, neighborElems);
                    for (int i = 0; i < counts.length; i++) {
                        this.neighborhoodCountsMatrix[shiftDBInt - minShift][3 + 4 + i] += counts[i];
                    }
                    // add this atom container index and atom index within it to belonging hash map
                    this.shiftIndicesInACSet.get(shiftDBInt).add(new Integer[]{k, atomIndexDB});
                }
            }
        }
    }
    
    
    
    
    
    
    
    
    
    
    
  
    
}
