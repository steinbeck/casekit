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

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Process extends ParseRawData {
    
    final private IAtomContainer mol;
    final private IMolecularFormula molFormula;
    private final HashMap<String, ArrayList<Integer>> atomTypeIndices = new HashMap<>();    
    private int[][] neighborhoodCountsMatrix;
    final private HashMap<Integer, ArrayList<Integer[]>> shiftIndicesInACSet = new HashMap<>(); // holding of all indices of each ac set (DB) entry [first value] and it's atom indices [second value] too
    
    
    public Process(){
        super();
        this.molFormula = super.getMolecularFormula();
        this.mol = super.getAtomContainer();
    }
    
    public Process(final IMolecularFormula molFormula){
        super(molFormula);
        this.molFormula = super.getMolecularFormula();
        this.mol = super.getAtomContainer();
        this.setAtomTypeIndices();
    }

    
     /**
     * Sets bonds from already set experiment information (H,H-COSY, INADEQUATE
     * and HMBC).
     * 
     * @param experiments
     */
    public void setBonds(final String[] experiments) {

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
                            if (experiments[e].equals(CDKConstants.NMRSPECTYPE_2D_HMBC)) {
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

    private void setBond(final int index1, final int index2) {

        if (this.mol.getBond(this.mol.getAtom(index1), this.mol.getAtom(index2)) != null) {
            this.mol.removeBond(this.mol.getAtom(index1), this.mol.getAtom(index2));
        }
        this.mol.addBond(index1, index2, casekit.NMR.Utils.getBondTypeFromHybridizations(this.mol.getAtom(index1), this.mol.getAtom(index2)));
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
