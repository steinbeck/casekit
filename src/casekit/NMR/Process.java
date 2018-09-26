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
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Process extends ParseRawData {
    
    private final IAtomContainer mol;
    private final IMolecularFormula molFormula;
    private final HashMap<String, ArrayList<Integer>> atomTypeIndices;    
    private int[][] neighborhoodCountsMatrix;
    private final HashMap<Integer, ArrayList<Integer[]>> shiftIndicesInACSet = new HashMap<>(); // holding of all indices of each ac set (DB) entry [first value] and it's atom indices [second value] too
    
    
    public Process(){
        super();
        this.molFormula = super.getMolecularFormula();
        this.mol = super.getAtomContainer();
        this.atomTypeIndices = super.getAtomTypeIndices();
    }
    
    public Process(final IMolecularFormula molFormula){
        super(molFormula);
        this.molFormula = super.getMolecularFormula();
        this.mol = super.getAtomContainer();
        this.atomTypeIndices = super.getAtomTypeIndices();
    }
    
    
    /**
     *
     * @param projectName
     * @param pathToOutputFile
     * @param pathsToFilters
     * @throws FileNotFoundException
     * @throws UnsupportedEncodingException
     */
    public void createLSDInputFile(final String projectName, final String pathToOutputFile, final String[] pathsToFilters) throws FileNotFoundException, UnsupportedEncodingException, IOException{        
      
        String wholeContent, hybrid, protons, MULT = "", HSQC = "", BOND = "", HMBC = "", COSY = "";       
        wholeContent = "; project name: " + projectName + "\n";
        if(this.molFormula != null){
            wholeContent += "; molecular formula: " + MolecularFormulaManipulator.getString(this.molFormula) + "\n\n";
        } else {
            wholeContent += "; molecular formula: unknown \n\n";
        }
        for (int i = 0; i < this.mol.getAtomCount(); i++) {
            // set MULT section in LSD input file
            // set hybridization level
            if(this.mol.getAtom(i).getHybridization() == null){
                hybrid = "-";
            } else {
                switch (this.mol.getAtom(i).getHybridization()) {
                    case SP1:
                    case S:
                        hybrid = "1"; break;
                    case SP2:
                        hybrid = "2"; break;
                    default:
                        hybrid = "3";
                }
            }
            // set implicit proton number
            if(this.mol.getAtom(i).getImplicitHydrogenCount() == null){
                protons = "-";
            } else {
                protons = String.valueOf(this.mol.getAtom(i).getImplicitHydrogenCount());
            }
            MULT += "MULT " + (i+1) + " " + this.mol.getAtom(i).getSymbol() + " " + hybrid + " " + protons;
            if(this.mol.getAtom(i).getProperty(Utils.getNMRShiftConstant(this.mol.getAtom(i).getSymbol())) != null){
                String hCount;
                if(this.mol.getAtom(i).getImplicitHydrogenCount() == null){
                    hCount = "x";
                } else {
                    hCount = String.valueOf(this.mol.getAtom(i).getImplicitHydrogenCount());
                }
                MULT += ";\t" + this.mol.getAtom(i).getProperty(Utils.getNMRShiftConstant(this.mol.getAtom(i).getSymbol())) + ",\t" + this.mol.getAtom(i).getSymbol() + "H" + hCount;
            }
            MULT += "\n";
            // set HSQC section in LSD input file
            if((this.mol.getAtom(i).getImplicitHydrogenCount() != null) && (this.mol.getAtom(i).getImplicitHydrogenCount() > 0)){
                HSQC += "HSQC " + (i+1) + " " + (i+1) + ";\t" + this.mol.getAtom(i).getSymbol() + "H" + this.mol.getAtom(i).getImplicitHydrogenCount() + "\n";
            }   
        }
        wholeContent += MULT + "\n";
        wholeContent += HSQC + "\n";
        
        // set BOND information in LSD input file by INADEQUATE
        for (IBond bond : this.mol.bonds()) {
            BOND += "BOND " + (bond.getAtom(0).getIndex()+1) + " " + (bond.getAtom(1).getIndex()+1) + ";\t" + this.mol.getAtom(bond.getAtom(0).getIndex()).getSymbol() + "H" + this.mol.getAtom(bond.getAtom(0).getIndex()).getImplicitHydrogenCount() + " - " + this.mol.getAtom(bond.getAtom(1).getIndex()).getSymbol() + "H" + this.mol.getAtom(bond.getAtom(1).getIndex()).getImplicitHydrogenCount() + "\n";
        }
        wholeContent += BOND + "\n";
             
        // set HMBC information to LSD input file
        ArrayList<Integer> indicesInAtomContainerDim1;
        ArrayList<Integer> indicesInAtomContainerDim2;
        final boolean [][] HMBCTable = new boolean[this.mol.getAtomCount()][this.mol.getAtomCount()];
        for (int i = 0; i < this.mol.getAtomCount(); i++) {
            for (int j = 0; j < this.mol.getAtomCount(); j++) {
                HMBCTable[i][j] = false;
            }
        }
        for (final Spectrum spectrum : this.getSpectra().values()) {
            if((spectrum.getDimCount() != 2) || !spectrum.getSpecType().startsWith(CDKConstants.NMRSPECTYPE_2D_HMBC)){
                continue;
            }
            indicesInAtomContainerDim1 = this.getAssignedAtomIndices(spectrum, 0);
            indicesInAtomContainerDim2 = this.getAssignedAtomIndices(spectrum, 1);
            HMBC += ";\t " + spectrum.getSpecType() + " " + Utils.getSpectrumNucleiAsString(spectrum) + "\n";                        
            for (int i = 0; i < spectrum.getSignalCount(); i++) {
                if((indicesInAtomContainerDim1.get(i) > -1) && (indicesInAtomContainerDim2.get(i) > -1)){
                    // set signal only if it is not already covered by BOND
                    // here reversed order (see LSD manual page): 1. heavy atom, 2. proton
                    if(this.mol.getBond(this.mol.getAtom(indicesInAtomContainerDim2.get(i)), this.mol.getAtom(indicesInAtomContainerDim1.get(i))) != null){
                        HMBC += ";";
                    }
                    HMBC += "HMBC " + (indicesInAtomContainerDim2.get(i) + 1) + " " + (indicesInAtomContainerDim1.get(i) + 1) + ";\t" + this.mol.getAtom(indicesInAtomContainerDim2.get(i)).getSymbol() + "H" + this.mol.getAtom(indicesInAtomContainerDim2.get(i)).getImplicitHydrogenCount() + " - " + this.mol.getAtom(indicesInAtomContainerDim1.get(i)).getSymbol() + "H" + this.mol.getAtom(indicesInAtomContainerDim1.get(i)).getImplicitHydrogenCount() + "\n";
                    HMBCTable[indicesInAtomContainerDim2.get(i)][indicesInAtomContainerDim1.get(i)] = true;
                }                   
            }                            
        }  
        wholeContent += HMBC + "\n";
        // set COSY information to LSD input file
        for (final Spectrum spectrum : this.getSpectra().values()) {
            if((spectrum.getDimCount() != 2) || !spectrum.getSpecType().startsWith(CDKConstants.NMRSPECTYPE_2D_HHCOSY)){
                continue;
            }
            indicesInAtomContainerDim1 = this.getAssignedAtomIndices(spectrum, 0);
            indicesInAtomContainerDim2 = this.getAssignedAtomIndices(spectrum, 1);
            COSY += ";\t " + spectrum.getSpecType() + " " + Utils.getSpectrumNucleiAsString(spectrum) + "\n";                        
            for (int i = 0; i < spectrum.getSignalCount(); i++) {
                if((indicesInAtomContainerDim1.get(i) > -1) && (indicesInAtomContainerDim2.get(i) > -1)){
                    // set signal only if it is not already covered by BOND or HMBC
                    if((this.mol.getBond(this.mol.getAtom(indicesInAtomContainerDim1.get(i)), this.mol.getAtom(indicesInAtomContainerDim2.get(i))) != null)
                            || HMBCTable[indicesInAtomContainerDim1.get(i)][indicesInAtomContainerDim2.get(i)]){
                        COSY += ";";
                    }
                    COSY += "COSY " + (indicesInAtomContainerDim1.get(i) + 1) + " " + (indicesInAtomContainerDim2.get(i) + 1) + ";\t" + this.mol.getAtom(indicesInAtomContainerDim1.get(i)).getSymbol() + "H" + this.mol.getAtom(indicesInAtomContainerDim1.get(i)).getImplicitHydrogenCount() + " - " + this.mol.getAtom(indicesInAtomContainerDim2.get(i)).getSymbol() + "H" + this.mol.getAtom(indicesInAtomContainerDim2.get(i)).getImplicitHydrogenCount() + "\n";                    
                }                   
            }                            
        }  
        wholeContent += COSY + "\n";
        // set filter definitions
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
        
        wholeContent += DEFF + "\n";
        wholeContent += FEXP + "\n";
        
        Utils.writeTextFile(pathToOutputFile, wholeContent);
    }
    
    
    public int[][] getNeighborhoodBondsCountMatrix(){
        
        return this.neighborhoodCountsMatrix;
    }
    
    
    
    
    
    public void countNeighborhoodBonds(final IAtomContainerSet acSet, final String[] bondsSet, final String elem, final ArrayList<String> neighborElems, final int minShift, final int maxShift, final int stepSize) throws FileNotFoundException, IOException{
        
        if (stepSize < 1) {
            System.err.println("stepSize < 1 not allowed!!!");
            return;
        }
        // creation of frequency counting matrix and shift indices holder
        this.neighborhoodCountsMatrix = new int[stepSize * (maxShift - minShift + 1)][3 + 4 + neighborElems.size() * bondsSet.length];
        this.shiftIndicesInACSet.clear();
        for (int i = 0; i < stepSize * maxShift; i++) {
            for (int j = 0; j < 3 + 4 + neighborElems.size() * bondsSet.length; j++) {
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
                String[][] shiftsDBvalues = casekit.NMR.Utils.parseShiftsInNMRShiftDBEntry(shiftsDB);
                for (String[] shiftsDBvalue : shiftsDBvalues) {
                    atomIndexDB = Integer.parseInt(shiftsDBvalue[2]);
                    // sometimes the DB atom index is wrong and out of array range 
                    if (atomIndexDB > acDB.getAtomCount() - 1) {
                        continue;
                    }
                    shiftDBDouble = Math.round(Double.parseDouble(shiftsDBvalue[0]) * stepSize) / (double) stepSize;
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
