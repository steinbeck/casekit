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
package casekit.nmr.convert;

import casekit.nmr.Utils;
import casekit.nmr.model.Spectrum;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.MolecularFormula;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class LSDConverter {

    /**
     *
     * @param projectName
     * @param pathToOutputFile
     * @param pathsToFilters
     * @param molecularFormula
     * @param ac
     * @param spectra
     *
     * @throws IOException
     */
    public static void ConvertToLSD(final String projectName, final String pathToOutputFile, final String[] pathsToFilters, final MolecularFormula molecularFormula, final IAtomContainer ac, final HashMap<String, Spectrum> spectra) throws IOException {

        String wholeContent, hybrid, protons, MULT = "", HSQC = "", BOND = "", HMBC = "", COSY = "";       
        wholeContent = "; project name: " + projectName + "\n";
        if(molecularFormula != null){
            wholeContent += "; molecular formula: " + MolecularFormulaManipulator.getString(molecularFormula) + "\n\n";
        } else {
            wholeContent += "; molecular formula: unknown \n\n";
        }
        for (int i = 0; i < ac.getAtomCount(); i++) {
            // set MULT section in LSD input file
            // set hybridization level
            if(ac.getAtom(i).getHybridization() == null){
                hybrid = "-";
            } else {
                switch (ac.getAtom(i).getHybridization()) {
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
            if(ac.getAtom(i).getImplicitHydrogenCount() == null){
                protons = "-";
            } else {
                protons = String.valueOf(ac.getAtom(i).getImplicitHydrogenCount());
            }
            MULT += "MULT " + (i+1) + " " + ac.getAtom(i).getSymbol() + " " + hybrid + " " + protons;
            if(ac.getAtom(i).getProperty(Utils.getNMRShiftConstant(ac.getAtom(i).getSymbol())) != null){
                String hCount;
                if(ac.getAtom(i).getImplicitHydrogenCount() == null){
                    hCount = "x";
                } else {
                    hCount = String.valueOf(ac.getAtom(i).getImplicitHydrogenCount());
                }
                MULT += ";\t" + ac.getAtom(i).getProperty(Utils.getNMRShiftConstant(ac.getAtom(i).getSymbol())) + ",\t" + ac.getAtom(i).getSymbol() + "H" + hCount;
            }
            MULT += "\n";
            // set HSQC section in LSD input file
            if((ac.getAtom(i).getImplicitHydrogenCount() != null) && (ac.getAtom(i).getImplicitHydrogenCount() > 0)){
                HSQC += "HSQC " + (i+1) + " " + (i+1) + ";\t" + ac.getAtom(i).getSymbol() + "H" + ac.getAtom(i).getImplicitHydrogenCount() + "\n";
            }   
        }
        wholeContent += MULT + "\n";
        wholeContent += HSQC + "\n";
        
        // set BOND information in LSD input file by INADEQUATE or general bond knowledge
        for (IBond bond : ac.bonds()) {
            BOND += "BOND " + (bond.getAtom(0).getIndex()+1) + " " + (bond.getAtom(1).getIndex()+1) + ";\t" + ac.getAtom(bond.getAtom(0).getIndex()).getSymbol() + "H" + ac.getAtom(bond.getAtom(0).getIndex()).getImplicitHydrogenCount() + " - " + ac.getAtom(bond.getAtom(1).getIndex()).getSymbol() + "H" + ac.getAtom(bond.getAtom(1).getIndex()).getImplicitHydrogenCount() + "\n";
        }
        wholeContent += BOND + "\n";


        // @TODO repair HMBC and COSY information output
             
//        // set HMBC information to LSD input file
//        ArrayList<Integer> indicesInAtomContainerDim1;
//        ArrayList<Integer> indicesInAtomContainerDim2;
//        final boolean [][] HMBCTable = new boolean[ac.getAtomCount()][ac.getAtomCount()];
//        for (int i = 0; i < ac.getAtomCount(); i++) {
//            for (int j = 0; j < ac.getAtomCount(); j++) {
//                HMBCTable[i][j] = false;
//            }
//        }
//        for (final Spectrum spectrum : spectra.values()) {
//            if((spectrum.getNDim() != 2) || !spectrum.getSpecType().startsWith(CDKConstants.NMRSPECTYPE_2D_HMBC)){
//                continue;
//            }
//            indicesInAtomContainerDim1 = this.getAssignedAtomIndices(spectrum, 0);
//            indicesInAtomContainerDim2 = this.getAssignedAtomIndices(spectrum, 1);
//            HMBC += ";\t " + spectrum.getSpecType() + " " + Utils.getSpectrumNucleiAsString(spectrum) + "\n";
//            for (int i = 0; i < spectrum.getSignalCount(); i++) {
//                if((indicesInAtomContainerDim1.get(i) > -1) && (indicesInAtomContainerDim2.get(i) > -1)){
//                    // set signal only if it is not already covered by BOND
//                    // here reversed order (see LSD manual page): 1. heavy atom, 2. proton
//                    if(ac.getBond(ac.getAtom(indicesInAtomContainerDim2.get(i)), ac.getAtom(indicesInAtomContainerDim1.get(i))) != null){
//                        HMBC += ";";
//                    }
//                    HMBC += "HMBC " + (indicesInAtomContainerDim2.get(i) + 1) + " " + (indicesInAtomContainerDim1.get(i) + 1) + ";\t" + ac.getAtom(indicesInAtomContainerDim2.get(i)).getSymbol() + "H" + ac.getAtom(indicesInAtomContainerDim2.get(i)).getImplicitHydrogenCount() + " - " + ac.getAtom(indicesInAtomContainerDim1.get(i)).getSymbol() + "H" + ac.getAtom(indicesInAtomContainerDim1.get(i)).getImplicitHydrogenCount() + "\n";
//                    HMBCTable[indicesInAtomContainerDim2.get(i)][indicesInAtomContainerDim1.get(i)] = true;
//                }
//            }
//        }
//        wholeContent += HMBC + "\n";
//        // set COSY information to LSD input file
//        for (final Spectrum spectrum : spectra.values()) {
//            if((spectrum.getNDim() != 2) || !spectrum.getSpecType().startsWith(CDKConstants.NMRSPECTYPE_2D_HHCOSY)){
//                continue;
//            }
//            indicesInAtomContainerDim1 = this.getAssignedAtomIndices(spectrum, 0);
//            indicesInAtomContainerDim2 = this.getAssignedAtomIndices(spectrum, 1);
//            COSY += ";\t " + spectrum.getSpecType() + " " + Utils.getSpectrumNucleiAsString(spectrum) + "\n";
//            for (int i = 0; i < spectrum.getSignalCount(); i++) {
//                if((indicesInAtomContainerDim1.get(i) > -1) && (indicesInAtomContainerDim2.get(i) > -1)){
//                    // set signal only if it is not already covered by BOND or HMBC
//                    if((ac.getBond(ac.getAtom(indicesInAtomContainerDim1.get(i)), ac.getAtom(indicesInAtomContainerDim2.get(i))) != null)
//                            || HMBCTable[indicesInAtomContainerDim1.get(i)][indicesInAtomContainerDim2.get(i)]){
//                        COSY += ";";
//                    }
//                    COSY += "COSY " + (indicesInAtomContainerDim1.get(i) + 1) + " " + (indicesInAtomContainerDim2.get(i) + 1) + ";\t" + ac.getAtom(indicesInAtomContainerDim1.get(i)).getSymbol() + "H" + ac.getAtom(indicesInAtomContainerDim1.get(i)).getImplicitHydrogenCount() + " - " + ac.getAtom(indicesInAtomContainerDim2.get(i)).getSymbol() + "H" + ac.getAtom(indicesInAtomContainerDim2.get(i)).getImplicitHydrogenCount() + "\n";
//                }
//            }
//        }
//        wholeContent += COSY + "\n";
        // set filter definitions
        String DEFF = "";
        String FEXP = "";
        if((pathsToFilters != null) && pathsToFilters.length > 0){
            int fragmentCounter = 1;
            for (final String pathToFilter : pathsToFilters) {
                File folder = new File(pathToFilter);
                File[] listOfFiles = folder.listFiles();
                for (final File file : listOfFiles) {
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

    
    
//    public static void countNeighborhoodBonds(final String pathToNMRShiftDB, final String[] bondsSet, final String nucleus, final ArrayList<String> neighborElems, final int minShift, final int maxShift, final int stepSize) throws FileNotFoundException, IOException, CDKException {
//
//        if (stepSize < 1) {
//            System.err.println("stepSize < 1 not allowed!!!");
//            return;
//        }
//        // creation of frequency counting matrix and shift indices holder
//        final int[][] neighborhoodCountsMatrix = new int[stepSize * (maxShift - minShift + 1)][3 + 4 + neighborElems.size() * bondsSet.length];
//        final IAtomContainerSet acSet = NMRShiftDB.getStructuresFromSDFile(pathToNMRShiftDB, true);
//        final HashMap<Integer, ArrayList<Integer[]>> shiftIndicesInACSet = new HashMap<>();
//        for (int i = 0; i < stepSize * maxShift; i++) {
//            for (int j = 0; j < 3 + 4 + neighborElems.size() * bondsSet.length; j++) {
//                neighborhoodCountsMatrix[i][j] = 0;
//            }
//            shiftIndicesInACSet.put(i, new ArrayList<>());
//        }
//        int atomIndexDB, shiftDBInt; double shiftDBDouble; IAtomContainer acDB;
//        // go through all molecules in MongoDB
//        for (int k = 0; k < acSet.getAtomContainerCount(); k++) {
//            acDB = acSet.getAtomContainer(k);
//            // for all MongoDB entries containing a spectrum for the current query atom type
//            for (final String shiftsDB : NMRShiftDB.getSpectraFromNMRShiftDB(pathToNMRShiftDB, nucleus)) {
//                if (shiftsDB == null) {
//                    continue;
//                }
//                String[][] shiftsDBvalues = casekit.casekit.nmr.dbservice.NMRShiftDB.parseNMRShiftDBSpectrum(shiftsDB);
//                for (String[] shiftsDBvalue : shiftsDBvalues) {
//                    atomIndexDB = Integer.parseInt(shiftsDBvalue[2]);
//                    // sometimes the MongoDB atom index is wrong and out of array range
//                    if (atomIndexDB > acDB.getAtomCount() - 1) {
//                        continue;
//                    }
//                    shiftDBDouble = Math.round(Double.parseDouble(shiftsDBvalue[0]) * stepSize) / (double) stepSize;
//                    // if MongoDB shift value out of min-max-range then skip this shift
//                    if(shiftDBDouble < minShift || shiftDBDouble > maxShift - 1){
//                        continue;
//                    }
//                    shiftDBInt = (int) (shiftDBDouble * stepSize);
//                    neighborhoodCountsMatrix[shiftDBInt - minShift][0] += 1; // increase number of this shift occurence
//                    neighborhoodCountsMatrix[shiftDBInt - minShift][1] += (acDB.getAtom(atomIndexDB).isInRing()) ? 1 : 0; // increase if atom is a ring member
//                    neighborhoodCountsMatrix[shiftDBInt - minShift][2] += (acDB.getAtom(atomIndexDB).isAromatic()) ? 1 : 0; // increase if atom is aromatic
//                    neighborhoodCountsMatrix[shiftDBInt - minShift][3] += ((acDB.getAtom(atomIndexDB).getImplicitHydrogenCount() != null) && (acDB.getAtom(atomIndexDB).getImplicitHydrogenCount() == 0)) ? 1 : 0; // qC count or equivalents, e.g. qN
//                    neighborhoodCountsMatrix[shiftDBInt - minShift][4] += ((acDB.getAtom(atomIndexDB).getImplicitHydrogenCount() != null) && (acDB.getAtom(atomIndexDB).getImplicitHydrogenCount() == 1)) ? 1 : 0; // CH count or equivalents, e.g. NH
//                    neighborhoodCountsMatrix[shiftDBInt - minShift][5] += ((acDB.getAtom(atomIndexDB).getImplicitHydrogenCount() != null) && (acDB.getAtom(atomIndexDB).getImplicitHydrogenCount() == 2)) ? 1 : 0; // CH2 count or equivalents, e.g. NH2
//                    neighborhoodCountsMatrix[shiftDBInt - minShift][6] += ((acDB.getAtom(atomIndexDB).getImplicitHydrogenCount() != null) && (acDB.getAtom(atomIndexDB).getImplicitHydrogenCount() == 3)) ? 1 : 0; // CH3 count or equivalents, e.g. NH3
//                    // add counts for a specific atom to matrix m
//                    int[] counts = casekit.casekit.nmr.Utils.getNeighborhoodBondsCount(acDB, atomIndexDB, bondsSet, neighborElems);
//                    for (int i = 0; i < counts.length; i++) {
//                        neighborhoodCountsMatrix[shiftDBInt - minShift][3 + 4 + i] += counts[i];
//                    }
//                    // add this atom container index and atom index within it to belonging hash map
//                    shiftIndicesInACSet.get(shiftDBInt).add(new Integer[]{k, atomIndexDB});
//                }
//            }
//        }
//    }

}
