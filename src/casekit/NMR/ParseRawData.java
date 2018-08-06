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
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import javax.xml.parsers.ParserConfigurationException;
import org.openscience.cdk.Atom;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.xml.sax.SAXException;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class ParseRawData {
    
    final private IAtomContainer mol;
    final private IMolecularFormula molFormula;
    private HashMap<String, ArrayList<Integer>> atomTypeIndices;
    
    public final static String PROP_EQUIVALENCE = "equivalences";
    

    /**
     * Creates an instances of this class with an empty class atom container.
     */
    public ParseRawData(){
        this.molFormula = null;
        this.mol = SilentChemObjectBuilder.getInstance().newAtomContainer();
        this.setAtomTypeIndices();
    }
    
    /**
     * Creates an instances of this class with a class atom container consisting 
     * of all heavy atoms in given molecular formula.
     * 
     * @param molFormula IMolecularFormula object for IAtomContainer creation
     */
    public ParseRawData(final IMolecularFormula molFormula){
        this.molFormula = molFormula;
        this.mol = Utils.removeAtoms(MolecularFormulaManipulator.getAtomContainer(this.molFormula), "H");
        this.setAtomTypeIndices();
    }
    
    
    /**
     * Returns used IMolecularFormula object for this class instance.
     * 
     * @return 
     */
    public final IMolecularFormula getMolecularFormula() {

        return this.molFormula;
    }

    
    /**
     * Returns used IAtomContainer object for this class instance.
     *
     * @return 
     */
    public final IAtomContainer getAtomContainer() {

        return this.mol;
    }

    
    /**
     * Returns a HashMap object with the indices of all atoms for all atom types
     * (elements) within the atom container of this class.
     *
     * @return
     */
    public final HashMap<String, ArrayList<Integer>> getAtomTypeIndices() {

        return this.atomTypeIndices;
    }
    
    
    /**
     * Sets the indices of all atoms in this class atom container.
     * @see Utils#getAtomTypeIndices(org.openscience.cdk.interfaces.IAtomContainer) 
     *
     */
    public final void setAtomTypeIndices(){
        
        this.atomTypeIndices = Utils.getAtomTypeIndices(this.mol);
    }
    
    
    /**
     * Copies all up to here set properties from an atom in atom container to its 
     * linked atoms with equivalent shift values.
     *
     */
    public final void setEquivalentProperties() {
        
        Map<Object, Object> properties;
        for (int i = 0; i < this.mol.getAtomCount(); i++) {
            if (this.mol.getAtom(i).getProperty(ParseRawData.PROP_EQUIVALENCE) != null) {
                properties = this.mol.getAtom(i).getProperties();
                for (final Object prop: properties.keySet()) {
                    if (this.mol.getAtom(i).getProperty(prop) != null && !prop.equals(ParseRawData.PROP_EQUIVALENCE)) {
                        for (final int k : (ArrayList<Integer>) this.mol.getAtom(i).getProperty(ParseRawData.PROP_EQUIVALENCE)) {
                            this.mol.getAtom(k).setProperty(prop, this.mol.getAtom(i).getProperty(prop));
                        }
                    }
                }
            }
        }
    }
    
    
    /**
     * Creates a Spectrum class object from 1D NMR peak list in CSV file format.
     * 
     * @param pathToCSV Path to peak list (Bruker's TopSpin csv file
     * format)
     * @param atomType Element name (e.g. "C") which also occurrs in
     * {@link Utils#getNMRShiftConstant(java.lang.String)}
     * @return Spectrum class object from given input file
     * @throws java.io.IOException
     */
    public final Spectrum parse1DNMRviaCSV(final String pathToCSV, final String atomType) throws IOException {
        
          final Spectrum spectrum = Utils.CSVtoSpectrum(pathToCSV, new int[]{4}, new String[]{atomType}, 6);
          
          return spectrum;
    }

    
    /**
     * Creates a Spectrum class object from 1D NMR peak list in XML file format.
     * 
     * @param pathToXML Path to peak list (Bruker's TopSpin csv file
     * format)
     * @param atomType Element name (e.g. "C") which also occurrs in
     * {@link Utils#getNMRShiftConstant(java.lang.String)}
     * @return Spectrum class object from given input file
     * @throws java.io.IOException
     * @throws javax.xml.parsers.ParserConfigurationException
     * @throws org.xml.sax.SAXException
     */
    public final Spectrum parse1DNMRviaXML(final String pathToXML, final String atomType) throws IOException, ParserConfigurationException, SAXException {

        final Spectrum spectrum = Utils.XMLtoSpectrum(pathToXML, 1, new int[]{1}, new String[]{atomType});

        return spectrum;
    }

    
    /**
     * Sets the 1D NMR shift values for given Spectrum object to atoms of the class IAtomContainer.
     * The shift values will be assigned sequentially.  
     * In case of a molecular formula is given in this class, the number of 
     * shifts must be equal to the number of atoms in this molecular formula. 
     * For less shifts in shift list you will be asked for entering equivalences.
     * Otherwise this function will return a false value.
     * In case of no molecular was given to this class, a new atom in the atom container 
     * will be created regarding to the input shift list.
     * Each shift value is set to {@link IAtomContainer#setProperty(java.lang.Object, java.lang.Object) 
     * as result of Utils#getNMRShiftConstant(java.lang.String)}, depending on 
     * the specified atom type (element).
     * After usage of this function, the input Spectrum class object might be extended during 
     * equivalent signal selection by user and contains the signal assignment indices 
     * in atom container.
     * 
     * @param spectrum Spectrum class object containing the 1D shift information
     * @return 
     * @throws java.io.IOException
     */
    public final boolean set1DNMR(final Spectrum spectrum) throws IOException{
        final String atomType = Utils.getElementIdentifier(spectrum.getNuclei()[0]);
        ArrayList<Double> shifts = spectrum.getShiftsByDim(0);
        // check whether indices for that atom type exist or the number of input signals are greater than the atom number in atom container for that atom type
        if (!this.atomTypeIndices.containsKey(atomType) || shifts.size() > this.atomTypeIndices.get(atomType).size()) {
            // if molecular formula is known and too much picked peaks are to be assigned
            if(this.atomTypeIndices.containsKey(atomType) || MolecularFormulaManipulator.getElementCount(this.molFormula, atomType) == 0){
                System.err.println("Too many peaks in peak list for \"" + atomType + "\" and molecular formula \"" + MolecularFormulaManipulator.getString(this.molFormula) + "\"!!!");
                return false;
            } else { // 
                // "fill up" the first peaks for that atom type from given peak list
                IAtom atom;
                for (final double shift : shifts) {
                    atom = new Atom(atomType);
                    atom.setProperty(casekit.NMR.Utils.getNMRShiftConstant(atomType), shift);
                    atom.setImplicitHydrogenCount(null);
                    this.mol.addAtom(atom);
                }
                this.setAtomTypeIndices();
            }
        }
        // assign shift values to atoms sequentially
        int assignedShiftCount = 0;
        for (final int i : this.atomTypeIndices.get(atomType)) {
            if(assignedShiftCount < shifts.size()){
                // shift assignment in atom
                this.mol.getAtom(i).setProperty(casekit.NMR.Utils.getNMRShiftConstant(atomType), shifts.get(assignedShiftCount));
                spectrum.getSignal(assignedShiftCount).setAssignedAtomIndex(i, 0);
            }
            assignedShiftCount++;
        }
        // "fill up" the missing equivalent peaks
        // check whether the number of input signals is smaller than the number of atoms in atom container from that atom type
        if (shifts.size() < this.atomTypeIndices.get(atomType).size()) {
            System.out.println("Not enough peaks in 1D peak list for \"" + atomType + "\"!!!");
            this.askForEquivalentPeaks(spectrum, atomType);
        }
        
        return true;
    }
    
    
    private void askForEquivalentPeaks(final Spectrum spectrum, final String atomType) throws IOException {

        BufferedReader br = new BufferedReader(new InputStreamReader(System.in)); int n;
        final ArrayList<Integer> validIndices = new ArrayList<>();
        // walk through all atoms of given atom type
        for (final int i : this.atomTypeIndices.get(atomType)) {
            // ignore atoms with already set NMR shift value
            if (this.mol.getAtom(i).getProperty(Utils.getNMRShiftConstant(atomType)) != null) {
                continue;
            }
            // display all selectable atom indices in atom container
            System.out.println("\nThe " + i + "th shift value is missing!\nWhich shift value is not unique?");
            for (final int k : this.atomTypeIndices.get(atomType)) {
                if(this.mol.getAtom(k).getProperty(Utils.getNMRShiftConstant(atomType)) != null){
                    System.out.println(k + "\t: " + this.mol.getAtom(k).getProperty(Utils.getNMRShiftConstant(atomType)));
                    validIndices.add(k);
                }
            }
            // get selected index by user input
            n = -1;
            while(!validIndices.contains(n)){
                System.out.println("Enter the index: ");
                n = Integer.parseInt(br.readLine());
            }
            
            spectrum.addSignal(spectrum.getSignal(validIndices.indexOf(n)));
            // copy NMR shift value from already set value in selected atom n to unset value in atom i
            this.mol.getAtom(i).setProperty(Utils.getNMRShiftConstant(atomType), this.mol.getAtom(n).getProperty(Utils.getNMRShiftConstant(atomType)));
            
            
            if(this.mol.getAtom(i).getProperty(ParseRawData.PROP_EQUIVALENCE) == null){
                this.mol.getAtom(i).setProperty(ParseRawData.PROP_EQUIVALENCE, new ArrayList<>());
            }
            if(this.mol.getAtom(n).getProperty(ParseRawData.PROP_EQUIVALENCE) == null){
                this.mol.getAtom(n).setProperty(ParseRawData.PROP_EQUIVALENCE, new ArrayList<>());
            }
            ((ArrayList<Integer>) this.mol.getAtom(i).getProperty(ParseRawData.PROP_EQUIVALENCE)).add(n);
            ((ArrayList<Integer>) this.mol.getAtom(n).getProperty(ParseRawData.PROP_EQUIVALENCE)).add(i);
        }
    }
    
    
    /**
     * Sets the number of implicit hydrogens from two carbon DEPT90 and DEPT135
     * peak
     * tables to carbon atoms. The meanwhile found matches are corrected,
     * see
     * {@link Utils#correctShiftMatches(IAtomContainer, ArrayList, ArrayList, double,String)}.
     *
     * @param pathToDEPT Path to one DEPT peak list (Bruker's TopSpin csv file
     * format)
     * @return 
     * @throws java.io.IOException
     */
    public final Spectrum parseDEPTviaCSV(final String pathToDEPT) throws IOException {
        
        final Spectrum spectrum = Utils.CSVtoSpectrum(pathToDEPT, new int[]{4}, new String[]{"C"}, 6);

        return spectrum;
    }

    /**
     * Sets the number of implicit hydrogens from two carbon DEPT90 and DEPT135
     * XML files to carbon atoms. The meanwhile found matches are corrected, see
     * {@link Utils#correctShiftMatches(IAtomContainer, ArrayList, ArrayList, double,String)}.
     *
     * @param pathToDEPT Path to one DEPT peak list (Bruker's TopSpin XML file
     * format)
     * @return 
     * @throws java.io.IOException
     * @throws javax.xml.parsers.ParserConfigurationException
     * @throws org.xml.sax.SAXException
     */
    public final Spectrum parseDEPTviaXML(final String pathToDEPT) throws IOException, ParserConfigurationException, SAXException {

        final Spectrum spectrum = Utils.XMLtoSpectrum(pathToDEPT, 1, new int[]{1}, new String[]{"C"});

        return spectrum;
    }
    
    /**
     * Sets the hydrogen count information of carbon atoms in atom conatiner
     * by usage of DEPT90 and DEPT135 information.
     *
     * @param spectrumDEPT90 DEPT90 spectrum
     * @param spectrumDEPT135 DEPT135 spectrum which has to contain intensity 
     * information
     * @param tol tolerance value [ppm] for carbon shift matching
     * @return false if one of the spectra is not set or the intensities are missing
     */
    public final int setDEPT(final Spectrum spectrumDEPT90, final Spectrum spectrumDEPT135, final double tol){
        
        if(spectrumDEPT90 == null || spectrumDEPT135 == null || spectrumDEPT135.getIntensities() == null){
            return 0;
        }
        final ArrayList<Double> shiftsDEPT90 = spectrumDEPT90.getShiftsByDim(0);
        final ArrayList<Double> shiftsDEPT135 = spectrumDEPT135.getShiftsByDim(0);
        final ArrayList<Double> intensitiesDEPT135 = spectrumDEPT135.getIntensities();
        ArrayList<Integer> matchesDEPT90 = casekit.NMR.Utils.findShiftMatches(this.mol, shiftsDEPT90, tol, "C");
        matchesDEPT90 = casekit.NMR.Utils.correctShiftMatches(this.mol, shiftsDEPT90, matchesDEPT90, tol, "C");
        ArrayList<Integer> matchesDEPT135 = casekit.NMR.Utils.findShiftMatches(this.mol, shiftsDEPT135, tol, "C");
        matchesDEPT135 = casekit.NMR.Utils.correctShiftMatches(this.mol, shiftsDEPT135, matchesDEPT135, tol, "C");
        
        spectrumDEPT90.setAssignedAtomIndicesByDim(matchesDEPT90, 0);
        spectrumDEPT135.setAssignedAtomIndicesByDim(matchesDEPT135, 0);
        
        return this.setImplicitHydrogenNumberFromDEPT(matchesDEPT90, matchesDEPT135, intensitiesDEPT135);
    } 
    

    /**
     *
     * @param matchesDEPT90
     * @param matchesDEPT135
     * @param intensitiesDEPT135
     */
    private int setImplicitHydrogenNumberFromDEPT(final ArrayList<Integer> matchesDEPT90, final ArrayList<Integer> matchesDEPT135, final ArrayList<Double> intensitiesDEPT135) {

        int matchDEPT90, matchDEPT135, hCount, hCountAll = 0;
        for (int i : this.atomTypeIndices.get("C")) {
            if ((this.mol.getAtom(i).getProperty(CDKConstants.NMRSHIFT_CARBON) != null) && (this.mol.getAtom(i).getImplicitHydrogenCount() == null)) {
                matchDEPT90 = matchesDEPT90.indexOf(i);
                matchDEPT135 = matchesDEPT135.indexOf(i);
                if (matchDEPT90 >= 0) {
                    // CH
                    hCount = 1;
                } else if (matchDEPT90 == -1 && matchDEPT135 >= 0) {
                    // CH2 or CH3
                    if (intensitiesDEPT135.get(matchDEPT135) < 0) {
                        hCount = 2;
                    } else if (intensitiesDEPT135.get(matchDEPT135) > 0) {
                        hCount = 3;
                    } else {
                        // qC
                        hCount = 0;
                    }
                } else {
                    // qC
                    hCount = 0;
                }
                this.mol.getAtom(i).setImplicitHydrogenCount(hCount);
                hCountAll += hCount;
                if (this.mol.getAtom(i).getProperty(ParseRawData.PROP_EQUIVALENCE) != null) {
                    for (Integer k : (ArrayList<Integer>) this.mol.getAtom(i).getProperty(ParseRawData.PROP_EQUIVALENCE)) {
                        this.mol.getAtom(k).setImplicitHydrogenCount(hCount);
                        hCountAll += hCount;
                    }
                }
            }
        }
        if(this.molFormula != null){
            System.out.println("assigned protons to carbons: " + hCountAll + " (" + MolecularFormulaManipulator.getElementCount(this.molFormula, "H") + ") -> " + (MolecularFormulaManipulator.getElementCount(this.molFormula, "H") - hCountAll) + " protons to be attached on hetero atoms!!!");
        } else {
            System.out.println("assigned protons to carbons: " + hCountAll+ "!!!");
        }
        
        return hCountAll;
    }
    
    
    /**
     * Creates a Spectrum class object from given HSQC input file in CSV format.
     *
     * @param pathToCSV path to HSQC peak table (Bruker's TopSpin csv file
     * format)
     * @param heavyAtomType Element name of H bonded heavy atom (e.g. "C") which also occurrs in
     * {@link Utils#getNMRShiftConstant(java.lang.String)}
     * @return 
     * @throws IOException
     */
    public final Spectrum parseHSQCviaCSV(final String pathToCSV, final String heavyAtomType) throws IOException {

        final Spectrum spectrum = new Spectrum( new String[]{Utils.getIsotopeIdentifier("H"),
                                                            Utils.getIsotopeIdentifier(heavyAtomType)},
                                                new ArrayList[]{casekit.NMR.Utils.parseCSV(pathToCSV, 5),
                                                                casekit.NMR.Utils.parseCSV(pathToCSV, 6)},
                                                casekit.NMR.Utils.parseCSV(pathToCSV, 9));
        
        return spectrum;
    }

    /**
     * Creates a Spectrum class object from given HSQC input file in XML format.
     *
     * @param pathToXML path to HSQC XML file
     * @param heavyAtomType Element name of H bonded heavy atom (e.g. "C") which also occurrs in
     * {@link Utils#getNMRShiftConstant(java.lang.String)}
     * @return 
     * @throws IOException
     * @throws javax.xml.parsers.ParserConfigurationException
     * @throws org.xml.sax.SAXException
     */
    public final Spectrum parseHSQCviaXML(final String pathToXML, final String heavyAtomType) throws IOException, ParserConfigurationException, SAXException {

        final Spectrum spectrum = new Spectrum( new String[]{   Utils.getIsotopeIdentifier("H"),
                                                                Utils.getIsotopeIdentifier(heavyAtomType)},
                                                new ArrayList[]{casekit.NMR.Utils.parseXML(pathToXML, 2, 2),
                                                                casekit.NMR.Utils.parseXML(pathToXML, 2, 1)},
                                                casekit.NMR.Utils.parseXML(pathToXML, 2, 3));
        
        return spectrum;
    }
    
    
     /**
     * Sets the proton shift(s) as list to belonging heavy atoms of an 
     * HSQC signal relationship.
     * The property is then set to {@link #CONST_PROP_PROTONSHIFTS} in 
     * {@link IAtomContainer#setProperty(java.lang.Object, java.lang.Object)}.
     * 
     * @param spectrum Spectrum class object consisting of Signal class objects 
     * where the proton values are given first and the heavy atom values as the second.
     * @param tolHeavyAtom tolerance value [ppm] for heavy atom shift matching
     */
    public final void setHSQC(final Spectrum spectrum, final double tolHeavyAtom) {
        
        final ArrayList<Double> shiftsHydrogen = spectrum.getShiftsByDim(0);
        final ArrayList<Double> shiftsHeavyAtom = spectrum.getShiftsByDim(1);
        ArrayList<Integer> matchesHeavyAtom = casekit.NMR.Utils.findShiftMatches(this.mol, shiftsHeavyAtom, tolHeavyAtom, Utils.getElementIdentifier(spectrum.getNuclei()[1]));
        matchesHeavyAtom = casekit.NMR.Utils.correctShiftMatches(this.mol, shiftsHeavyAtom, matchesHeavyAtom, tolHeavyAtom, Utils.getElementIdentifier(spectrum.getNuclei()[1]));
        
        spectrum.setAssignedAtomIndicesByDim(matchesHeavyAtom, 1);
        
        this.setImplicitHydrogenShifts(shiftsHydrogen, matchesHeavyAtom);
         
    }
    

    private void setImplicitHydrogenShifts(final ArrayList<Double> shiftsHydrogen, final ArrayList<Integer> matchesHeavyAtomType) {
        
        IAtom matchAtom;
        ArrayList<Double> assignedHydrogensShifts;
        for (int i = 0; i < matchesHeavyAtomType.size(); i++) {
            if (matchesHeavyAtomType.get(i) >= 0) {
                matchAtom = this.mol.getAtom(matchesHeavyAtomType.get(i));
                if (matchAtom.getImplicitHydrogenCount() == null || matchAtom.getImplicitHydrogenCount() == 0) {
                    continue;
                }
                if (matchAtom.getProperty(CDKConstants.NMRSPECTYPE_2D_HSQC) == null) {
                    matchAtom.setProperty(CDKConstants.NMRSPECTYPE_2D_HSQC, new ArrayList<>(matchAtom.getImplicitHydrogenCount()));
                }
                assignedHydrogensShifts = matchAtom.getProperty(CDKConstants.NMRSPECTYPE_2D_HSQC);
                if (assignedHydrogensShifts.size() < matchAtom.getImplicitHydrogenCount()) {
                    assignedHydrogensShifts.add(shiftsHydrogen.get(i));
                }
            }
        }
    }
    

    /**
     * Finds the matches with the lowest deviations between a given hydrogen
     * shift value set and implicit hydrogens of heavy atoms in the atom
     * container.
     *
     * @param shiftList shift value list to match
     * @param tol tolerance value [ppm]
     * @return
     */
    private ArrayList<Integer> findImplicitHydrogenShiftMatches(final ArrayList<Double> shiftList, final double tol) {

        final ArrayList<Integer> matches = new ArrayList<>();
        for (int i = 0; i < shiftList.size(); i++) {
            matches.add(this.findSingleImplicitHydrogenShiftMatch(shiftList.get(i), tol)[0]);
        }

        return matches;
    }

    /**
     * Finds a match with the lowest deviations between a given hydrogen
     * shift value and implicit hydrogens of heavy atoms in the atom
     * container.
     *
     * @param queryShift hydrogen shift value [ppm] to match
     * @param tol tolerance value [ppm] for matching
     * @return int array of two values: 1. index of matched heavy atom in
     * atom container, 2. index of matched hydrogen in hydrogen shift list
     * of corresponding found heavy atom
     */
    private int[] findSingleImplicitHydrogenShiftMatch(final double queryShift, final double tol) {

        int matchIndexAtom = -1;
        int matchIndexProton = -1;
        double minDiff = tol;
        ArrayList<Double> protonShiftList;
        for (int i = 0; i < this.mol.getAtomCount(); i++) {
            // skip atoms without implicit hydrogens
            if (this.mol.getAtom(i).getProperty(CDKConstants.NMRSPECTYPE_2D_HSQC) == null) {
                continue;
            }
            protonShiftList = this.mol.getAtom(i).getProperty(CDKConstants.NMRSPECTYPE_2D_HSQC);
            for (int j = 0; j < protonShiftList.size(); j++) {
                // figure out the atom with lowest shift deviation
                if ((queryShift - tol <= protonShiftList.get(j)) && (protonShiftList.get(j) <= queryShift + tol) && (Math.abs(queryShift - protonShiftList.get(j)) < minDiff)) {
                    minDiff = Math.abs(queryShift - protonShiftList.get(j));
                    matchIndexProton = j;
                    matchIndexAtom = i;
                }
            }
        }

        return new int[]{matchIndexAtom, matchIndexProton};
    }

    /**
     * Corrects a hydrogen match list regarding a given shift list and an atom
     * container.
     * This is useful when two ore more hydrogen shift values match
     * with the same hydrogen shift (actually heavy atom) in the atom container.
     * So the purpose here is to enable more unambiguous matches. This method
     * first looks for unambiguous matches and calculates the median of the
     * difference values between the shift list values and the shifts of atom
     * container. Then, all shift list values are adjusted (+/-) with this
     * median value.
     *
     * @param shifts Shift value list to match
     * @param matches Match list to correct
     * @param tol Tolerance value [ppm] for hydrogen rematching
     * @return
     */
    private ArrayList<Integer> correctHydrogenShiftMatches(final ArrayList<Double> shifts, ArrayList<Integer> matches, final double tol) {

        int matchIndex, middle;
        double diff, median;
        int[] singleMatchIndex;
        ArrayList<Double> singleMatchShifts;
        ArrayList<Double> diffs = new ArrayList<>();
        final HashSet<Integer> uniqueMatchIndicesSet = new HashSet<>(matches);
        for (Integer matchIndexAtomContainer : uniqueMatchIndicesSet) {
            if (Collections.frequency(matches, matchIndexAtomContainer) == 1) {
                matchIndex = matches.indexOf(matchIndexAtomContainer);
                if (matches.get(matchIndex) >= 0) {
                    singleMatchIndex = this.findSingleImplicitHydrogenShiftMatch(shifts.get(matchIndex), tol);
                    singleMatchShifts = this.mol.getAtom(singleMatchIndex[0]).getProperty(CDKConstants.NMRSPECTYPE_2D_HSQC);
                    diff = shifts.get(matchIndex) - singleMatchShifts.get(singleMatchIndex[1]);
                    diffs.add(diff);
                }
            }
        }
        if (diffs.size() > 0) {
            middle = diffs.size() / 2;
            if (diffs.size() % 2 == 1) {
                median = diffs.get(middle);
            } else {
                median = (diffs.get(middle - 1) + diffs.get(middle)) / 2.0;
            }
            // add or subtract the median of the differences to all shift list values (input) and match again then
            for (int i = 0; i < shifts.size(); i++) {
                shifts.set(i, shifts.get(i) - median);
            }
            // rematch
            matches = this.findImplicitHydrogenShiftMatches(shifts, tol);
        }

        return matches;
    }
    

    /**
     * Creates a Spectrum class object from given H,H-COSY input file in CSV format.
     *
     * @param pathToCSV path to H,H-COSY peak table (Bruker's TopSpin csv
     * file
     * format)
     * @return
     * @throws IOException
     */
    public final Spectrum parseHHCOSYviaCSV(final String pathToCSV) throws IOException {

        final Spectrum spectrum = new Spectrum( new String[]{   Utils.getIsotopeIdentifier("H"), 
                                                                Utils.getIsotopeIdentifier("H")}, 
                                                new ArrayList[]{casekit.NMR.Utils.parseCSV(pathToCSV, 5), 
                                                                casekit.NMR.Utils.parseCSV(pathToCSV, 6)},
                                                casekit.NMR.Utils.parseCSV(pathToCSV, 9));

        return spectrum;
    }

    /**
     * Creates a Spectrum class object from given H,H-COSY input file in XML format.
     *
     * @param pathToXML path to H,H-COSY peak XML file (Bruker's TopSpin XML
     * file format)
     * @return
     * @throws IOException
     * @throws javax.xml.parsers.ParserConfigurationException
     * @throws org.xml.sax.SAXException
     */
    public final Spectrum parseHHCOSYviaXML(final String pathToXML) throws IOException, ParserConfigurationException, SAXException {

        final Spectrum spectrum = new Spectrum( new String[]{   Utils.getIsotopeIdentifier("H"), 
                                                                Utils.getIsotopeIdentifier("H")}, 
                                                new ArrayList[]{casekit.NMR.Utils.parseXML(pathToXML, 2, 2), 
                                                                casekit.NMR.Utils.parseXML(pathToXML, 2, 1)},
                                                casekit.NMR.Utils.parseXML(pathToXML, 2, 3));

        return spectrum;
    }

    /**
     * Sets links between two heavy atoms of H,H-COSY signals. The property
     * is then set to {@link #CONST_PROP_HHCOSY} in 
     * {@link IAtomContainer#setProperty(java.lang.Object, java.lang.Object)}
     *
     * @param spectrum Spectrum class object containing the 2D spectrum proton shift information
     * @param tol tolerance value [ppm] for matching belonging protons 
     * of heavy atom 
     * @return true if the links could be set; otherwise false
     */
    public final boolean setHHCOSY(final Spectrum spectrum, final double tol) {

        final ArrayList<Integer> hydrogenShiftMatches1 = this.findImplicitHydrogenShiftMatches(spectrum.getShiftsByDim(0), tol);
        final ArrayList<Integer> hydrogenShiftMatches2 = this.findImplicitHydrogenShiftMatches(spectrum.getShiftsByDim(1), tol);
        // are all signals bidirectional?
        if (!casekit.NMR.Utils.isBidirectional(hydrogenShiftMatches1, hydrogenShiftMatches2)) {
            return false;
        }
        casekit.NMR.Utils.setBidirectionalLinks(this.mol, hydrogenShiftMatches1, hydrogenShiftMatches2, CDKConstants.NMRSPECTYPE_2D_HHCOSY);

        spectrum.setAssignedAtomIndicesByDim(hydrogenShiftMatches1, 0);
        spectrum.setAssignedAtomIndicesByDim(hydrogenShiftMatches2, 1);
        
        return true;
    }

    
    /**
     * Creates a Spectrum class object from given INADEQUATE input file in CSV format.
     *
     * @param pathToCSV path to INADEQUATE peak table (Bruker's TopSpin csv
     * file format)
     * @return
     * @throws IOException
     */
    public final Spectrum parseINADEQUATEviaCSV(final String pathToCSV) throws IOException {

        final Spectrum spectrum = new Spectrum( new String[]{   Utils.getIsotopeIdentifier("C"), 
                                                                Utils.getIsotopeIdentifier("C")}, 
                                                new ArrayList[]{casekit.NMR.Utils.parseCSV(pathToCSV, 5), 
                                                                casekit.NMR.Utils.parseCSV(pathToCSV, 6)},
                                                casekit.NMR.Utils.parseCSV(pathToCSV, 9));
        
        return spectrum;
    }

    /**
     * Creates a Spectrum class object from given INADEQUATE input file in XML format.
     *
     * @param pathToXML path to INADEQUATE peak XML file (Bruker's TopSpin XML
     * file format)
     * @return
     * @throws IOException
     * @throws javax.xml.parsers.ParserConfigurationException
     * @throws org.xml.sax.SAXException
     */
    public final Spectrum parseINADEQUATEviaXML(final String pathToXML) throws IOException, ParserConfigurationException, SAXException {

        final Spectrum spectrum = new Spectrum( new String[]{   Utils.getIsotopeIdentifier("C"), 
                                                                Utils.getIsotopeIdentifier("C")}, 
                                                new ArrayList[]{casekit.NMR.Utils.parseXML(pathToXML, 2, 2), 
                                                                casekit.NMR.Utils.parseXML(pathToXML, 2, 1)},
                                                casekit.NMR.Utils.parseXML(pathToXML, 2, 3));
        
        return spectrum;
    }

    
    /**
     * Sets links between two carbon atoms in an INADEQUATE signal relationship.
     * The property is then set to {@link #CONST_PROP_INADEQUATE} in 
     * {@link IAtomContainer#setProperty(java.lang.Object, java.lang.Object)}.
     * Returns true if all signals are bidirectional, so that atom A has a
     * signal according to atom B and vice versa.
     * 
     * @param spectrum Spectrum class object consisting of Signal class objects 
     * @param tol tolerance value [ppm] for carbon atom shift matching
     * @return 
     */
    public final boolean setINADEQUATE(final Spectrum spectrum, final double tol) {

        final ArrayList<Integer> carbonShiftMatches1 = casekit.NMR.Utils.findShiftMatches(this.mol, spectrum.getShiftsByDim(0), tol, "C");
        final ArrayList<Integer> carbonShiftMatches2 = casekit.NMR.Utils.findShiftMatches(this.mol, spectrum.getShiftsByDim(1), tol, "C");
        // are all signals bidirectional?
        if (!casekit.NMR.Utils.isBidirectional(carbonShiftMatches1, carbonShiftMatches2)) {
            return false;
        }
        casekit.NMR.Utils.setBidirectionalLinks(this.mol, carbonShiftMatches1, carbonShiftMatches2, CDKConstants.NMRSPECTYPE_2D_INADEQUATE);

        spectrum.setAssignedAtomIndicesByDim(carbonShiftMatches1, 0);
        spectrum.setAssignedAtomIndicesByDim(carbonShiftMatches2, 1);
        
        return true;
    }
    
    
    /**
     * Creates a Spectrum class object from given HMBC input file in CSV format.
     *
     * @param pathToCSV path to HMBC peak table (Bruker's TopSpin csv
     * file format)
     * @param heavyAtomType Element name (e.g. "C") which also occurrs in
     * {@link Utils#getNMRShiftConstant(java.lang.String)}
     * @return 
     * @throws IOException
     */
    public final Spectrum parseHMBCviaCSV(final String pathToCSV, final String heavyAtomType) throws IOException {
        
        final Spectrum spectrum = new Spectrum( new String[]{Utils.getIsotopeIdentifier("H"),
                                                            Utils.getIsotopeIdentifier(heavyAtomType)},
                                                new ArrayList[]{casekit.NMR.Utils.parseCSV(pathToCSV, 5),
                                                                casekit.NMR.Utils.parseCSV(pathToCSV, 6)},
                                                casekit.NMR.Utils.parseCSV(pathToCSV, 9));
        
        return spectrum;
    }

    
    /**
     * Creates a Spectrum class object from given HMBC input file in XML format.
     * 
     * @param pathToXML path to HMBC peak XML file (Bruker's TopSpin XML file
     * format)
     * @param heavyAtomType Element name (e.g. "C") which also occurrs in
     * {@link Utils#getNMRShiftConstant(java.lang.String)}
     * @return 
     * @throws IOException
     * @throws javax.xml.parsers.ParserConfigurationException
     * @throws org.xml.sax.SAXException
     */
    public final Spectrum parseHMBCviaXML(final String pathToXML, final String heavyAtomType) throws IOException, ParserConfigurationException, SAXException {

        final Spectrum spectrum = new Spectrum( new String[]{   Utils.getIsotopeIdentifier("H"), 
                                                                Utils.getIsotopeIdentifier(heavyAtomType)}, 
                                                new ArrayList[]{casekit.NMR.Utils.parseXML(pathToXML, 2, 2), 
                                                                casekit.NMR.Utils.parseXML(pathToXML, 2, 1)},
                                                casekit.NMR.Utils.parseXML(pathToXML, 2, 3));
        
        return spectrum;
    }
    
    
    /**
     * Sets links between heavy atoms which are in HMBC signal relationship.
     * The property is then set to {@link #CONST_PROP_HMBC} in 
     * {@link IAtomContainer#setProperty(java.lang.Object, java.lang.Object)}.
     * 
     * @param spectrum Spectrum class object consisting of Signal class objects 
     * where the proton shift values is given first and the heavy atom shifts as the second.
     * @param tolHydrogen tolerance value [ppm] for hydrogen shift matching
     * @param tolHeavy tolerance value [ppm] for heavy atom shift matching
     */
    public final void setHMBC(final Spectrum spectrum, final double tolHydrogen, final double tolHeavy) {
        
        final ArrayList<Double> shiftsHydrogen = spectrum.getShiftsByDim(0);
        final ArrayList<Integer> matchesHydrogen = this.correctHydrogenShiftMatches(shiftsHydrogen, this.findImplicitHydrogenShiftMatches(shiftsHydrogen, tolHydrogen), tolHydrogen);        
        final ArrayList<Double> shiftsHeavyAtom = spectrum.getShiftsByDim(1);
        ArrayList<Integer> matchesHeavyAtom = casekit.NMR.Utils.findShiftMatches(this.mol, shiftsHeavyAtom, tolHeavy, Utils.getElementIdentifier(spectrum.getNuclei()[1]));
        matchesHeavyAtom = casekit.NMR.Utils.correctShiftMatches(this.mol, shiftsHeavyAtom, matchesHeavyAtom, tolHeavy, Utils.getElementIdentifier(spectrum.getNuclei()[1]));
        
        ArrayList<Integer> HMBCList;
        for (int i = 0; i < matchesHydrogen.size(); i++) {
            if (matchesHydrogen.get(i) >= 0 && matchesHeavyAtom.get(i) >= 0) {
                if (this.mol.getAtom(matchesHydrogen.get(i)).getProperty(CDKConstants.NMRSPECTYPE_2D_HMBC) == null) {
                    this.mol.getAtom(matchesHydrogen.get(i)).setProperty(CDKConstants.NMRSPECTYPE_2D_HMBC, new ArrayList<>());
                }
                HMBCList = this.mol.getAtom(matchesHydrogen.get(i)).getProperty(CDKConstants.NMRSPECTYPE_2D_HMBC);
                if (!HMBCList.contains(matchesHeavyAtom.get(i))) {
                    HMBCList.add(matchesHeavyAtom.get(i));
                }
            }
        }
        
        spectrum.setAssignedAtomIndicesByDim(matchesHydrogen, 0);
        spectrum.setAssignedAtomIndicesByDim(matchesHeavyAtom, 1);
    }
}
