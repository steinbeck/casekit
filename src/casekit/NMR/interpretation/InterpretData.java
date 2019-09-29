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
package casekit.NMR.interpretation;

import casekit.NMR.Utils;
import casekit.NMR.model.Assignment;
import casekit.NMR.model.Spectrum;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import org.openscience.cdk.silent.Atom;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class InterpretData {
    
    final private IAtomContainer mol;
    final private IMolecularFormula molFormula;
    private HashMap<String, ArrayList<Integer>> atomTypeIndices;
    final private HashMap<String, Spectrum> spectra = new HashMap<>();
    final private HashMap<String, Assignment> assignments = new HashMap<>();

    /**
     * Creates an instances of this class with an empty class atom container.
     */
    public InterpretData(){
        this.molFormula = null;
        this.mol = SilentChemObjectBuilder.getInstance().newAtomContainer();
        this.updateAtomTypeIndices();
    }
    
    /**
     * Creates an instances of this class with a class atom container consisting 
     * of all heavy atoms in given molecular formula.
     * 
     * @param molFormula IMolecularFormula object for IAtomContainer creation
     */
    public InterpretData(final IMolecularFormula molFormula){
        this.molFormula = molFormula;
        this.mol = Utils.removeAtoms(MolecularFormulaManipulator.getAtomContainer(this.molFormula), "H");
        this.updateAtomTypeIndices();
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
    private void updateAtomTypeIndices(){
        
        this.atomTypeIndices = Utils.getAtomTypeIndices(this.mol);
    }
    
    /**
     * Returns all given and used spectra.
     *
     * @return
     */
    public final HashMap<String, Spectrum> getSpectra(){
        
        return this.spectra;
    } 
    
    
    /**
     * Returns all created and used Assignment objects. The assigned indices 
     * refer to atom indices in class atom container.
     *
     * @return
     */
    public final HashMap<String, Assignment> getAssignments(){
        
        return this.assignments;
    } 
    
    
    /**
     * Returns one specific created and used Assignment object. 
     * The assigned indices refer to atom indices in class atom container.
     *
     * @param spectrum
     * @return
     */
    public final Assignment getAssignment(final Spectrum spectrum){
        
        if (spectrum.getSpecType().equals(CDKConstants.NMRSPECTYPE_1D_DEPT90) || spectrum.getSpecType().equals(CDKConstants.NMRSPECTYPE_1D_DEPT135)) {
            
            return this.getAssignments().get(spectrum.getSpecType());
        }
         
        return this.assignments.get(spectrum.getSpecType() + "_" + Utils.getSpectrumNucleiAsString(spectrum));
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
     * equivalent signal selection by user.
     * 
     * @param spectrum Spectrum class object containing the 1D shift information
     * @throws java.io.IOException
     * @throws org.openscience.cdk.exception.CDKException
     */
    public final void assign1DSpectrum(final Spectrum spectrum) throws Exception {
        // checks whether number of signals is equal to molecular formula if given
        // if not equal then edit signal list in spectrum
        this.check1DSpectrum(spectrum);
        // assign shift values to atoms sequentially
        this.assignShiftValuesToAtoms(spectrum);
        
        final Assignment assignment = new Assignment(spectrum);
        if(this.atomTypeIndices.get(Utils.getAtomTypeFromSpectrum(spectrum, 0)) != null){
            assignment.setAssignments(0, this.atomTypeIndices.get(Utils.getAtomTypeFromSpectrum(spectrum, 0)));
        }        
        
        this.spectra.put(CDKConstants.NMRSPECTYPE_1D + "_" + Utils.getSpectrumNucleiAsString(spectrum), spectrum);
        this.assignments.put(CDKConstants.NMRSPECTYPE_1D + "_" + Utils.getSpectrumNucleiAsString(spectrum), assignment);
    }
    
    /**
     * Checks the number of signals in a spectrum against the number of atoms
     * in molecular formula of class, if given. In case of different numbers, 
     * a user input for spectrum editing will be requested. 
     *
     * @param spectrum
     * @throws IOException
     * @see Utils#editSignalsInSpectrum(Spectrum, IMolecularFormula, int)
     */
    private void check1DSpectrum(final Spectrum spectrum) throws Exception {
        if(this.molFormula != null) {
            final int diff = Utils.getDifferenceSpectrumSizeAndMolecularFormulaCount(spectrum, this.molFormula, 0);
            if (diff != 0) {
                // adjust Spectrum size by user
                Utils.editSignalsInSpectrum(spectrum, this.molFormula, 0);
            }
        }
    }
    
    
    /**
     * Sets shift values in atoms of class atom container as property (see below), sequentially.
     *
     * @param spectrum Spectrum class object which contains shifts in first 
     * dimension
     * @see  Utils#getNMRShiftConstant(java.lang.String)
     *
     */
    private void assignShiftValuesToAtoms(final Spectrum spectrum){
        final String atomType = Utils.getAtomTypeFromSpectrum(spectrum, 0);
        final ArrayList<Double> shifts = spectrum.getShifts(0);
        if((this.molFormula == null) && !atomType.equals("H")){
            // (re-)filling up of peaks for that atom type from given peak list in spectrum
            this.removeAtoms(atomType); 
            IAtom atom;
            for (final double shift : shifts) {                
                atom = new Atom(atomType);
                atom.setProperty(Utils.getNMRShiftConstant(atomType), shift);
                atom.setImplicitHydrogenCount(null);
                this.mol.addAtom(atom);                
            }
            this.updateAtomTypeIndices();
        }
        // assign shifts to atoms as property
        if(this.atomTypeIndices.get(atomType) != null){
            int assignedShiftCount = 0;
            for (final int i : this.atomTypeIndices.get(atomType)) {
                if (assignedShiftCount < shifts.size()) {
                    // shift assignment in atom
                    this.mol.getAtom(i).setProperty(Utils.getNMRShiftConstant(atomType), shifts.get(assignedShiftCount));
                }
                assignedShiftCount++;
            }
        }
    }
    
    
    /**
     * Removes atoms from a given atom type from the class' atom container.
     *
     * @param atomType Atom type (element's name, e.g. C or Br)
     * @return IAtomContainer where the atoms were removed
     */
    private void removeAtoms(final String atomType) {
        if(this.getAtomTypeIndices().get(atomType) == null){
            return;
        }
        final ArrayList<IAtom> toRemoveList = new ArrayList<>();
        for (final int i: this.getAtomTypeIndices().get(atomType)) {
            toRemoveList.add(this.mol.getAtom(i));
        }
        for (IAtom iAtom : toRemoveList) {
            this.mol.removeAtom(iAtom);
        }
        
        this.updateAtomTypeIndices();
    }

    /**
     * Sets the assignments of carbon atoms in class atom container
     * by usage of DEPT90 and DEPT135 information. The implicit hydrogen count
     * property is set too.
     * 
     * @see InterpretData#setImplicitHydrogenCountsFromDEPT()
     * 
     * @param spectrum1D_DEPT90 DEPT90 spectrum
     * @param spectrum1D_DEPT135 DEPT135 spectrum which has to contain intensity 
     * information
     * @param tol tolerance value [ppm] for carbon shift matching
     * @return false if 1-dimensional 13C spectrum is missing (not set beforehand)
     * or something is missing in one of the two input spectra
     * 
     */
    public final boolean assignDEPT(final Spectrum spectrum1D_DEPT90, final Spectrum spectrum1D_DEPT135, final double tol){        
        if((spectrum1D_DEPT90 == null) || (spectrum1D_DEPT135 == null) || (spectrum1D_DEPT135.getIntensities() == null) 
                || (this.getSpectra().get(CDKConstants.NMRSPECTYPE_1D + "_13C") == null)){
            return false;
        }        

        final Assignment assignment1D_DEPT90 = new Assignment(spectrum1D_DEPT90);
        final Assignment assignment1D_DEPT135 = new Assignment(spectrum1D_DEPT135);
        final ArrayList<Integer> matchesIn1DSpectrum_DEPT90 = this.findMatchesIn1DSpectra(spectrum1D_DEPT90, 0, tol);
        final ArrayList<Integer> matchesIn1DSpectrum_DEPT135 = this.findMatchesIn1DSpectra(spectrum1D_DEPT135, 0, tol);
        final Assignment assignment1D_13C = this.getAssignment(this.getSpectra().get(CDKConstants.NMRSPECTYPE_1D + "_13C"));
               
        for (int i = 0; i < assignment1D_DEPT90.getAssignmentsCount(); i++) {            
            if (assignment1D_13C.getAssignment(0, matchesIn1DSpectrum_DEPT90.get(i)) >= 0) {
                assignment1D_DEPT90.setAssignment(0, i, assignment1D_13C.getAssignment(0, matchesIn1DSpectrum_DEPT90.get(i)));
            }
        }
        for (int i = 0; i < assignment1D_DEPT135.getAssignmentsCount(); i++) {            
            if (assignment1D_13C.getAssignment(0, matchesIn1DSpectrum_DEPT135.get(i)) >= 0) {
                assignment1D_DEPT135.setAssignment(0, i, assignment1D_13C.getAssignment(0, matchesIn1DSpectrum_DEPT135.get(i)));
            }
        }
        
        this.spectra.put(CDKConstants.NMRSPECTYPE_1D_DEPT90, spectrum1D_DEPT90);
        this.assignments.put(CDKConstants.NMRSPECTYPE_1D_DEPT90, assignment1D_DEPT90);
        this.spectra.put(CDKConstants.NMRSPECTYPE_1D_DEPT135, spectrum1D_DEPT135);
        this.assignments.put(CDKConstants.NMRSPECTYPE_1D_DEPT135, assignment1D_DEPT135);
        
        this.setImplicitHydrogenCountsFromDEPT();
        
        return true;
    } 
    
    
    /**
     * Sets the implicitHydrogenCount() property in atoms of class atom container 
     * by using the already set DEPT information.
     * @see InterpretData#assignDEPT(casekit.NMR.model.Spectrum, casekit.NMR.model.Spectrum, double)
     */
    private void setImplicitHydrogenCountsFromDEPT() {
        
        final ArrayList<Double> intensitiesDEPT135 = this.getSpectra().get(CDKConstants.NMRSPECTYPE_1D_DEPT135).getIntensities();
        final ArrayList<Integer> matchesDEPT90InAtomContainer = this.getAssignedAtomIndices(this.getSpectra().get(CDKConstants.NMRSPECTYPE_1D_DEPT90), 0);
        final ArrayList<Integer> matchesDEPT135InAtomContainer = this.getAssignedAtomIndices(this.getSpectra().get(CDKConstants.NMRSPECTYPE_1D_DEPT135), 0);
        
        int matchDEPT90, matchDEPT135, hCount, hCountAll = 0;
        for (int i : this.atomTypeIndices.get("C")) {
            if ((this.mol.getAtom(i).getProperty(CDKConstants.NMRSHIFT_CARBON) != null) && (this.mol.getAtom(i).getImplicitHydrogenCount() == null)) {
                matchDEPT90 = matchesDEPT90InAtomContainer.indexOf(i);
                matchDEPT135 = matchesDEPT135InAtomContainer.indexOf(i);
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
                if( this.mol.getAtom(i).getImplicitHydrogenCount() >= 3){
                    this.mol.getAtom(i).setHybridization(IAtomType.Hybridization.SP3);
                }
                hCountAll += hCount;
            }
        }
        if(this.molFormula != null){
            System.out.println("assigned protons to carbons: " + hCountAll + " (" + MolecularFormulaManipulator.getElementCount(this.molFormula, "H") + ") -> " + (MolecularFormulaManipulator.getElementCount(this.molFormula, "H") - hCountAll) + " protons to be attached on hetero atoms!!!");
        } else {
            System.out.println("assigned protons to carbons: " + hCountAll+ "!!!");
        }
        
    }


    /**
     *  
     * @param spectrum Spectrum class object consisting of Signal class objects 
     * where the proton shifts values are given in first dimension and the 
     * heavy atom shifts in the second.
     * @param tolProton tolerance value [ppm] for proton shift matching
     * @param tolHeavyAtom tolerance value [ppm] for heavy atom shift matching
     */
    public final void assignHSQC(final Spectrum spectrum, final double tolProton, final double tolHeavyAtom) {
        // assign index of matching atoms to both dimensions and save the Spectrum and Assignment objects in class
        this.assign2DSpectrum(spectrum, tolProton, tolHeavyAtom);
        // in case the 1H spectrum is given, then assign protons to same indices from belonging carbon atoms
        if(this.getAssignments().get(CDKConstants.NMRSPECTYPE_1D + "_1H") != null){
            final Assignment assignment1D_1H = this.getAssignments().get(CDKConstants.NMRSPECTYPE_1D + "_1H");
            final Assignment assignment2D_HSQC = this.getAssignments().get(spectrum.getSpecType() + "_" + Utils.getSpectrumNucleiAsString(spectrum));
            final ArrayList<Integer> matchesIn1DSpectrum_1H = this.findMatchesIn1DSpectra(spectrum, 0, tolProton);
            
            for (int i = 0; i < assignment2D_HSQC.getAssignmentsCount(); i++) {
                // if heavy atom i has an assignment in class atom container then assign that index i to belonging protons as index
                if (assignment2D_HSQC.getAssignment(1, i) >= 0) {
                    assignment1D_1H.setAssignment(0, matchesIn1DSpectrum_1H.get(i), assignment2D_HSQC.getAssignment(1, i));
                    assignment2D_HSQC.setAssignment(0, i, assignment1D_1H.getAssignment(0, matchesIn1DSpectrum_1H.get(i)));
                }
            }
        }
        // attach protons on other heavy atoms than carbons via HSQC assignment counting
        if(!spectrum.getNuclei()[1].equals("13C")){
            final Assignment assignment2D_HSQC = this.getAssignments().get(spectrum.getSpecType() + "_" + Utils.getSpectrumNucleiAsString(spectrum));
            for (int i = 0; i < assignment2D_HSQC.getAssignmentsCount(); i++) {
                if((assignment2D_HSQC.getAssignment(1, i) > -1)){
                    if(this.mol.getAtom(assignment2D_HSQC.getAssignment(1, i)).getImplicitHydrogenCount() == null){
                        this.mol.getAtom(assignment2D_HSQC.getAssignment(1, i)).setImplicitHydrogenCount(0);
                    }
                    this.mol.getAtom(assignment2D_HSQC.getAssignment(1, i)).setImplicitHydrogenCount(this.mol.getAtom(assignment2D_HSQC.getAssignment(1, i)).getImplicitHydrogenCount() + 1);
                }
            }
        }
    }


    private void assign2DSpectrum(final Spectrum spectrum, final double tolDim1, final double tolDim2){

        final ArrayList<Integer> matchesQueryIn1DSpectrumDim1 = this.findMatchesIn1DSpectra(spectrum, 0, tolDim1);
        final ArrayList<Integer> matchesQueryIn1DSpectrumDim2 = this.findMatchesIn1DSpectra(spectrum, 1, tolDim2);
        final ArrayList<Integer> matches1DInAtomContainerDim1 = this.getAssignedAtomIndices(this.getSpectra().get(CDKConstants.NMRSPECTYPE_1D + "_" + spectrum.getNuclei()[0]), 0);
        final ArrayList<Integer> matches1DInAtomContainerDim2 = this.getAssignedAtomIndices(this.getSpectra().get(CDKConstants.NMRSPECTYPE_1D + "_" + spectrum.getNuclei()[1]), 0);

        final Assignment assignment = new Assignment(spectrum);
        for (int i = 0; i < matchesQueryIn1DSpectrumDim1.size(); i++) {
            if((matches1DInAtomContainerDim1 != null) && (matchesQueryIn1DSpectrumDim1.get(i) >= 0)){
                assignment.setAssignment(0, i, matches1DInAtomContainerDim1.get(matchesQueryIn1DSpectrumDim1.get(i)));
            }
            if((matches1DInAtomContainerDim2 != null) && (matchesQueryIn1DSpectrumDim2.get(i) >= 0)){
                assignment.setAssignment(1, i, matches1DInAtomContainerDim2.get(matchesQueryIn1DSpectrumDim2.get(i)));
            }
        }

        this.spectra.put(spectrum.getSpecType() + "_" + Utils.getSpectrumNucleiAsString(spectrum), spectrum);
        this.assignments.put(spectrum.getSpecType() + "_" + Utils.getSpectrumNucleiAsString(spectrum), assignment);
    }


    private ArrayList<Integer> findMatchesIn1DSpectra(final Spectrum spectrum, final int dim, final double tol){

        ArrayList<Integer> matchesQueryInOrigin1DSpectrum = new ArrayList<>();
//        final ArrayList<Double> shiftsQuery = spectrum.getShifts(dim);
//        if(this.getSpectra().get(CDKConstants.NMRSPECTYPE_1D + "_" + spectrum.getNuclei()[dim]) != null){
//            final ArrayList<Double> shiftsOrigin1DSpectrum = this.getSpectra().get(CDKConstants.NMRSPECTYPE_1D + "_" + spectrum.getNuclei()[dim]).getShifts(0);
//            matchesQueryInOrigin1DSpectrum = Utils.findShiftMatches(shiftsOrigin1DSpectrum, shiftsQuery, tol);
//            matchesQueryInOrigin1DSpectrum = Utils.correctShiftMatches(shiftsOrigin1DSpectrum, shiftsQuery, matchesQueryInOrigin1DSpectrum, tol);
//        } else {
//            for (int i = 0; i < spectrum.getSignalCount(); i++) {
//                matchesQueryInOrigin1DSpectrum.add(-1);
//            }
//        }

        return matchesQueryInOrigin1DSpectrum;
    }

    /**
     * Returns the indices of atoms within the class atom container which match
     * to the shifts of given spectrum and dimension.
     *
     * @param spectrum
     * @param dim
     * @return
     */
    public final ArrayList<Integer> getAssignedAtomIndices(final Spectrum spectrum, final int dim){
        
        if(spectrum == null){
            return null;
        } else if(this.getAssignment(spectrum) == null){    
            final ArrayList<Integer> atomIndices = new ArrayList<>();
            for (int i = 0; i < spectrum.getSignalCount(); i++) {
                atomIndices.add(-1);
            }
            return atomIndices;
        }
        
        return new ArrayList<>(this.getAssignment(spectrum).getAssignments(dim));
    }


    /**
     * Sets links between two heavy atoms of H,H-COSY signals.
     *
     * @param spectrum Spectrum class object containing the 2D spectrum proton shift information
     * @param tolProton tolerance value [ppm] for matching belonging protons 
     * of heavy atom 
     * @return 
     */
    public final boolean assignHHCOSY(final Spectrum spectrum, final double tolProton) {
        
        final ArrayList<Integer> protonShiftMatches1 = this.findMatchesIn1DSpectra(spectrum, 0, tolProton);        
        final ArrayList<Integer> protonShiftMatches2 = this.findMatchesIn1DSpectra(spectrum, 1, tolProton);
        // are all signals bidirectional?
        if (!Utils.isBidirectional(protonShiftMatches1, protonShiftMatches2)) {
            return false;
        }
        this.assign2DSpectrum(spectrum, tolProton, tolProton);
        
        return true;
    }


    /**
     * Sets links between two carbon atoms in an INADEQUATE signal relationship.
     * Returns true if all signals are bidirectional, so that atom A has a
     * signal according to atom B and vice versa.
     * 
     * @param spectrum Spectrum class object consisting of Signal class objects 
     * @param tolCarbon tolerance value [ppm] for carbon atom shift matching
     * @return 
     */
    public final boolean assignINADEQUATE(final Spectrum spectrum, final double tolCarbon) {

        final ArrayList<Integer> carbonShiftMatches1 = this.findMatchesIn1DSpectra(spectrum, 0, tolCarbon);
        final ArrayList<Integer> carbonShiftMatches2 = this.findMatchesIn1DSpectra(spectrum, 1, tolCarbon);        
        // are all signals bidirectional?
        if (!casekit.NMR.Utils.isBidirectional(carbonShiftMatches1, carbonShiftMatches2)) {
            return false;
        }
        this.assign2DSpectrum(spectrum, tolCarbon, tolCarbon);
        
        final ArrayList<Integer> indicesInAtomContainerDim1 = this.getAssignedAtomIndices(spectrum, 0);
        final ArrayList<Integer> indicesInAtomContainerDim2 = this.getAssignedAtomIndices(spectrum, 1);
        for (int i = 0; i < spectrum.getSignalCount(); i++) {
            if((indicesInAtomContainerDim1.get(i) > -1) && (indicesInAtomContainerDim2.get(i) > -1)){                    
                this.setBond(indicesInAtomContainerDim1.get(i), indicesInAtomContainerDim2.get(i));
            }                
        }
        
        return true;
    }
    
    
    private void setBond(final int index1, final int index2) {

        if (this.mol.getBond(this.mol.getAtom(index1), this.mol.getAtom(index2)) != null) {
            this.mol.removeBond(this.mol.getAtom(index1), this.mol.getAtom(index2));
        }
        this.mol.addBond(index1, index2, IBond.Order.UNSET);
    }


    /**
     * Sets links between heavy atoms which are in HMBC signal relationship.
     *
     * @param spectrum Spectrum class object consisting of Signal class objects 
     * where the proton shift values is given first and the heavy atom shifts as the second.
     * @param tolProton tolerance value [ppm] for hydrogen shift matching
     * @param tolHeavy tolerance value [ppm] for heavy atom shift matching
     */
    public final void assignHMBC(final Spectrum spectrum, final double tolProton, final double tolHeavy) {
     
       this.assign2DSpectrum(spectrum, tolProton, tolHeavy);
    }
}
