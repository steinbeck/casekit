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

package casekit.NMR.dbservice;

import casekit.NMR.Utils;
import casekit.NMR.model.Assignment;
import casekit.NMR.model.Signal;
import casekit.NMR.model.Spectrum;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.AtomContainerSet;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.*;

public class NMRShiftDB {
    /**
     * Returns the molecules of a given MOL/SDF file.
     * This function sets the molecule aromaticity (with allowed exocyclic pi
     * bonds) by using the
     * {@link Utils#setAromaticitiesInAtomContainer(org.openscience.cdk.interfaces.IAtomContainer)}
     * function.
     *
     * @param pathToNMRShiftDB path to NMRShiftDB file
     * @param setAromaticity whether to set aromaticities in structures or not
     * @return
     * @throws FileNotFoundException
     * @throws CDKException
     * @deprecated
     */
    public static IAtomContainerSet getStructuresFromSDFile(final String pathToNMRShiftDB, final boolean setAromaticity) throws FileNotFoundException, CDKException {
        final IAtomContainerSet acSet = new AtomContainerSet();
        final IteratingSDFReader iterator = new IteratingSDFReader(
                new FileReader(pathToNMRShiftDB),
                SilentChemObjectBuilder.getInstance()
        );
        IAtomContainer ac;
        while (iterator.hasNext()) {
            ac = iterator.next();
            if(setAromaticity){
                Utils.setAromaticitiesInAtomContainer(ac);
            }
            acSet.addAtomContainer(ac);
        }

        return acSet;
    }

    /**
     * Returns all spectra for each molecule and a given nucleus which exist as
     * property in a NMRSHiftDB SDF.
     *
     * @param pathToNMRShiftDB path to NMRShiftDB file
     * @param nucleus nucleus of requested spectra
     * @return
     * @throws FileNotFoundException
     * @throws CDKException
     *
     */
    public static ArrayList<ArrayList<Object[]>> getSpectraFromNMRShiftDB(final String pathToNMRShiftDB, final String nucleus) throws FileNotFoundException, CDKException {
        final ArrayList<ArrayList<Object[]>> spectraSet = new ArrayList<>();
        final IteratingSDFReader iterator = new IteratingSDFReader(
                new FileReader(pathToNMRShiftDB),
                SilentChemObjectBuilder.getInstance()
        );
        IAtomContainer ac;
        Spectrum spectrum;
        ArrayList<Object[]> spectra;
        HashMap<String, String> spectraStrings;
        String spectrumIndexInRecord, solvent;
        while (iterator.hasNext()) {
            ac = iterator.next();
            if(ac == null){
                continue;
            }
            spectraStrings = getSpectraStrings(ac, nucleus);
            if(spectraStrings.isEmpty() || (ac.getProperty("Solvent") == null)){
                continue;
            }
            spectra = new ArrayList<>();
            for (final String spectrumPropertyString : spectraStrings.keySet()) {
                spectrum = NMRShiftDBSpectrumToSpectrum(spectraStrings.get(spectrumPropertyString), nucleus);
                if(spectrum == null){
                    continue;
                }
                spectrumIndexInRecord = spectrumPropertyString.split("\\s")[spectrumPropertyString.split("\\s").length - 1];
                solvent = getSolvent(ac.getProperty("Solvent"), spectrumIndexInRecord);
                if(solvent == null){
                    continue;
                }
                spectrum.setSolvent(solvent);

                if(Utils.getAtomTypeIndicesByElement(ac, nucleus.replaceAll("\\d", "")).size() != spectrum.getSignalCount()){
                    continue;
                }

                spectra.add(new Object[]{spectrum, NMRShiftDBSpectrumToAssignment(spectraStrings.get(spectrumPropertyString), nucleus)});
            }
            spectraSet.add(spectra);
        }

        return spectraSet;
    }

    public static String getSolvent(final String solventPropertyString, final String spectrumIndexInRecord){
        final String[] solventPropertyStringSplit = solventPropertyString.split(":");
        String solvent;
        for (int i = 0; i < solventPropertyStringSplit.length; i++) {
            if (solventPropertyStringSplit[i].endsWith(spectrumIndexInRecord)) {
                solvent = solventPropertyStringSplit[i + 1];
                if(solvent.substring(solvent.length() - 1).matches("\\d")){
                    solvent = solvent.substring(0, solvent.length() - 1);
                }
                if(solvent.substring(solvent.length() - 1).matches("\\d")){
                    solvent = solvent.substring(0, solvent.length() - 1);
                }
                solvent = solvent.substring(0, solvent.length() - 1);

                return solvent;
            }
        }

        return null;
    }

    /**
     * Returns 3-tuples consisting of structure, spectrum and assignments
     * for each valid molecule record in the given NMRShiftDB file. Valid means
     * here that each molecule record has to contain the given spectrum
     * property string as well as the number of signals in that spectrum has to
     * be the same as atoms of that atom type in molecule.
     *
     * @param pathToNMRShiftDB path to NMRShiftDB file
     * @param NMRShiftDBSpectrumProperty spectrum property string to use
     * @return
     * @throws FileNotFoundException
     * @throws CDKException
     */
    public static HashMap<Integer, Object[]> getSSCComponentsFromNMRShiftDB(final String pathToNMRShiftDB, final String NMRShiftDBSpectrumProperty) throws FileNotFoundException, CDKException {
        final HashMap<Integer, Object[]> structureSetWithSpectra = new HashMap<>();
        final IteratingSDFReader iterator = new IteratingSDFReader(
                new FileReader(pathToNMRShiftDB),
                SilentChemObjectBuilder.getInstance()
        );
        IAtomContainer ac;
        Spectrum spectrum;
        Assignment assignment;
        final String nucleus = getNucleusFromNMRShiftDBSpectrumProperty(NMRShiftDBSpectrumProperty);
        final String spectrumIndexInRecord = NMRShiftDBSpectrumProperty.split("\\s")[NMRShiftDBSpectrumProperty.split("\\s").length - 1];
        while (iterator.hasNext()) {
            ac = iterator.next();
            // skip molecules which not contain any of requested spectrum information
            if(ac.getProperty(NMRShiftDBSpectrumProperty) == null){
                continue;
            }
            spectrum = NMRShiftDBSpectrumToSpectrum(ac.getProperty(NMRShiftDBSpectrumProperty), nucleus);
            // if no spectrum could be built or the number of signals in spectrum is different than the atom number in molecule
            if((spectrum == null) || Utils.getAtomTypeIndicesByElement(ac, nucleus.replaceAll("\\d", "")).size() != spectrum.getSignalCount()){
                continue;
            }
            if(ac.getProperty("Solvent") != null){
                spectrum.setSolvent(getSolvent(ac.getProperty("Solvent"), spectrumIndexInRecord));
            }
            if(ac.getProperty("Field Strength [MHz]") != null){
                for (final String fieldStrength : ac.getProperty("Field Strength [MHz]").toString().split("\\s")) {
                    if (fieldStrength.startsWith(spectrumIndexInRecord + ":")) {
                        try {
                            spectrum.setSpectrometerFrequency(Double.parseDouble(fieldStrength.split(spectrumIndexInRecord + ":")[1]));
                        } catch (NumberFormatException e) {
//                            spectrum.setSpectrometerFrequency(null);
                        }
                        break;
                    }
                }
            }

            assignment = NMRShiftDBSpectrumToAssignment(ac.getProperty(NMRShiftDBSpectrumProperty), nucleus);
//            if ((ac != null) && (spectrum != null)) {
                structureSetWithSpectra.put(structureSetWithSpectra.size(), new Object[]{ac, spectrum, assignment});
//            }

            Utils.setAromaticitiesInAtomContainer(ac);
        }

        return structureSetWithSpectra;
    }

    /**
     * Returns a hashmap containing combined keys (by "_") of solvents
     * and lists of calculated deviations between all given spectra for a
     * nucleus in molecule record as values. <br>
     * Here, only molecule records in NMRShiftDB file are considered which have
     * at least two different spectra for same nucleus. <br>
     * Example: "Spectrum 13C 0", "Spectrum 13C 1" will be used for given
     * nucleus 13C.
     *
     *
     * @param pathToNMRShiftDB
     * @param nucleus
     * @return
     * @throws FileNotFoundException
     * @throws CDKException
     */
    public static HashMap<String, ArrayList<Double>> getSolventDeviations(final String pathToNMRShiftDB, final String nucleus) throws FileNotFoundException, CDKException{
        int signalCount;
        Spectrum spectrum;
        Assignment assignment;
        final ArrayList<ArrayList<Object[]>> spectraSets = getSpectraFromNMRShiftDB(pathToNMRShiftDB, nucleus);
        HashMap<Integer, ArrayList<Double>> shiftsPerAtom;
        HashMap<Integer, ArrayList<String>> solventsPerAtom;
        ArrayList<String> solvents;
        String[] solventsToSort;

        final HashMap<String, ArrayList<Double>> deviations = new HashMap<>();
        String combiKey;

        for (final ArrayList<Object[]> spectraSetInRecord : spectraSets) {
            shiftsPerAtom = new HashMap<>();
            solventsPerAtom = new HashMap<>();
            signalCount = -1;
            for (final Object[] spectrumAndAssignment : spectraSetInRecord) {
                spectrum = (Spectrum) spectrumAndAssignment[0];
                assignment = (Assignment) spectrumAndAssignment[1];
                if (signalCount == -1) {
                    signalCount = spectrum.getSignalCount();
                } else if (signalCount != spectrum.getSignalCount()) {
                    continue;
                }
                for (final int atomIndex : assignment.getAtomIndices(0)) {
                    if (!shiftsPerAtom.containsKey(atomIndex)) {
                        shiftsPerAtom.put(atomIndex, new ArrayList<>());
                        solventsPerAtom.put(atomIndex, new ArrayList<>());
                    }
                    shiftsPerAtom.get(atomIndex).add(spectrum.getSignal(assignment.getSignalIndex(0, atomIndex)).getShift(0));
                    solventsPerAtom.get(atomIndex).add(spectrum.getSolvent());
                }
            }
            if (shiftsPerAtom.isEmpty() || (shiftsPerAtom.get(Collections.min(shiftsPerAtom.keySet())).size() < 2)) {
                continue;
            }
            solvents = new ArrayList<>(solventsPerAtom.get(Collections.min(solventsPerAtom.keySet())));
//                if(Collections.frequency(solvents, "Unreported") + Collections.frequency(solvents, "Unknown") > solvents.size() - 2){
//                    continue;
//                }

            for (final int atomIndex : shiftsPerAtom.keySet()) {
                for (int s1 = 0; s1 < solvents.size(); s1++) {
//                        if(solvents.get(s1).equals("Unreported") || solvents.get(s1).equals("Unknown")){
//                            continue;
//                        }
                    for (int s2 = s1 + 1; s2 < solvents.size(); s2++) {
//                            if (solvents.get(s2).equals("Unreported") || solvents.get(s2).equals("Unknown")) {
//                                continue;
//                            }
                        solventsToSort = new String[2];
                        solventsToSort[0] = solvents.get(s1);
                        solventsToSort[1] = solvents.get(s2);
                        Arrays.sort(solventsToSort);
                        combiKey = solventsToSort[0] + "_" + solventsToSort[1];
                        if (!deviations.containsKey(combiKey)) {
                            deviations.put(combiKey, new ArrayList<>());
                        }
                        deviations.get(combiKey).add(Math.abs(shiftsPerAtom.get(atomIndex).get(s1) - shiftsPerAtom.get(atomIndex).get(s2)));
                    }
                }
            }
        }

        return deviations;
    }

    /**
     *
     * @param pathToDB
     * @return
     * @throws FileNotFoundException
     * @deprecated
     */
    public static HashSet<String> getAtomTypesInDB(final String pathToDB) throws FileNotFoundException{
        final HashSet<String> atomTypes = new HashSet<>();
        final IteratingSDFReader iterator = new IteratingSDFReader(
                new FileReader(pathToDB),
                SilentChemObjectBuilder.getInstance()
        );
        while (iterator.hasNext()) {
            atomTypes.addAll(Utils.getAtomTypesInAtomContainer(iterator.next()));
        }

        return atomTypes;
    }

    // currently only for 1D spectra
    public static HashMap<String, String> getSpectraStrings(final IAtomContainer ac, final String nucleus) {
        final ArrayList<String> props = (ArrayList<String>) (ArrayList<?>) (new ArrayList<>(ac.getProperties().keySet()));
        final HashMap<String, String> spectra = new HashMap<>();
        for (final String prop : props) {
            if (prop.startsWith("Spectrum " + nucleus)) {
                spectra.put(prop, ac.getProperty(prop));
            }
        }

        return spectra;
    }

    /**
     * Creates a two dimensional array of a given NMRShiftDB NMR entry
     * with all signal shift values, intensities, multiplicities and atom indices.
     *
     * @param NMRShiftDBSpectrum
     * @return two dimensional array:
     * 1. dimension: signal index (row);
     * 2. dimension: signal shift value (column 1), signal intensity (column 2),
     * signal multiplicity (column 3), atom index in structure (column 4)
     */
    public static String[][] parseNMRShiftDBSpectrum(final String NMRShiftDBSpectrum){
        if(NMRShiftDBSpectrum.trim().isEmpty()){
            return new String[][]{};
        }
        String[] signalSplit;
        final String[] shiftsSplit = NMRShiftDBSpectrum.split("\\|");
        final String[][] values = new String[shiftsSplit.length][4];
        for (int i = 0; i < shiftsSplit.length; i++) {
            signalSplit = shiftsSplit[i].split(";");
            values[i][0] = signalSplit[0]; // shift value
            values[i][1] = signalSplit[1].substring(0, signalSplit[1].length() - 1); // intensity
            values[i][2] = signalSplit[1].substring(signalSplit[1].length() - 1); // multiplicity
            values[i][3] = signalSplit[2]; // atom index
        }

        return values;
    }

    /**
     * Sets shifts, intensities and implicit hydrogen counts in atoms of an atom container
     * by means of given spectrum property string.
     *
     * @param ac IAtomContainer to set
     * @param NMRShiftDBSpectrum Property string of spectrum in NMRShiftDB format.
     * @return
     *
     * @see MongoDB#parseNMRShiftDBSpectrum(String)
     * @see Utils#getHydrogenCountFromMultiplicity(String)
     * @deprecated
     */
    public static boolean setNMRShiftDBShiftsToAtomContainer(final IAtomContainer ac, final String NMRShiftDBSpectrum){
        if (ac.getProperty(NMRShiftDBSpectrum) == null) {
            return false;
        }
        final String[][] spectrumStringArray = parseNMRShiftDBSpectrum(ac.getProperty(NMRShiftDBSpectrum));

        Integer atomIndexSpectrum;
//        String multiplicity;
        Double shift;

        for (int i = 0; i < spectrumStringArray.length; i++) {
            atomIndexSpectrum = Integer.parseInt(spectrumStringArray[i][3]);
            shift = Double.parseDouble(spectrumStringArray[i][0]);
//            multiplicity = spectrumStringArray[i][3];
            if(Utils.checkIndexInAtomContainer(ac, atomIndexSpectrum)){
                ac.getAtom(atomIndexSpectrum).setProperty(Utils.getNMRShiftConstant(ac.getAtom(atomIndexSpectrum).getSymbol()), shift);
//                ac.getAtom(atomIndexSpectrum).setImplicitHydrogenCount(Utils.getHydrogenCountFromMultiplicity(multiplicity));
            }
        }

        return true;
    }

    public static String getNucleusFromNMRShiftDBSpectrumProperty(final String NMRShiftDBSpectrumProperty){
        return NMRShiftDBSpectrumProperty.split(" ")[1];
    }

    public static Spectrum NMRShiftDBSpectrumToSpectrum(final String NMRShiftDBSpectrum, final String nucleus){
        if ((NMRShiftDBSpectrum == null) || NMRShiftDBSpectrum.trim().isEmpty()) {
            return null;
        }
        final String[][] spectrumStringArray = parseNMRShiftDBSpectrum(NMRShiftDBSpectrum);
        final Spectrum spectrum = new Spectrum(new String[]{nucleus});
        String multiplicity;
        Double shift, intensity;
        try {
            for (int i = 0; i < spectrumStringArray.length; i++) {
                shift = Double.parseDouble(spectrumStringArray[i][0]);
                intensity = Double.parseDouble(spectrumStringArray[i][1]);
                multiplicity = spectrumStringArray[i][2];
                spectrum.addSignal(new Signal(new String[]{nucleus}, new Double[]{shift}, multiplicity, intensity));
            }
            spectrum.detectEquivalences();
        } catch (Exception e) {

            return null;
        }

        return spectrum;
    }

    public static Assignment NMRShiftDBSpectrumToAssignment(final String NMRShiftDBSpectrum, final String nucleus) {
        if ((NMRShiftDBSpectrum == null) || NMRShiftDBSpectrum.trim().isEmpty()) {
            return null;
        }
        final String[][] NMRShiftDBSpectrumStringArray = parseNMRShiftDBSpectrum(NMRShiftDBSpectrum);
        final Spectrum spectrum = NMRShiftDBSpectrumToSpectrum(NMRShiftDBSpectrum, nucleus);
        final Assignment assignment = new Assignment(spectrum);
        for (int i = 0; i < NMRShiftDBSpectrumStringArray.length; i++) {
            assignment.setAssignment(0, i, new Integer(NMRShiftDBSpectrumStringArray[i][3]));
        }

        return assignment;
    }
}
