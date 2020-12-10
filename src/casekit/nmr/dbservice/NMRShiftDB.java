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

package casekit.nmr.dbservice;

import casekit.nmr.Utils;
import casekit.nmr.model.Assignment;
import casekit.nmr.model.DataSet;
import casekit.nmr.model.Signal;
import casekit.nmr.model.Spectrum;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.*;

public class NMRShiftDB {

    public static String getSolvent(final String solventPropertyString, final String spectrumIndexInRecord) {
        final String[] solventPropertyStringSplit = solventPropertyString.split(":");
        String solvent;
        for (int i = 0; i < solventPropertyStringSplit.length; i++) {
            if (solventPropertyStringSplit[i].endsWith(spectrumIndexInRecord)) {
                solvent = solventPropertyStringSplit[i + 1];
                if (solvent.substring(solvent.length() - 1).matches("\\d")) {
                    solvent = solvent.substring(0, solvent.length() - 1);
                }
                if (solvent.substring(solvent.length() - 1).matches("\\d")) {
                    solvent = solvent.substring(0, solvent.length() - 1);
                }
                //                solvent = solvent.substring(0, solvent.length() - 1);
                solvent = solvent.trim();

                return solvent;
            }
        }

        return null;
    }

    public static List<String> getSpectraProperties1D(final IAtomContainer ac, final String nucleus) {
        final List<String> spectraProperties1D = new ArrayList<>();
        for (final Object obj : ac.getProperties().keySet()) {
            if (obj instanceof String && ((String) obj).startsWith("Spectrum " + nucleus)) {
                spectraProperties1D.add((String) obj);
            }
        }

        return spectraProperties1D;
    }

    /**
     * Returns a {@link DataSet} class object
     * for each valid molecule record in the given NMRShiftDB file. Valid means
     * here that each molecule record has to contain the given spectrum
     * property string as well as the number of signals in that spectrum has to
     * be the same as atoms of that atom type in molecule.
     *
     * @param pathToNMRShiftDB path to NMRShiftDB file
     * @param nuclei           nuclei to get the spectra for
     *
     * @return
     *
     * @throws FileNotFoundException
     * @throws CDKException
     * @see DataSet
     */
    public static Collection<DataSet> getDataSetsFromNMRShiftDB(final String pathToNMRShiftDB, final String[] nuclei) throws FileNotFoundException, CDKException {
        final Collection<DataSet> dataSets = new ArrayList<>();
        final IteratingSDFReader iterator = new IteratingSDFReader(new FileReader(pathToNMRShiftDB), SilentChemObjectBuilder.getInstance());
        IAtomContainer structure;
        Spectrum spectrum;
        Assignment assignment;
        HashMap<String, String> meta;
        final CDKHydrogenAdder hydrogenAdder = CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance());

        List<String> spectraProperties1D;
        String[] split;
        String spectrumIndexInRecord;
        IMolecularFormula mf;
        List<Integer> explicitHydrogenIndices;

        while (iterator.hasNext()) {
            structure = iterator.next();
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(structure);
            explicitHydrogenIndices = Utils.getExplicitHydrogenIndices(structure);
            Collections.sort(explicitHydrogenIndices);
            if (!explicitHydrogenIndices.isEmpty()) {
                // remove explicit hydrogens
                Utils.removeAtoms(structure, "H");
            }
            hydrogenAdder.addImplicitHydrogens(structure);
            Utils.setAromaticityAndKekulize(structure);

            meta = new HashMap<>();
            meta.put("title", structure.getTitle());
            meta.put("id", structure.getProperty("nmrshiftdb2 ID"));
            mf = Utils.getMolecularFormulaFromAtomContainer(structure);
            meta.put("mf", Utils.molecularFormularToString(mf));

            for (final String nucleus : nuclei) {
                spectraProperties1D = getSpectraProperties1D(structure, nucleus);
                for (final String spectrumProperty1D : spectraProperties1D) {

                    split = spectrumProperty1D.split("\\s");
                    spectrumIndexInRecord = split[split.length - 1];

                    // skip molecules which do not contain any of requested spectrum information
                    spectrum = NMRShiftDBSpectrumToSpectrum(structure.getProperty(spectrumProperty1D), nucleus);
                    // if no spectrum could be built or the number of signals in spectrum is different than the atom number in molecule
                    if ((spectrum == null) || Utils.getDifferenceSpectrumSizeAndMolecularFormulaCount(spectrum, mf, 0) != 0) {
                        continue;
                    }
                    if (structure.getProperty("Solvent") != null) {
                        spectrum.setSolvent(getSolvent(structure.getProperty("Solvent"), spectrumIndexInRecord));
                    }
                    if (structure.getProperty("Field Strength [MHz]") != null) {
                        for (final String fieldStrength : structure.getProperty("Field Strength [MHz]").toString().split("\\s")) {
                            if (fieldStrength.startsWith(spectrumIndexInRecord + ":")) {
                                try {
                                    spectrum.setSpectrometerFrequency(Double.parseDouble(fieldStrength.split(spectrumIndexInRecord + ":")[1]));
                                } catch (NumberFormatException e) {
                                    //                                    e.printStackTrace();
                                }
                                break;
                            }
                        }
                    }

                    assignment = NMRShiftDBSpectrumToAssignment(structure.getProperty(spectrumProperty1D), nucleus);
                    if (assignment != null && !explicitHydrogenIndices.isEmpty()) {
                        int hCount;
                        for (int i = 0; i < assignment.getAssignmentsCount(); i++) {
                            hCount = 0;
                            for (int j = 0; j < explicitHydrogenIndices.size(); j++) {
                                if (explicitHydrogenIndices.get(j) >= assignment.getAssignment(0, i)) {
                                    break;
                                }
                                hCount++;
                            }
                            assignment.setAssignment(0, i, assignment.getAssignment(0, i) - hCount);
                        }
                    }

                    dataSets.add(new DataSet(structure, spectrum, assignment, meta));
                }
            }


        }

        return dataSets;
    }

    //    /**
    //     * Returns a hashmap containing combined keys (by "_") of solvents
    //     * and lists of calculated deviations between all given spectra for a
    //     * nucleus in molecule record as values. <br>
    //     * Here, only molecule records in NMRShiftDB file are considered which have
    //     * at least two different spectra for same nucleus. <br>
    //     * Example: "Spectrum 13C 0", "Spectrum 13C 1" will be used for given
    //     * nucleus 13C.
    //     *
    //     * @param pathToNMRShiftDB
    //     * @param nucleus
    //     *
    //     * @return
    //     *
    //     * @throws FileNotFoundException
    //     * @throws CDKException
    //     */
    //    public static HashMap<String, ArrayList<Double>> getSolventDeviations(final String pathToNMRShiftDB, final String nucleus) throws FileNotFoundException, CDKException {
    //        int signalCount;
    //        Spectrum spectrum;
    //        Assignment assignment;
    //        final ArrayList<ArrayList<Object[]>> spectraSets = getSpectraFromNMRShiftDB(pathToNMRShiftDB, nucleus);
    //        HashMap<Integer, ArrayList<Double>> shiftsPerAtom;
    //        HashMap<Integer, ArrayList<String>> solventsPerAtom;
    //        ArrayList<String> solvents;
    //        String[] solventsToSort;
    //
    //        final HashMap<String, ArrayList<Double>> deviations = new HashMap<>();
    //        String combiKey;
    //
    //        for (final ArrayList<Object[]> spectraSetInRecord : spectraSets) {
    //            shiftsPerAtom = new HashMap<>();
    //            solventsPerAtom = new HashMap<>();
    //            signalCount = -1;
    //            for (final Object[] spectrumAndAssignment : spectraSetInRecord) {
    //                spectrum = (Spectrum) spectrumAndAssignment[0];
    //                assignment = (Assignment) spectrumAndAssignment[1];
    //                if (signalCount == -1) {
    //                    signalCount = spectrum.getSignalCount();
    //                } else if (signalCount != spectrum.getSignalCount()) {
    //                    continue;
    //                }
    //                for (final int atomIndex : assignment.getAssignments(0)) {
    //                    if (!shiftsPerAtom.containsKey(atomIndex)) {
    //                        shiftsPerAtom.put(atomIndex, new ArrayList<>());
    //                        solventsPerAtom.put(atomIndex, new ArrayList<>());
    //                    }
    //                    shiftsPerAtom.get(atomIndex).add(spectrum.getSignal(assignment.getIndex(0, atomIndex)).getShift(0));
    //                    solventsPerAtom.get(atomIndex).add(spectrum.getSolvent());
    //                }
    //            }
    //            if (shiftsPerAtom.isEmpty() || (shiftsPerAtom.get(Collections.min(shiftsPerAtom.keySet())).size() < 2)) {
    //                continue;
    //            }
    //            solvents = new ArrayList<>(solventsPerAtom.get(Collections.min(solventsPerAtom.keySet())));
    //            //                if(Collections.frequency(solvents, "Unreported") + Collections.frequency(solvents, "Unknown") > solvents.size() - 2){
    //            //                    continue;
    //            //                }
    //
    //            for (final int atomIndex : shiftsPerAtom.keySet()) {
    //                for (int s1 = 0; s1 < solvents.size(); s1++) {
    //                    //                        if(solvents.get(s1).equals("Unreported") || solvents.get(s1).equals("Unknown")){
    //                    //                            continue;
    //                    //                        }
    //                    for (int s2 = s1 + 1; s2 < solvents.size(); s2++) {
    //                        //                            if (solvents.get(s2).equals("Unreported") || solvents.get(s2).equals("Unknown")) {
    //                        //                                continue;
    //                        //                            }
    //                        solventsToSort = new String[2];
    //                        solventsToSort[0] = solvents.get(s1);
    //                        solventsToSort[1] = solvents.get(s2);
    //                        Arrays.sort(solventsToSort);
    //                        combiKey = solventsToSort[0] + "_" + solventsToSort[1];
    //                        if (!deviations.containsKey(combiKey)) {
    //                            deviations.put(combiKey, new ArrayList<>());
    //                        }
    //                        deviations.get(combiKey).add(Math.abs(shiftsPerAtom.get(atomIndex).get(s1) - shiftsPerAtom.get(atomIndex).get(s2)));
    //                    }
    //                }
    //            }
    //        }
    //
    //        return deviations;
    //    }
    //
    //    /**
    //     * @param pathToDB
    //     *
    //     * @return
    //     *
    //     * @throws FileNotFoundException
    //     * @deprecated
    //     */
    //    public static Set<String> getAtomTypesInDB(final String pathToDB) throws FileNotFoundException {
    //        final HashSet<String> atomTypes = new HashSet<>();
    //        final IteratingSDFReader iterator = new IteratingSDFReader(new FileReader(pathToDB), SilentChemObjectBuilder.getInstance());
    //        while (iterator.hasNext()) {
    //            atomTypes.addAll(Utils.getAtomTypesInAtomContainer(iterator.next()));
    //        }
    //
    //        return atomTypes;
    //    }

    /**
     * Creates a two dimensional array of a given NMRShiftDB casekit.nmr entry
     * with all signal shift values, intensities, multiplicities and atom indices.
     *
     * @param NMRShiftDBSpectrum
     *
     * @return two dimensional array:
     * 1. dimension: signal index (row);
     * 2. dimension: signal shift value (column 1), signal intensity (column 2),
     * signal multiplicity (column 3), atom index in structure (column 4)
     */
    public static String[][] parseNMRShiftDBSpectrum(final String NMRShiftDBSpectrum) {
        if (NMRShiftDBSpectrum.trim().isEmpty()) {
            return new String[][]{};
        }
        String[] signalSplit;
        final String[] shiftsSplit = NMRShiftDBSpectrum.split("\\|");
        final String[][] values = new String[shiftsSplit.length][4];
        for (int i = 0; i < shiftsSplit.length; i++) {
            signalSplit = shiftsSplit[i].split(";");
            values[i][0] = signalSplit[0]; // shift value
            values[i][1] = signalSplit[1].toLowerCase().split("[a-z]")[0]; // intensity
            values[i][2] = signalSplit[1].split("\\d+\\.\\d+").length > 0 ? signalSplit[1].split("\\d+\\.\\d+")[1].toLowerCase() : ""; // multiplicity
            values[i][3] = signalSplit[2]; // atom index
        }

        return values;
    }

    public static String NMRShiftDBSpectrumToBasicTextSpectrum(final String NMRShiftDBSpectrum, final String nucleus, final String description) {
        if ((NMRShiftDBSpectrum == null) || NMRShiftDBSpectrum.trim().isEmpty()) {
            return null;
        }
        final StringBuilder basicSpectrum = new StringBuilder();
        // append description
        if (!description.trim().startsWith("//")) {
            basicSpectrum.append("// ");
        }
        basicSpectrum.append(description).append("\n");
        final String[][] spectrumStringArray = NMRShiftDB.parseNMRShiftDBSpectrum(NMRShiftDBSpectrum);
        try {
            for (int i = 0; i < spectrumStringArray.length; i++) {
                // append nucleus
                basicSpectrum.append(nucleus).append(", ");
                // append chemical shift
                basicSpectrum.append(Double.parseDouble(spectrumStringArray[i][0])).append(", ");
                // append multiplicity
                basicSpectrum.append(spectrumStringArray[i][2]).append(", ");
                // append intensity
                basicSpectrum.append(Double.parseDouble(spectrumStringArray[i][1])).append("\n");
            }
        } catch (Exception e) {
            return null;
        }

        return basicSpectrum.toString();
    }

    public static Spectrum NMRShiftDBSpectrumToSpectrum(final String NMRShiftDBSpectrum, final String nucleus) {
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
                spectrum.addSignal(new Signal(new String[]{nucleus}, new Double[]{shift}, multiplicity, "signal", intensity, 0));
            }
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

    //    public static Map<String, Map<String, Map<String, ArrayList<Double>>>> buildHybridizationDistributions(final String pathToDB) {
    //        // for atom type -> hybridization -> multiplicity -> shift list
    //        final Map<String, Map<String, Map<String, ArrayList<Double>>>> hybridizationDistributions = new HashMap<>();
    //
    //        try (final IteratingSDFReader iterator = new IteratingSDFReader(new FileReader(pathToDB), SilentChemObjectBuilder.getInstance())) {
    //            IAtom atom;
    //            String nucleus;
    //            IAtomContainer structure;
    //            List<String> spectraProperties13C, spectraProperties1H, spectraProperties15N, spectraProperties1D;
    //            Spectrum spectrum;
    //            Assignment assignment;
    //            Signal signal;
    //            while (iterator.hasNext()) {
    //                structure = iterator.next();
    //                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(structure);
    //                Utils.setAromaticity(structure);
    //
    //                spectraProperties13C = getSpectraProperties1D(structure, "13C");
    //                spectraProperties1H = getSpectraProperties1D(structure, "1H");
    //                spectraProperties15N = getSpectraProperties1D(structure, "15N");
    //
    //                for (int i = 0; i < structure.getAtomCount(); i++) {
    //                    atom = structure.getAtom(i);
    //                    if (!hybridizationDistributions.containsKey(atom.getSymbol())) {
    //                        hybridizationDistributions.put(atom.getSymbol(), new HashMap<>());
    //                    }
    //                    if (!hybridizationDistributions.get(atom.getSymbol()).containsKey(atom.getHybridization().name())) {
    //                        hybridizationDistributions.get(atom.getSymbol()).put(atom.getHybridization().name(), new HashMap<>());
    //                    }
    //
    //                    switch (atom.getSymbol()) {
    //                        case "C":
    //                            spectraProperties1D = spectraProperties13C;
    //                            nucleus = "13C";
    //                            break;
    //                        case "H":
    //                            spectraProperties1D = spectraProperties1H;
    //                            nucleus = "1H";
    //                            break;
    //                        case "N":
    //                            spectraProperties1D = spectraProperties15N;
    //                            nucleus = "15N";
    //                            break;
    //                        default:
    //                            spectraProperties1D = new ArrayList<>();
    //                            nucleus = "";
    //                            break;
    //                    }
    //
    //                    for (final String spectrumProperty1D : spectraProperties1D) {
    //                        spectrum = NMRShiftDBSpectrumToSpectrum(structure.getProperty(spectrumProperty1D), nucleus);
    //                        assignment = NMRShiftDBSpectrumToAssignment(structure.getProperty(spectrumProperty1D), nucleus);
    //                        signal = spectrum.getSignal(assignment.getIndex(0, i));
    //
    //                        if (signal != null && signal.getMultiplicity() != null) {
    //                            if (!hybridizationDistributions.get(atom.getSymbol()).get(atom.getHybridization().name()).containsKey(signal.getMultiplicity())) {
    //                                hybridizationDistributions.get(atom.getSymbol()).get(atom.getHybridization().name()).put(signal.getMultiplicity(), new ArrayList<>());
    //                            }
    //                            hybridizationDistributions.get(atom.getSymbol()).get(atom.getHybridization().name()).get(signal.getMultiplicity()).add(signal.getShift(0));
    //                        }
    //                    }
    //                }
    //            }
    //        } catch (IOException | CDKException e) {
    //            e.printStackTrace();
    //        }
    //
    //        System.out.println(hybridizationDistributions);
    //
    //        return hybridizationDistributions;
    //    }
}
