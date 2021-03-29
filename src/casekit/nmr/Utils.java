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
package casekit.nmr;


import casekit.nmr.model.Spectrum;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * @author Michael Wenk [https://github.com/michaelwenk]
 * @deprecated
 */
public class Utils {

    /**
     * Returns a hashmap consisting of lists of atom indices in an atom container.
     * This is done for all atom types (e.g. C or Br) in given atom container.
     *
     * @param ac IAtomContainer to look in
     *
     * @return
     *
     * @see #getAtomTypeIndicesByElement(org.openscience.cdk.interfaces.IAtomContainer, java.lang.String)
     */
    public static Map<String, List<Integer>> getAtomTypeIndices(final IAtomContainer ac) {

        final Map<String, List<Integer>> atomTypeIndices = new HashMap<>();
        final Set<String> atomTypes = new HashSet<>();
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
     * @param ac       IAtomContainer to use for search
     * @param atomType Atom type to find in atom container
     *
     * @return
     */
    public static List<Integer> getAtomTypeIndicesByElement(final IAtomContainer ac, final String atomType) {

        final ArrayList<Integer> indices = new ArrayList<>();
        for (int i = 0; i
                < ac.getAtomCount(); i++) {
            if (ac.getAtom(i)
                  .getSymbol()
                  .equals(atomType)) {
                indices.add(i);
            }
        }

        return indices;
    }


    public static IMolecularFormula getMolecularFormulaFromAtomContainer(final IAtomContainer ac) {
        return MolecularFormulaManipulator.getMolecularFormula(ac);
    }


    public static String molecularFormularToString(final IMolecularFormula molecularFormula) {
        return MolecularFormulaManipulator.getString(molecularFormula);
    }

    public static int getDifferenceSpectrumSizeAndMolecularFormulaCount(final Spectrum spectrum,
                                                                        final IMolecularFormula molFormula,
                                                                        final int dim) throws CDKException {
        if (!spectrum.containsDim(dim)) {
            throw new CDKException(Thread.currentThread()
                                         .getStackTrace()[2].getClassName()
                                           + "."
                                           + Thread.currentThread()
                                                   .getStackTrace()[2].getMethodName()
                                           + ": invalid dimension in spectrum given");
        }
        final String atomType = casekit.nmr.utils.Utils.getAtomTypeFromSpectrum(spectrum, dim);
        int atomsInMolFormula = 0;
        if (molFormula
                != null) {
            atomsInMolFormula = MolecularFormulaManipulator.getElementCount(molFormula, atomType);
        }
        return atomsInMolFormula
                - spectrum.getSignalCountWithEquivalences();
    }

    /**
     * Specified for carbons only -> not generic!!!
     *
     * @param protonsCount
     *
     * @return
     */
    public static String getMultiplicityFromProtonsCount(final int protonsCount) {
        switch (protonsCount) {
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
     * Returns the casekit.nmr isotope identifier for a given element, e.g. C -> 13C.
     * Elements defined so far: C, H, N, P, F, D, O, S, Si, B, Pt.
     *
     * @param element element's symbol (e.g. "C")
     *
     * @return
     */
    public static String getIsotopeIdentifier(final String element) {
        switch (element) {
            case "C":
                return "13C";
            case "H":
                return "1H";
            case "N":
                return "15N";
            case "P":
                return "31P";
            case "F":
                return "19F";
            case "O":
                return "17O";
            case "S":
                return "33S";
            case "Si":
                return "29Si";
            case "B":
                return "11B";
            case "Pt":
                return "195Pt";
            default:
                return element;
        }
    }


    public static Set<String> getAtomTypesInAtomContainer(final IAtomContainer ac) {
        final HashSet<String> atomTypes = new HashSet<>();
        for (final IAtom atom : ac.atoms()) {
            atomTypes.add(atom.getSymbol());
        }

        return atomTypes;
    }


    /**
     * Detects outliers in given array list of input values and removes them. <br>
     * Here, outliers are those which are outside of a calculated lower and upper bound (whisker).
     * The interquartile range (IQR) of the input values is therefore multiplied with a given value
     * for whisker creation.
     *
     * @param input         list of values to process
     * @param multiplierIQR multiplier for IQR to use for lower and upper bound creation
     *
     * @return new array list without values outside the generated boundaries
     */
    public static List<Double> removeOutliers(final List<Double> input, final double multiplierIQR) {
        final ArrayList<Double> inputWithoutOutliers = new ArrayList<>(input);
        inputWithoutOutliers.removeAll(Utils.getOutliers(inputWithoutOutliers, multiplierIQR));

        return inputWithoutOutliers;
    }

    /**
     * @param input
     *
     * @return
     */
    public static List<Double> getOutliers(final List<Double> input, final double multiplierIQR) {
        final ArrayList<Double> outliers = new ArrayList<>();
        if (input.size()
                <= 1) {
            return outliers;
        }
        Collections.sort(input);
        final ArrayList<Double> data1 = new ArrayList<>(input.subList(0, input.size()
                / 2));
        final ArrayList<Double> data2;
        if (input.size()
                % 2
                == 0) {
            data2 = new ArrayList<>(input.subList(input.size()
                                                          / 2, input.size()));
        } else {
            data2 = new ArrayList<>(input.subList(input.size()
                                                          / 2
                                                          + 1, input.size()));
        }
        final double q1 = getMedian(new ArrayList<>(data1));
        final double q3 = getMedian(new ArrayList<>(data2));
        final double iqr = q3
                - q1;
        final double lowerBound = q1
                - multiplierIQR
                * iqr;
        final double upperBound = q3
                + multiplierIQR
                * iqr;
        for (int i = 0; i
                < input.size(); i++) {
            if ((input.get(i)
                    < lowerBound)
                    || (input.get(i)
                    > upperBound)) {
                outliers.add(input.get(i));
            }
        }
        //        System.out.println("input size: " + input.size());
        //        System.out.println("output size: " + outliers.size());
        return outliers;
    }


    /**
     * @param data
     *
     * @return
     */
    public static Double getMedian(final List<Double> data) {
        if ((data
                == null)
                || data.isEmpty()) {
            return null;
        }
        if (data.size()
                == 1) {
            return data.get(0);
        }
        Collections.sort(data);
        if (data.size()
                % 2
                == 1) {
            return data.get(data.size()
                                    / 2);
        } else {
            return (data.get(data.size()
                                     / 2
                                     - 1)
                    + data.get(data.size()
                                       / 2))
                    / 2.0;
        }
    }


    /**
     * @param data
     *
     * @return
     */
    public static Double getMean(final Collection<Double> data) {
        if ((data
                == null)
                || data.isEmpty()) {
            return null;
        }
        double sum = 0;
        int nullCounter = 0;
        for (final Double d : data) {
            if (d
                    != null) {
                sum += d;
            } else {
                nullCounter++;
            }
        }
        return ((data.size()
                - nullCounter)
                != 0)
               ? (sum
                / (data.size()
                - nullCounter))
               : null;
    }

    /**
     * @param data
     *
     * @return
     */
    public static Double getStandardDeviation(final List<Double> data) {
        if ((data
                == null)
                || data.isEmpty()) {
            return null;
        }
        final Double variance = Utils.getVariance(data);

        return (variance
                != null)
               ? Math.sqrt(variance)
               : null;
    }

    public static Double getVariance(final Collection<Double> data) {
        if ((data
                == null)
                || data.isEmpty()) {
            return null;
        }
        final int nullCounter = Collections.frequency(data, null);
        double quadrSum = 0.0;
        final Double mean = Utils.getMean(data);
        if (mean
                == null) {
            return null;
        }
        for (final Double d : data) {
            if (d
                    != null) {
                quadrSum += Math.pow(d
                                             - mean, 2);
            }
        }

        return ((data.size()
                - nullCounter)
                != 0)
               ? (quadrSum
                / (data.size()
                - nullCounter))
               : null;
    }


    /**
     * @param data
     *
     * @return
     */
    public static Double getMean(final Double[] data) {
        if ((data
                == null)
                || (data.length
                == 0)) {
            return null;
        }
        double sum = 0;
        int nullCounter = 0;
        for (final Double d : data) {
            if (d
                    != null) {
                sum += d;
            } else {
                nullCounter++;
            }
        }
        return ((data.length
                - nullCounter)
                != 0)
               ? (sum
                / (data.length
                - nullCounter))
               : null;
    }

    public static Map<String, Double> getMean(final Map<String, ArrayList<Double>> lookup) {

        final HashMap<String, Double> means = new HashMap<>();
        Double meanInList;
        for (final String key : lookup.keySet()) {
            meanInList = Utils.getMean(lookup.get(key));
            if (meanInList
                    != null) {
                means.put(key, meanInList);
            }
        }

        return means;
    }

    public static boolean isValidBondAddition(final IAtomContainer ac, final int atomIndex, final IBond bondToAdd) {
        float bondOrderSum = Utils.getBondOrderSum(ac, atomIndex, true);
        bondOrderSum += Utils.getBondOrderAsNumeric(bondToAdd);

        //        System.out.print(atomIndex + " --> " + Utils.getBondOrderSum(ac, atomIndex, true) + " + " + Utils.getBondOrderAsNumeric(bondToAdd));
        final IAtom atom = ac.getAtom(atomIndex);
        // -1 for cases with heterocyclic aromatics, like the N in the small aromatic ring in coffein if we want to add the bond to the CH3 group
        if (atom.isAromatic()
                && (!atom.getSymbol()
                         .equals("C"))) {
            //            System.out.print("[ -1 ]");
            bondOrderSum -= 1;
        }
        //        System.out.print(" = " + bondOrderSum + " <= " + atom.getValency() + " ? -> " + (bondOrderSum <= atom.getValency()) + "\n");

        // @TODO including charges
        return bondOrderSum
                <= atom.getValency();
    }


    /**
     * @param pathToFile
     *
     * @return
     */
    public static String getFileFormat(final String pathToFile) {

        if (pathToFile
                == null
                || pathToFile.trim()
                             .isEmpty()) {
            return "";
        }
        final String[] split = pathToFile.split("\\.");

        return split[split.length
                - 1];
    }


    //    /**
    //     * @param lookup
    //     *
    //     * @return
    //     */
    //    public static Map<String, Double> getRMSD(final Map<String, ArrayList<Double>> lookup) {
    //        final HashMap<String, Double> rmsd = new HashMap<>();
    //        Double rmsdInList;
    //        for (final String key : lookup.keySet()) {
    //            rmsdInList = casekit.nmr.utils.Utils.getRMSD(lookup.get(key));
    //            if (rmsdInList
    //                    != null) {
    //                rmsd.put(key, rmsdInList);
    //            }
    //        }
    //
    //        return rmsd;
    //    }

    public static Boolean isSaturated(final IAtomContainer ac, final int atomIndex) {
        if (!Utils.checkIndexInAtomContainer(ac, atomIndex)) {
            return null;
        }
        return Utils.getBondOrderSum(ac, atomIndex, true)
                    .intValue()
                >= ac.getAtom(atomIndex)
                     .getValency();
    }

    public static void addImplicitHydrogens(final IAtomContainer ac) throws CDKException {
        final CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(ac.getBuilder());
        IAtomType type;
        for (final IAtom atom : ac.atoms()) {
            type = matcher.findMatchingAtomType(ac, atom);
            AtomTypeManipulator.configure(atom, type);
        }
        final CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(ac.getBuilder());
        adder.addImplicitHydrogens(ac);
    }

    /**
     * @param lookup
     *
     * @return
     */
    public static Map<String, Double> getMedian(final Map<String, List<Double>> lookup) {

        final Map<String, Double> medians = new HashMap<>();
        Double medianInList;
        for (final String key : lookup.keySet()) {
            medianInList = Utils.getMedian(lookup.get(key));
            if (medianInList
                    != null) {
                medians.put(key, medianInList);
            }
        }

        return medians;
    }

    /**
     * @param hoseLookupToExtend
     * @param hoseLookup
     *
     * @deprecated
     */
    public static void combineHashMaps(final Map<String, List<Double>> hoseLookupToExtend,
                                       final Map<String, List<Double>> hoseLookup) {
        for (final String hose : hoseLookup.keySet()) {
            if (!hoseLookupToExtend.containsKey(hose)) {
                hoseLookupToExtend.put(hose, new ArrayList<>());
            }
            hoseLookupToExtend.get(hose)
                              .addAll(hoseLookup.get(hose));
        }
    }

    public static Double roundDouble(final Double value, final int decimalPlaces) {
        if (value
                == null) {
            return null;
        }
        final int decimalFactor = (int) (Math.pow(10, decimalPlaces));

        return (Math.round(value
                                   * decimalFactor)
                / (double) decimalFactor);
    }

    /**
     * Checks whether a structure contains explicit hydrogen atoms or not.
     *
     * @param ac structure to check
     *
     * @return
     */
    public static boolean containsExplicitHydrogens(final IAtomContainer ac) {
        return getExplicitHydrogenCount(ac)
                > 0;
    }

    /**
     * Stores all explicit hydrogens as implicit counter for the bonded heavy
     * atoms and removes those from the atom container. <br>
     * Also, a HashMap containing non-hydrogen atoms and its indices
     * before the removals will be returned which one can use for atom index
     * comparison (before and after the removals).
     *
     * @param ac the structure to lsd
     *
     * @return
     *
     * @see #containsExplicitHydrogens(org.openscience.cdk.interfaces.IAtomContainer)
     */
    public static Map<IAtom, Integer> convertExplicitToImplicitHydrogens(final IAtomContainer ac) {
        // create a list of atom indices which one can use for index comparison (before vs. after) after removing the explicit hydrogens
        final Map<IAtom, Integer> atomIndices = new HashMap<>();
        final List<IAtom> toRemoveList = new ArrayList<>();
        IAtom atomB;
        for (final IAtom atomA : ac.atoms()) {
            // check each atom whether it is an hydrogen;
            // if yes then store (increase) the number of implicit hydrogens 
            // for its bonded heavy atom
            if (atomA.getSymbol()
                     .equals("H")) {
                atomB = ac.getConnectedAtomsList(atomA)
                          .get(0);
                if (atomB.getImplicitHydrogenCount()
                        == null) {
                    atomB.setImplicitHydrogenCount(0);
                }
                atomB.setImplicitHydrogenCount(atomB.getImplicitHydrogenCount()
                                                       + 1);
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
     * @param ac
     *
     * @return
     */
    public static List<Integer> getExplicitHydrogenIndices(final IAtomContainer ac) {
        final List<Integer> explicitHydrogenIndicesList = new ArrayList<>();
        for (int i = 0; i
                < ac.getAtomCount(); i++) {
            if (ac.getAtom(i)
                  .getSymbol()
                  .equals("H")) {
                explicitHydrogenIndicesList.add(i);
            }
        }

        return explicitHydrogenIndicesList;
    }

    /**
     * @param ac
     *
     * @return
     */
    public static int getExplicitHydrogenCount(final IAtomContainer ac) {
        return getExplicitHydrogenIndices(ac).size();
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
     * @param ac       IAtomContainer object where to remove the atoms
     * @param atomType Atom type (element's name, e.g. C or Br)
     *
     * @return IAtomContainer where the atoms were removed
     */
    public static IAtomContainer removeAtoms(final IAtomContainer ac, final String atomType) {

        final ArrayList<IAtom> toRemoveList = new ArrayList<>();
        for (final IAtom atomA : ac.atoms()) {
            if (atomA.getSymbol()
                     .equals(atomType)) {// detect whether the current atom A is a from the given atom type
                toRemoveList.add(atomA);
            }
        }
        for (final IAtom iAtom : toRemoveList) {
            ac.removeAtom(iAtom);
        }

        return ac;
    }

    //    public static String getSpectrumNucleiAsString(final Spectrum spectrum) {
    //        String specID = "";
    //        for (int i = 0; i < spectrum.getNDim(); i++) {
    //            specID += spectrum.getNuclei()[i];
    //            if (i < spectrum.getNDim() - 1) {
    //                specID += "-";
    //            }
    //        }
    //
    //        return specID;
    //    }

    public static boolean checkIndexInAtomContainer(final IAtomContainer ac, final int atomIndex) {
        return ((atomIndex
                >= 0)
                && atomIndex
                < ac.getAtomCount());
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
        } catch (final InterruptedException e) {
            System.err.println("killing non-finished tasks!");
            executor.shutdownNow();
        }
    }

    /**
     * Returns the bond order for a numeric order value.
     *
     * @param orderAsNumeric
     *
     * @return
     */
    public static IBond.Order getBondOrder(final int orderAsNumeric) {
        for (final IBond.Order order : IBond.Order.values()) {
            if (order.numeric()
                    == orderAsNumeric) {
                return order;
            }
        }

        return null;
    }

    public static Float getBondOrderAsNumeric(final IBond bond) {
        if (bond
                == null) {
            return null;
        }
        final float bondOrderAsNumeric;
        if (bond.isAromatic()) {
            bondOrderAsNumeric = (float) 1.5;
        } else {
            bondOrderAsNumeric = bond.getOrder()
                                     .numeric();
        }

        return bondOrderAsNumeric;
    }

    public static Float getBondOrderSum(final IAtomContainer ac, final int atomIndex,
                                        final boolean includeImplicitHydrogenCount) {
        if (!Utils.checkIndexInAtomContainer(ac, atomIndex)) {
            return null;
        }
        float bondsOrderSum = 0;
        final IAtom atom = ac.getAtom(atomIndex);
        for (final IBond bond : ac.getConnectedBondsList(atom)) {
            bondsOrderSum += Utils.getBondOrderAsNumeric(bond);
        }
        if (includeImplicitHydrogenCount
                && (atom.getImplicitHydrogenCount()
                != null)) {
            bondsOrderSum += atom.getImplicitHydrogenCount();
        }

        return bondsOrderSum;
    }

}
