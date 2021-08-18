package casekit.nmr.utils;

import casekit.nmr.model.DataSet;
import casekit.nmr.model.ExtendedAdjacencyList;
import casekit.nmr.model.Spectrum;
import casekit.nmr.model.nmrdisplayer.Correlation;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Utils {

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
                return "s";
            case 1:
                return "d";
            case 2:
                return "t";
            case 3:
                return "q";
            default:
                return null;
        }
    }

    public static String getMultiplicityFromProtonsCount(final Correlation correlation) {
        if (correlation.getAtomType()
                       .equals("C")
                && correlation.getProtonsCount()
                              .size()
                == 1) {
            return getMultiplicityFromProtonsCount(correlation.getProtonsCount()
                                                              .get(0));
        }
        return null;
    }

    public static String getAtomTypeFromSpectrum(final Spectrum spectrum, final int dim) {
        if (spectrum.containsDim(dim)) {
            return getAtomTypeFromNucleus(spectrum.getNuclei()[dim]);
        }

        return null;
    }

    public static String getAtomTypeFromNucleus(final String nucleus) {
        final String[] nucleusSplit = nucleus.split("\\d");
        return nucleusSplit[nucleusSplit.length
                - 1];
    }

    public static IMolecularFormula getMolecularFormulaFromString(final String mf) {
        return MolecularFormulaManipulator.getMolecularFormula(mf, SilentChemObjectBuilder.getInstance());
    }

    public static String getSmilesFromAtomContainer(final IAtomContainer ac) throws CDKException {
        // SmiFlavor.Unique instead of SmiFlavor.Absolute because current errors with InChI generator
        final SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Unique);

        return smilesGenerator.create(ac);
    }

    public static String getAlphabeticMF(final String mf) {
        final StringBuilder mfAlphabeticStringBuilder = new StringBuilder();
        final Map<String, Integer> mfAlphabeticMap = new TreeMap<>(getMolecularFormulaElementCounts(mf));
        for (final Map.Entry<String, Integer> entry : mfAlphabeticMap.entrySet()) {
            mfAlphabeticStringBuilder.append(entry.getKey());
            if (entry.getValue()
                    > 1) {
                mfAlphabeticStringBuilder.append(entry.getValue());
            }
        }

        return mfAlphabeticStringBuilder.toString();
    }

    public static Map<String, Integer> getMolecularFormulaElementCounts(final String mf) {
        final LinkedHashMap<String, Integer> counts = new LinkedHashMap<>();
        final List<String> elements = new ArrayList<>();
        Matcher matcher = Pattern.compile("([A-Z][a-z]{0,1})")
                                 .matcher(mf);
        while (matcher.find()) {
            elements.add(matcher.group(1));
        }
        int count;
        for (final String element : elements) {
            matcher = Pattern.compile("("
                                              + element
                                              + "\\d+)")
                             .matcher(mf);
            count = 1;
            if (matcher.find()) {
                count = Integer.parseInt(matcher.group(1)
                                                .split(element)[1]);
            }
            counts.put(element, count);
        }

        return counts;
    }

    public static int getAtomTypeCount(final IAtomContainer structure, final String atomType) {
        return getAtomTypeIndicesByElement(structure, atomType).size();
    }

    public static int getAtomTypeCount(final String mf, final String atomType) {
        return MolecularFormulaManipulator.getElementCount(getMolecularFormulaFromString(mf), atomType);
    }

    public static boolean compareWithMolecularFormulaLessOrEqual(final IAtomContainer structure, final String mf) {
        if (mf
                != null
                && !mf.trim()
                      .isEmpty()) {
            for (final String atomType : getAtomTypesInAtomContainer(structure)) {
                if (getAtomTypeCount(structure, atomType)
                        > getAtomTypeCount(mf, atomType)) {
                    return false;
                }
            }
            return AtomContainerManipulator.getImplicitHydrogenCount(structure)
                    <= getAtomTypeCount(mf, "H");
        }

        return true;
    }

    public static boolean compareWithMolecularFormulaEqual(final IAtomContainer structure, final String mf) {
        if (mf
                != null
                && !mf.trim()
                      .isEmpty()) {
            for (final String atomType : getAtomTypesInAtomContainer(structure)) {
                if (getAtomTypeCount(structure, atomType)
                        != getAtomTypeCount(mf, atomType)) {
                    return false;
                }
            }
            return AtomContainerManipulator.getImplicitHydrogenCount(structure)
                    == Utils.getAtomTypeCount(mf, "H");
        }

        return true;
    }

    /**
     * Returns a hashmap consisting of lists of atom indices in an atom container.
     * This is done for all atom types (e.g. C or Br) in given atom container.
     *
     * @param ac IAtomContainer to look in
     *
     * @return
     *
     * @see #getAtomTypeIndicesByElement(IAtomContainer, String)
     */
    public static Map<String, List<Integer>> getAtomTypeIndices(final IAtomContainer ac) {

        final Map<String, List<Integer>> atomTypeIndices = new HashMap<>();
        final Set<String> atomTypes = new HashSet<>();
        for (final IAtom heavyAtom : AtomContainerManipulator.getHeavyAtoms(ac)) {
            atomTypes.add(heavyAtom.getSymbol());
        }
        for (final String atomType : atomTypes) {
            atomTypeIndices.put(atomType, getAtomTypeIndicesByElement(ac, atomType));
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
        final String atomType = getAtomTypeFromSpectrum(spectrum, dim);
        int atomsInMolFormula = 0;
        if (molFormula
                != null) {
            atomsInMolFormula = MolecularFormulaManipulator.getElementCount(molFormula, atomType);
        }
        return atomsInMolFormula
                - spectrum.getSignalCountWithEquivalences();
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

    public static boolean isValidBondAddition(final IAtomContainer ac, final int atomIndex, final IBond bondToAdd) {
        float bondOrderSum = getBondOrderSum(ac, atomIndex, true);
        bondOrderSum += getBondOrderAsNumeric(bondToAdd);

        //        System.out.print(atomIndex + " --> " + HOSECodeUtilities.getBondOrderSum(ac, atomIndex, true) + " + " + HOSECodeUtilities.getBondOrderAsNumeric(bondToAdd));
        final IAtom atom = ac.getAtom(atomIndex);
        // @TODO include different valencies: N3, N5, S2, S4, S6 etc.
        // -1 for cases with heterocyclic aromatics, like the N in the small aromatic ring in coffein if we want to add the bond to the CH3 group
        if (atom.isAromatic()
                && (atom.getSymbol()
                        .equals("N")
                || atom.getSymbol()
                       .equals("S")
                || atom.getSymbol()
                       .equals("P"))) {
            //            System.out.print("[ -1 ]");
            bondOrderSum -= 1;
        }
        //        System.out.print(" = " + bondOrderSum + " <= " + atom.getValency() + " ? -> " + (bondOrderSum <= atom.getValency()) + "\n");

        // @TODO include charges
        return bondOrderSum
                <= atom.getValency();
    }

    public static Boolean isSaturated(final IAtomContainer ac, final int atomIndex) {
        if (!checkIndexInAtomContainer(ac, atomIndex)) {
            return null;
        }
        return ac.getAtom(atomIndex)
                 .getValency()
                != null
                && getBondOrderSum(ac, atomIndex, true).intValue()
                >= ac.getAtom(atomIndex)
                     .getValency();
    }

    public static List<Integer> getUnsaturatedAtomIndices(final IAtomContainer ac) {
        final List<Integer> unsaturatedAtomIndices = new ArrayList<>();
        for (int i = 0; i
                < ac.getAtomCount(); i++) {
            // set the indices of unsaturated atoms in substructure
            if (!isSaturated(ac, i)) {
                unsaturatedAtomIndices.add(i);
            }
        }
        return unsaturatedAtomIndices;
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

    public static void addExplicitHydrogens(final IAtomContainer ac) throws CDKException {
        addImplicitHydrogens(ac);
        convertImplicitToExplicitHydrogens(ac);
    }

    public static void convertImplicitToExplicitHydrogens(final IAtomContainer ac) {
        AtomContainerManipulator.convertImplicitToExplicitHydrogens(ac);
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
     * @see #containsExplicitHydrogens(IAtomContainer)
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
        setAromaticity(ac);
        Kekulization.kekulize(ac);
    }

    public static void setAromaticity(final IAtomContainer ac, final Aromaticity aromaticity) throws CDKException {
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(ac);
        aromaticity.apply(ac);
    }

    public static void setAromaticityAndKekulize(final IAtomContainer ac,
                                                 final Aromaticity aromaticity) throws CDKException {
        setAromaticity(ac, aromaticity);
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

    public static boolean checkIndexInAtomContainer(final IAtomContainer ac, final int atomIndex) {
        return ((atomIndex
                >= 0)
                && atomIndex
                < ac.getAtomCount());
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
        if (!checkIndexInAtomContainer(ac, atomIndex)) {
            return null;
        }
        float bondsOrderSum = 0;
        final IAtom atom = ac.getAtom(atomIndex);
        for (final IBond bond : ac.getConnectedBondsList(atom)) {
            bondsOrderSum += getBondOrderAsNumeric(bond);
        }
        if (includeImplicitHydrogenCount
                && (atom.getImplicitHydrogenCount()
                != null)) {
            bondsOrderSum += atom.getImplicitHydrogenCount();
        }

        return bondsOrderSum;
    }

    public static DataSet atomContainerToDataSet(final IAtomContainer structure) throws CDKException {
        final CDKHydrogenAdder hydrogenAdder = CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance());
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(structure);
        hydrogenAdder.addImplicitHydrogens(structure);
        setAromaticityAndKekulize(structure);
        final Map<String, String> meta = new HashMap<>();
        meta.put("title", structure.getTitle());
        meta.put("mf", molecularFormularToString(getMolecularFormulaFromAtomContainer(structure)));
        try {
            final String smiles = getSmilesFromAtomContainer(structure);
            meta.put("smiles", smiles);
        } catch (final CDKException e) {
            e.printStackTrace();
        }
        final DataSet dataSet = new DataSet();
        dataSet.setStructure(new ExtendedAdjacencyList(structure));
        dataSet.setMeta(meta);

        return dataSet;
    }
}
