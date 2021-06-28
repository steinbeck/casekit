package casekit.nmr.fragmentation.functionalgroup;

import casekit.nmr.dbservice.NMRShiftDB;
import casekit.nmr.fragmentation.Fragmentation;
import casekit.nmr.fragmentation.FragmentationUtils;
import casekit.nmr.fragmentation.model.ConnectionTree;
import casekit.nmr.model.DataSet;
import casekit.nmr.model.Spectrum;
import casekit.nmr.utils.Match;
import casekit.nmr.utils.Utils;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class ErtlFunctionalGroupsUtilities {

    public static final List<DataSet> buildFunctionalGroupDataSets(final String pathToNMRShiftDB,
                                                                   final String[] nuclei) {
        final List<DataSet> functionalGroupDataSets = new ArrayList<>();
        try {
            final ErtlFunctionalGroupsFinder ertlFunctionalGroupsFinder = new ErtlFunctionalGroupsFinder(
                    ErtlFunctionalGroupsFinder.Mode.NO_GENERALIZATION);
            final List<DataSet> dataSetsFromNMRShiftDB = NMRShiftDB.getDataSetsFromNMRShiftDB(pathToNMRShiftDB, nuclei);
            List<DataSet> dataSetList;
            List<IAtomContainer> groups;
            List<ConnectionTree> fragmentTrees;
            ConnectionTree fragmentTree;
            IAtomContainer structure;
            String atomTypeInSpectrum;
            Aromaticity[] aromaticities;
            for (final DataSet dataSet : dataSetsFromNMRShiftDB) {
                structure = dataSet.getStructure()
                                   .toAtomContainer();
                aromaticities = buildDefaultAromaticities(structure);
                fragmentTrees = new ArrayList<>();
                for (final Aromaticity aromaticity : aromaticities) {
                    try {
                        Utils.setAromaticityAndKekulize(structure, aromaticity);
                        groups = ertlFunctionalGroupsFinder.find(structure, false);
                    } catch (final IllegalArgumentException | CDKException e) {
                        e.printStackTrace();
                        continue;
                    }
                    restoreOriginalEnvironmentalCarbons(groups, structure);
                    fragmentTrees = new ArrayList<>();
                    for (final IAtomContainer group : groups) {
                        // each group has to contain at least one atom of specific spectrum
                        atomTypeInSpectrum = casekit.nmr.utils.Utils.getAtomTypeFromNucleus(dataSet.getSpectrum()
                                                                                                   .getNuclei()[0]);
                        if (atomTypeInSpectrum.equals("H")) {
                            if (AtomContainerManipulator.getImplicitHydrogenCount(group)
                                    == 0) {
                                continue;
                            }
                        } else if (casekit.nmr.utils.Utils.getAtomTypeIndicesByElement(group, atomTypeInSpectrum)
                                                          .isEmpty()) {
                            continue;
                        }
                        fragmentTree = Fragmentation.buildFragmentTree(group, 0, null, new HashSet<>(), false);
                        FragmentationUtils.adjustNodeKeys(fragmentTree, structure);
                        FragmentationUtils.closeRings(fragmentTree, structure);

                        fragmentTrees.add(fragmentTree);
                    }
                }
                FragmentationUtils.removeDuplicates(fragmentTrees);
                dataSetList = Fragmentation.fragmentTreesToSubDataSets(dataSet, fragmentTrees);
                if (dataSetList
                        != null) {
                    functionalGroupDataSets.addAll(dataSetList);
                }
            }
        } catch (final IOException | CDKException e) {
            e.printStackTrace();
        }

        return functionalGroupDataSets;
    }

    public static LinkedHashMap<String, List<DataSet>> sortByFrequencies(
            final Map<String, List<DataSet>> functionalGroupDataSetsMap) {
        final LinkedHashMap<String, List<DataSet>> sortedCollection = new LinkedHashMap<>();
        final List<Map.Entry<String, List<DataSet>>> sortedFrequencies = getSortedFrequencies(
                functionalGroupDataSetsMap);
        for (final Map.Entry<String, List<DataSet>> frequency : sortedFrequencies) {
            sortedCollection.put(frequency.getKey(), frequency.getValue());
        }

        return sortedCollection;
    }

    public static List<Map.Entry<String, List<DataSet>>> getSortedFrequencies(
            final Map<String, List<DataSet>> functionalGroupDataSets) {
        return functionalGroupDataSets.entrySet()
                                      .stream()
                                      .sorted(Map.Entry.comparingByValue((list1, list2) -> -1
                                              * Integer.compare(list1.size(), list2.size())))
                                      .collect(Collectors.toList());
    }

    public static Map<String, Integer> countFrequencies(final List<DataSet> functionalGroupDataSets) {
        final Map<String, Integer> frequencies = new HashMap<>();
        String smiles;
        for (final DataSet functionalGroupDataSet : functionalGroupDataSets) {
            smiles = functionalGroupDataSet.getMeta()
                                           .get("smiles");
            if (smiles
                    != null) {
                frequencies.putIfAbsent(smiles, 0);
                frequencies.put(smiles, frequencies.get(smiles)
                        + 1);
            }
        }

        return frequencies;
    }

    public static Map<String, List<DataSet>> collectBySmiles(final List<DataSet> functionalGroupDataSets) {
        final Map<String, List<DataSet>> collection = new HashMap<>();
        String smiles;
        for (final DataSet functionalGroupDataSet : functionalGroupDataSets) {
            smiles = functionalGroupDataSet.getMeta()
                                           .get("smiles");
            if (smiles
                    != null) {
                collection.putIfAbsent(smiles, new ArrayList<>());
                collection.get(smiles)
                          .add(functionalGroupDataSet);
            }
        }

        return collection;
    }

    public static List<DataSet> findMatches(final List<DataSet> functionalGroupDataSets, final Spectrum querySpectrum,
                                            final String mf, final double shiftTol, final double maxAverageDeviation,
                                            final boolean checkMultiplicity) {
        final List<DataSet> matches = new ArrayList<>();
        for (final DataSet dataSet : functionalGroupDataSets) {
            if (isValidMatch(dataSet, querySpectrum, mf, shiftTol, maxAverageDeviation, checkMultiplicity)) {
                matches.add(dataSet);
            }
        }

        return matches;
    }

    public static boolean isValidMatch(final DataSet dataSet, final Spectrum querySpectrum, final String mf,
                                       final double shiftTol, final double maxAverageDeviation,
                                       final boolean checkMultiplicity) {
        final IMolecularFormula iMolecularFormula = MolecularFormulaManipulator.getMolecularFormula(mf,
                                                                                                    SilentChemObjectBuilder.getInstance());
        final IAtomContainer group = dataSet.getStructure()
                                            .toAtomContainer();

        if (!dataSet.getSpectrum()
                    .getNuclei()[0].equals(querySpectrum.getNuclei()[0])) {
            return false;
        }
        final String atomTypeInSpectrum = casekit.nmr.utils.Utils.getAtomTypeFromNucleus(dataSet.getSpectrum()
                                                                                                .getNuclei()[0]);
        if (atomTypeInSpectrum.equals("H")) {
            if (AtomContainerManipulator.getImplicitHydrogenCount(dataSet.getStructure()
                                                                         .toAtomContainer())
                    > MolecularFormulaManipulator.getElementCount(iMolecularFormula, atomTypeInSpectrum)) {
                return false;
            }
        } else {
            // check molecular formula with atom types in group
            if (!casekit.nmr.utils.Utils.compareWithMolecularFormulaLessOrEqual(group, mf)) {
                return false;
            }
            // do not allow unsaturated fragments with different size than given molecular formula
            if (casekit.nmr.utils.Utils.getUnsaturatedAtomIndices(group)
                                       .isEmpty()
                    && !casekit.nmr.utils.Utils.compareWithMolecularFormulaEqual(group, mf)) {
                return false;
            }
        }
        // check average deviation
        final Double averageDeviation = Match.calculateAverageDeviation(dataSet.getSpectrum(), querySpectrum, 0, 0,
                                                                        shiftTol, checkMultiplicity, true, true);

        if (averageDeviation
                == null
                || averageDeviation
                > maxAverageDeviation) {
            return false;
        }
        final Double rmsd = Match.calculateRMSD(dataSet.getSpectrum(), querySpectrum, 0, 0, shiftTol, checkMultiplicity,
                                                true, true);
        dataSet.getMeta()
               .put("avgDev", Double.toString(averageDeviation));
        dataSet.getMeta()
               .put("rmsd", Double.toString(rmsd));

        return true;
    }

    /**
     * Replaces the inserted environmental carbon atoms (new IAtom objects) by the
     * carbon IAtom objects from original structure.
     *
     * @param groups    groups created by ErtlFunctionalGroupFinder
     * @param structure original structure used to create the groups
     */
    private static void restoreOriginalEnvironmentalCarbons(final List<IAtomContainer> groups,
                                                            final IAtomContainer structure) {
        IAtomContainer group;
        IAtom connectedAtomInGroup;
        List<IAtom> atomList;
        for (int i = 0; i
                < groups.size(); i++) {
            group = groups.get(i);
            // convert explicit hydrogens back to implicit
            casekit.nmr.utils.Utils.convertExplicitToImplicitHydrogens(group);
            atomList = new ArrayList<>();
            // create a list (copy) of all atoms of the group because of atom removals and additions in group atom container
            group.atoms()
                 .spliterator()
                 .forEachRemaining(atomList::add);
            for (final IAtom atom : atomList) {
                // detect whether the current atom is an "unknown" one, inserted as new environmental IAtom object
                if (!structure.contains(atom)) {
                    // take its single parent from which should be in original
                    connectedAtomInGroup = group.getConnectedAtomsList(atom)
                                                .get(0);
                    // remove the inserted atom and the bond to it
                    group.removeBond(atom, connectedAtomInGroup);
                    group.removeAtom(atom);
                    // from the parent node search for neighboring carbons which are not already in the group
                    // and add them
                    for (final IAtom connectedAtomInOriginalStructure : structure.getConnectedAtomsList(
                            connectedAtomInGroup)) {
                        if (connectedAtomInOriginalStructure.getSymbol()
                                                            .equals("C")
                                && !group.contains(connectedAtomInOriginalStructure)) {
                            group.addAtom(connectedAtomInOriginalStructure);
                            group.addBond(structure.getBond(connectedAtomInGroup, connectedAtomInOriginalStructure));
                        }
                    }
                }
            }
        }
    }

    public static Aromaticity[] buildDefaultAromaticities(final IAtomContainer structure) {
        final CycleFinder cycles = Cycles.all(structure.getAtomCount());
        final ElectronDonation[] models = new ElectronDonation[]{ElectronDonation.cdk(),
                                                                 ElectronDonation.cdkAllowingExocyclic(),
                                                                 ElectronDonation.daylight(),
                                                                 ElectronDonation.piBonds()};
        final Aromaticity[] aromaticities = new Aromaticity[models.length];
        for (int i = 0; i
                < models.length; i++) {
            aromaticities[i] = new Aromaticity(models[i], cycles);
        }

        return aromaticities;
    }
}
