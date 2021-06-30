package casekit.nmr.fragments.functionalgroup;

import casekit.nmr.dbservice.NMRShiftDB;
import casekit.nmr.fragments.fragmentation.Fragmentation;
import casekit.nmr.fragments.fragmentation.FragmentationUtilities;
import casekit.nmr.fragments.model.ConnectionTree;
import casekit.nmr.model.DataSet;
import casekit.nmr.utils.Utils;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

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
                        FragmentationUtilities.adjustNodeKeys(fragmentTree, structure);
                        FragmentationUtilities.closeRings(fragmentTree, structure);

                        fragmentTrees.add(fragmentTree);
                    }
                }
                FragmentationUtilities.removeDuplicates(fragmentTrees);
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
