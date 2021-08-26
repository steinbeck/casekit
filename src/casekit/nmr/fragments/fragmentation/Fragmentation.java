package casekit.nmr.fragments.fragmentation;

import casekit.nmr.fragments.model.ConnectionTree;
import casekit.nmr.model.*;
import casekit.nmr.utils.Utils;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.RingSearch;

import java.util.*;
import java.util.stream.Collectors;

public class Fragmentation {

    /**
     * Builds connection trees starting at each atom in structure, including ring atoms.
     * Duplicates will be removed and subspectra and assignments set.
     *
     * @param dataSet         dataset with structure to build the fragments from
     * @param maxSphere       maximum spherical limit for single atom fragment creation
     * @param maxSphereRing   maximum spherical limit for ring atom fragment creation
     * @param withPseudoAtoms whether to place pseudo atoms in "outer" sphere
     *
     * @return
     *
     * @see #buildFragmentTrees(IAtomContainer, Integer, Integer, boolean)
     */
    public static List<DataSet> buildFragmentDataSets(final DataSet dataSet, final Integer maxSphere,
                                                      final Integer maxSphereRing, final boolean withPseudoAtoms) {

        final List<ConnectionTree> fragmentTrees = buildFragmentTrees(dataSet.getStructure()
                                                                             .toAtomContainer(), maxSphere,
                                                                      maxSphereRing, withPseudoAtoms);
        return fragmentTreesToSubDataSets(dataSet, fragmentTrees);
    }

    public static List<DataSet> fragmentTreesToSubDataSets(final DataSet dataSet,
                                                           final List<ConnectionTree> fragmentTrees) {
        final List<DataSet> fragmentDataSetList = new ArrayList<>();
        final IAtomContainer structure = dataSet.getStructure()
                                                .toAtomContainer();
        final Spectrum spectrum = dataSet.getSpectrum()
                                         .toSpectrum();
        final String spectrumAtomType = Utils.getAtomTypeFromSpectrum(spectrum, 0);
        List<Integer> substructureAtomIndices, signalIndices;
        IAtomContainer substructure;
        Spectrum subspectrum;
        Assignment subassignment;
        IAtom atomInStructure;
        Signal signal;
        DataSet subDataSet;
        Map<String, String> meta;
        String smiles;
        for (final ConnectionTree fragmentTree : fragmentTrees) {
            substructureAtomIndices = fragmentTree.getKeys();
            subspectrum = new Spectrum();
            subspectrum.setNuclei(dataSet.getSpectrum()
                                         .getNuclei());
            subspectrum.setSignals(new ArrayList<>());
            subassignment = new Assignment();
            subassignment.setNuclei(subspectrum.getNuclei());
            subassignment.initAssignments(0);
            for (int j = 0; j
                    < substructureAtomIndices.size(); j++) {
                if (substructureAtomIndices.get(j)
                        >= structure.getAtomCount()) {
                    // current node/atom is pseudo
                    continue;
                }
                atomInStructure = structure.getAtom(substructureAtomIndices.get(j));
                if (atomInStructure.getSymbol()
                                   .equals(spectrumAtomType)
                        || (spectrumAtomType.equals("H")
                        && atomInStructure.getImplicitHydrogenCount()
                        > 0)) {
                    signalIndices = dataSet.getAssignment()
                                           .getIndices(0, substructureAtomIndices.get(j));
                    if (signalIndices
                            == null
                            || signalIndices.isEmpty()) {
                        return null;
                    }

                    for (final int index : signalIndices) {
                        signal = spectrum.getSignal(index)
                                         .buildClone();
                        final int atomIndex = j;
                        final List<Integer> closestSignalIndexList = subspectrum.checkForEquivalences(signal,
                                                                                                      new double[]{0.0},
                                                                                                      true);
                        if (closestSignalIndexList.isEmpty()) {
                            signal.setEquivalencesCount(1);
                            subspectrum.addSignal(signal);
                            subassignment.addAssignment(0, new int[]{atomIndex});
                        } else {
                            final int signalIndex = closestSignalIndexList.get(0);
                            if (Arrays.stream(subassignment.getAssignment(0, signalIndex))
                                      .noneMatch(equiv -> equiv
                                              == atomIndex)) {
                                subspectrum.getSignal(signalIndex)
                                           .setEquivalencesCount(subspectrum.getSignal(signalIndex)
                                                                            .getEquivalencesCount()
                                                                         + 1); // + 1 because we add one atom only
                                subassignment.addAssignmentEquivalence(0, signalIndex, atomIndex);
                            }
                        }
                    }
                }
            }
            subspectrum.addMetaInfo("solvent", dataSet.getSpectrum()
                                                      .getMeta()
                                                      .get("solvent"));
            subspectrum.addMetaInfo("spectrometerFrequency", dataSet.getSpectrum()
                                                                    .getMeta()
                                                                    .get("spectrometerFrequency"));

            substructure = FragmentationUtilities.toAtomContainer(fragmentTree);
            subDataSet = new DataSet();
            subDataSet.setStructure(new StructureCompact(substructure));
            subDataSet.setSpectrum(new SpectrumCompact(subspectrum));
            subDataSet.setAssignment(subassignment);

            meta = new HashMap<>();
            try {
                smiles = casekit.nmr.utils.Utils.getSmilesFromAtomContainer(substructure);
                meta.put("smiles", smiles);

            } catch (final CDKException e) {
                e.printStackTrace();
            }
            //            meta.put("title", dataSet.getMeta()
            //                                     .get("title"));
            meta.put("id", dataSet.getMeta()
                                  .get("id"));
            meta.put("mf", Utils.molecularFormularToString(Utils.getMolecularFormulaFromAtomContainer(substructure)));
            subDataSet.setMeta(meta);

            fragmentDataSetList.add(subDataSet);
        }

        return fragmentDataSetList;
    }

    /**
     * Builds fragments as atom containers starting at each atom in structure, including ring atoms.
     * Duplicates will be removed.
     *
     * @param structure       structure to build the fragments from
     * @param maxSphere       maximum spherical limit for single atom fragment creation
     * @param maxSphereRing   maximum spherical limit for ring atom fragment creation
     * @param withPseudoAtoms whether to place pseudo atoms in "outer" sphere
     *
     * @return
     *
     * @see #buildFragmentTrees(IAtomContainer, Integer, Integer, boolean)
     * @see FragmentationUtilities#toAtomContainer(ConnectionTree)
     */
    public static List<IAtomContainer> buildFragments(final IAtomContainer structure, final Integer maxSphere,
                                                      final Integer maxSphereRing, final boolean withPseudoAtoms) {
        final List<ConnectionTree> fragmentTrees = buildFragmentTrees(structure, maxSphere, maxSphereRing,
                                                                      withPseudoAtoms);
        return fragmentTrees.stream()
                            .map(FragmentationUtilities::toAtomContainer)
                            .collect(Collectors.toList());
    }

    public static List<ConnectionTree> buildRingFragmentTrees(final IAtomContainer structure,
                                                              final Integer maxSphereRing,
                                                              final boolean withPseudoAtoms) {
        final List<ConnectionTree> ringFragmentTrees = new ArrayList<>();
        try {
            // build ring fragment trees from detected rings and extend by given maximum sphere for rings
            ConnectionTree connectionTreeRing, connectionTreeOuterSphere, subtreeToAdd;
            final IRingSet ringSet = Cycles.all(structure)//essential(structure)
                                           .toRingSet();
            final List<IAtomContainer> ringFragments = new ArrayList<>();
            for (int i = 0; i
                    < ringSet.getAtomContainerCount(); i++) {
                ringFragments.add(ringSet.getAtomContainer(i));
            }
            // add missing fused rings
            final RingSearch ringSearch = new RingSearch(structure);
            ringFragments.addAll(ringSearch.fusedRingFragments());
            List<Integer> atomIndicesInRing;
            Set<Integer> atomIndicesOutOfRing;
            IAtomContainer ringAtomContainer;
            for (int i = 0; i
                    < ringFragments.size(); i++) {
                ringAtomContainer = ringFragments.get(i);
                // add already "visited" ring nodes
                atomIndicesInRing = new ArrayList<>();
                for (int k = 0; k
                        < ringAtomContainer.getAtomCount(); k++) {
                    atomIndicesInRing.add(structure.indexOf(ringAtomContainer.getAtom(k)));
                }
                atomIndicesOutOfRing = new HashSet<>();
                for (int j = 0; j
                        < structure.getAtomCount(); j++) {
                    if (!atomIndicesInRing.contains(j)) {
                        atomIndicesOutOfRing.add(j);
                    }
                }
                connectionTreeRing = buildFragmentTree(structure, atomIndicesInRing.get(0), null, atomIndicesOutOfRing,
                                                       false);
                // add missing outer sphere nodes to ring
                for (int k = 0; k
                        < ringAtomContainer.getAtomCount(); k++) {
                    connectionTreeOuterSphere = Fragmentation.buildFragmentTree(structure, structure.indexOf(
                            ringAtomContainer.getAtom(k)), maxSphereRing, new HashSet<>(atomIndicesInRing), false);
                    if (connectionTreeOuterSphere.getMaxSphere(false)
                            == 0) {
                        continue;
                    }
                    for (final int key : connectionTreeOuterSphere.getNodeKeysInSphere(1)) {
                        subtreeToAdd = ConnectionTree.buildSubtree(connectionTreeOuterSphere, key);
                        if (!FragmentationUtilities.addToConnectionTree(connectionTreeRing,
                                                                        connectionTreeOuterSphere.getRootNode()
                                                                                                 .getKey(),
                                                                        subtreeToAdd, connectionTreeOuterSphere.getBond(
                                        connectionTreeOuterSphere.getRootNode()
                                                                 .getKey(), key))) {
                            continue;
                        }
                        atomIndicesInRing.addAll(subtreeToAdd.getKeys());
                    }
                }
                // close rings
                FragmentationUtilities.closeRings(connectionTreeRing, structure);
                // attach pseudo atoms if desired
                if (withPseudoAtoms) {
                    FragmentationUtilities.attachPseudoAtoms(connectionTreeRing, structure);
                }
                ringFragmentTrees.add(connectionTreeRing);
            }
        } catch (final CDKException e) {
            e.printStackTrace();
        }

        return ringFragmentTrees;
    }

    /**
     * Builds connection trees starting at each atom in structure, including ring atoms.
     * Duplicates are removed.
     *
     * @param structure       structure to build the fragments from
     * @param maxSphere       maximum spherical limit for single atom fragment creation
     * @param maxSphereRing   maximum spherical limit for ring atom fragment creation
     * @param withPseudoAtoms whether to place pseudo atoms in "outer" sphere
     *
     * @return
     *
     * @see #buildRingFragmentTrees(IAtomContainer, Integer, boolean)
     * @see #buildFragmentTree(IAtomContainer, int, Integer, Set, boolean)
     * @see FragmentationUtilities#removeDuplicates(List)
     */
    public static List<ConnectionTree> buildFragmentTrees(final IAtomContainer structure, final Integer maxSphere,
                                                          final Integer maxSphereRing, final boolean withPseudoAtoms) {
        // build fragment trees for rings
        final List<ConnectionTree> fragmentTrees = buildRingFragmentTrees(structure, maxSphereRing, withPseudoAtoms);
        // build fragment for each single atom
        for (int i = 0; i
                < structure.getAtomCount(); i++) {
            fragmentTrees.add(
                    Fragmentation.buildFragmentTree(structure, i, maxSphere, new HashSet<>(), withPseudoAtoms));
        }
        FragmentationUtilities.removeDuplicates(fragmentTrees);

        return fragmentTrees;
    }

    /**
     * Creates an atom container from a given connection tree built by using {@link #buildFragmentTree(IAtomContainer, int, Integer, Set, boolean)}.
     *
     * @param ac              atom container to go through
     * @param rootAtomIndex   root atom index to start from
     * @param maxSphere       spherical limit, a null value means no limit
     * @param exclude         atom indices which to exclude from search
     * @param withPseudoAtoms places pseudo atoms in the "outer" sphere
     *
     * @return connection tree
     *
     * @see #BFS(IAtomContainer, ConnectionTree, Queue, Set, Set, Integer)
     * @see #buildFragmentTree(IAtomContainer, int, Integer, Set, boolean)
     */
    public static IAtomContainer buildFragment(final IAtomContainer ac, final int rootAtomIndex,
                                               final Integer maxSphere, final Set<Integer> exclude,
                                               final boolean withPseudoAtoms) {
        return FragmentationUtilities.toAtomContainer(
                buildFragmentTree(ac, rootAtomIndex, maxSphere, exclude, withPseudoAtoms));
    }

    /**
     * Builds a fragment as connection tree from a given structure. <br>
     * Until a certain maximum sphere, each reachable next neighbor atom
     * is stored in a parent-child-relationship. <br>
     * And in addition, bonds between hetero atoms or carbon-hetero bonds will be kept. In such cases
     * the maximum spherical limit will be ignored.
     *
     * @param structure       atom container to go through
     * @param rootAtomIndex   root atom index to start from
     * @param maxSphere       spherical limit, a null value means no limit
     * @param exclude         atom indices which to exclude from search
     * @param withPseudoAtoms places pseudo atoms in the "outer" sphere
     *
     * @return connection tree
     */
    public static ConnectionTree buildFragmentTree(final IAtomContainer structure, final int rootAtomIndex,
                                                   final Integer maxSphere, final Set<Integer> exclude,
                                                   final boolean withPseudoAtoms) {
        // create queue and connection tree for BFS
        final Queue<int[]> queue = new LinkedList<>();
        queue.add(new int[]{rootAtomIndex, 0});
        final ConnectionTree connectionTree = new ConnectionTree(structure.getAtom(rootAtomIndex), rootAtomIndex);

        BFS(structure, connectionTree, queue, new HashSet<>(), exclude, maxSphere);

        // close rings
        FragmentationUtilities.closeRings(connectionTree, structure);
        // add pseudo atoms
        if (withPseudoAtoms) {
            FragmentationUtilities.attachPseudoAtoms(connectionTree, structure);
        }

        return connectionTree;
    }

    /**
     * Function for extending a given connection tree only containing
     * its root node (0th sphere) by means of Breadth-First-Search (BFS).
     * Until a certain maximum sphere, each reachable next neighbor atom
     * is stored in a parent-child-relationship.
     * And in addition, bonds between hetero atoms, carbon-hetero bonds and triple bonds will be kept.
     * In such cases the maximum spherical limit will be ignored.
     * <br> <br>
     * This method follows the rules given in section 7.4.1 of <br>
     * "Contemporary Computer-Assisted Approaches to Molecular Structure Elucidation (New Developments in NMR) by Mikhail E Elyashberg, Antony Williams and Kirill Blinov. Edited by William Price. RSC Publishing, 2012. ISBN: 978-1-84973-432-5; eISBN: 978-1-84973-457-8"
     *
     * @param ac             atom container to go through
     * @param connectionTree connection tree to expand, incl. the root node
     * @param queue          queue to use containing the atom index of the root node and start sphere
     * @param visited        atom indices which are already "visited" and
     *                       should be ignored
     * @param exclude        atom indices which to exclude from search
     * @param maxSphere      spherical limit, a null value means no limit
     */
    private static void BFS(final IAtomContainer ac, final ConnectionTree connectionTree, final Queue<int[]> queue,
                            final Set<Integer> visited, final Set<Integer> exclude, final Integer maxSphere) {
        // all nodes visited?
        if (queue.isEmpty()) {
            return;
        }
        final int[] queueValue = queue.remove();
        final int atomIndex = queueValue[0];
        final int sphere = queueValue[1];
        final IAtom atom = ac.getAtom(atomIndex);
        // mark atom as visited
        visited.add(atomIndex);

        IBond bond;
        // add nodes and bonds in lower spheres
        // go to all child nodes
        int connectedAtomIndex;
        for (final IAtom connectedAtom : ac.getConnectedAtomsList(atom)) {
            connectedAtomIndex = ac.indexOf(connectedAtom);
            bond = ac.getBond(atom, connectedAtom);
            // add children to queue if not already visited and connection is allowed or maxSphere is not reached yet
            if (!exclude.contains(connectedAtomIndex)) {
                if (keepBond(atom, connectedAtom, bond)
                        || maxSphere
                        == null
                        || sphere
                        < maxSphere) {
                    // add children to queue if not already visited and not already waiting in queue
                    if (!visited.contains(connectedAtomIndex)
                            && !queue.contains(connectedAtomIndex)) {
                        queue.add(new int[]{connectedAtomIndex, sphere
                                + 1});
                        connectionTree.addNode(connectedAtom, connectedAtomIndex, atomIndex, bond);
                    }
                }
            }
        }

        // further extension of connection tree
        BFS(ac, connectionTree, queue, visited, exclude, maxSphere);
    }

    private static boolean keepBond(final IAtom atom1, final IAtom atom2, final IBond bond) {
        // hetero-hetero or carbon-hetero
        if ((isHeteroAtom(atom1)
                && isHeteroAtom(atom2))
                || (isCarbonAtom(atom1)
                && isHeteroAtom(atom2))
                || (isHeteroAtom(atom1)
                && isCarbonAtom(atom2))) {
            return true;
        }

        //        // do not cut ring bonds
        //        if (bond.isInRing()) {
        //            return true;
        //        }

        // carbon-carbon or carbon-hetero with higher bond order
        if (
            //                ((isCarbonAtom(atom1)
            //                && isHeteroAtom(atom2))
            //                || (isHeteroAtom(atom1)
            //                && isCarbonAtom(atom2))
            //                || (isCarbonAtom(atom1)
            //                && isCarbonAtom(atom2)))
            //                &&
                bond.getOrder()
                    .numeric()
                        >= 3
            //                && !bond.isAromatic()
        ) {
            return true;
        }

        // one carbon has bonds to multiple hetero atoms
        if (isCarbonAtom(atom1)
                && isHeteroAtom(atom2)) {
            int heteroAtomCount = 0;
            for (final IAtom atom3 : atom1.getContainer()
                                          .getConnectedAtomsList(atom1)) {
                if (isHeteroAtom(atom3)) {
                    heteroAtomCount++;
                }
            }
            return heteroAtomCount
                    >= 2;
        } else if (isHeteroAtom(atom1)
                && isCarbonAtom(atom2)) {
            int heteroAtomCount = 0;
            for (final IAtom atom3 : atom2.getContainer()
                                          .getConnectedAtomsList(atom2)) {
                if (isHeteroAtom(atom3)) {
                    heteroAtomCount++;
                }
            }
            return heteroAtomCount
                    >= 2;
        }

        return false;
    }

    private static boolean isHeteroAtom(final IAtom atom) {
        return !atom.getSymbol()
                    .equals("H")
                && !isCarbonAtom(atom);
    }

    private static boolean isCarbonAtom(final IAtom atom) {
        return atom.getSymbol()
                   .equals("C");
    }

    public static List<Integer> buildFragmentAtomIndicesList(final IAtomContainer structure, final int rootAtomIndex,
                                                             final Integer maxSphere, final Set<Integer> exclude,
                                                             final boolean withPseudoAtoms) {
        return buildFragmentTree(structure, rootAtomIndex, maxSphere, exclude, withPseudoAtoms).getKeys();
    }
}
