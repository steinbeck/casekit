package casekit.nmr.fragmentation;

import casekit.nmr.fragmentation.model.ConnectionTree;
import casekit.nmr.fragmentation.model.ConnectionTreeNode;
import casekit.nmr.model.*;
import casekit.nmr.utils.Utils;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.silent.Bond;
import org.openscience.cdk.silent.PseudoAtom;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

import java.util.*;
import java.util.stream.Collectors;

public class Fragmentation {

    public static List<DataSet> buildFragments(final DataSet dataSet, final Integer maxSphere,
                                               final Integer maxSphereRing, final boolean withPseudoAtoms) {
        final List<DataSet> fragmentDataSetList = new ArrayList<>();
        final List<ConnectionTree> fragmentTrees = buildFragmentTrees(dataSet.getStructure()
                                                                             .toAtomContainer(), maxSphere,
                                                                      maxSphereRing, withPseudoAtoms);
        final IAtomContainer structure = dataSet.getStructure()
                                                .toAtomContainer();
        final String spectrumAtomType = Utils.getAtomTypeFromSpectrum(dataSet.getSpectrum(), 0);
        List<Integer> substructureAtomIndices;
        Spectrum subspectrum;
        Assignment subassignment;
        IAtom atomInStructure;
        Signal signal;
        DataSet subDataSet;
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
                atomInStructure = structure.getAtom(substructureAtomIndices.get(j));
                if (atomInStructure.getSymbol()
                                   .equals(spectrumAtomType)) {
                    if (dataSet.getAssignment()
                               .getIndex(0, substructureAtomIndices.get(j))
                            == null
                            || dataSet.getSpectrum()
                                      .getSignal(dataSet.getAssignment()
                                                        .getIndex(0, substructureAtomIndices.get(j)))
                            == null) {
                        return null;
                    }

                    signal = dataSet.getSpectrum()
                                    .getSignal(dataSet.getAssignment()
                                                      .getIndex(0, substructureAtomIndices.get(j)));
                    if (signal
                            != null) {
                        signal = signal.buildClone();
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
            subspectrum.setSolvent(dataSet.getSpectrum()
                                          .getSolvent());
            subspectrum.setSpectrometerFrequency(dataSet.getSpectrum()
                                                        .getSpectrometerFrequency());

            subDataSet = new DataSet();
            subDataSet.setStructure(new ExtendedConnectionMatrix(toAtomContainer(fragmentTree)));
            subDataSet.setSpectrum(subspectrum);
            subDataSet.setAssignment(subassignment);

            fragmentDataSetList.add(subDataSet);
        }

        return fragmentDataSetList;
    }

    public static List<IAtomContainer> buildFragments(final IAtomContainer structure, final Integer maxSphere,
                                                      final Integer maxSphereRing, final boolean withPseudoAtoms) {
        final List<ConnectionTree> fragmentTrees = buildFragmentTrees(structure, maxSphere, maxSphereRing,
                                                                      withPseudoAtoms);
        return fragmentTrees.stream()
                            .map(Fragmentation::toAtomContainer)
                            .collect(Collectors.toList());
    }

    public static List<ConnectionTree> buildFragmentTrees(final IAtomContainer structure, final Integer maxSphere,
                                                          final Integer maxSphereRing, final boolean withPseudoAtoms) {
        final List<ConnectionTree> fragments = new ArrayList<>();
        try {
            // build fragments from detected rings and extend by given maximum sphere for rings
            ConnectionTree connectionTreeRing, connectionTreeOuterSphere, subtreeToAdd;
            final IRingSet ringSet = Cycles.all(structure)//essential(structure)
                                           .toRingSet();
            List<Integer> atomIndicesInRing;
            Set<Integer> atomIndicesOutOfRing;
            IAtomContainer ringAtomContainer;
            for (int i = 0; i
                    < ringSet.getAtomContainerCount(); i++) {
                ringAtomContainer = ringSet.getAtomContainer(i);
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
                        if (!addToConnectionTree(connectionTreeRing, connectionTreeOuterSphere.getRootNode()
                                                                                              .getKey(), subtreeToAdd,
                                                 connectionTreeOuterSphere.getBond(
                                                         connectionTreeOuterSphere.getRootNode()
                                                                                  .getKey(), key))) {
                            continue;
                        }
                        atomIndicesInRing.addAll(subtreeToAdd.getKeys());
                    }
                }
                // close rings
                closeRings(connectionTreeRing, structure);
                // attach pseudo atoms if desired
                if (withPseudoAtoms) {
                    attachPseudoAtoms(connectionTreeRing, structure);
                }
                fragments.add(connectionTreeRing);
            }
            // build fragment for each non-ring atom
            for (int i = 0; i
                    < structure.getAtomCount(); i++) {
                fragments.add(
                        Fragmentation.buildFragmentTree(structure, i, maxSphere, new HashSet<>(), withPseudoAtoms));
            }
        } catch (final CDKException e) {
            e.printStackTrace();
        }

        return fragments;
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
     */
    public static IAtomContainer buildFragment(final IAtomContainer ac, final int rootAtomIndex,
                                               final Integer maxSphere, final Set<Integer> exclude,
                                               final boolean withPseudoAtoms) {
        return toAtomContainer(buildFragmentTree(ac, rootAtomIndex, maxSphere, exclude, withPseudoAtoms));
    }

    /**
     * Function for extending a given connection tree only containing
     * its root node (0th sphere) by means of Breadth-First-Search (BFS).
     * Until a certain maximum sphere, each reachable next neighbor atom
     * is stored in a parent-child-relationship.
     * In addition, bonds within rings or between hetero atoms will be kept.
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
        closeRings(connectionTree, structure);
        // add pseudo atoms
        if (withPseudoAtoms) {
            attachPseudoAtoms(connectionTree, structure);
        }

        return connectionTree;
    }

    public static void closeRings(final ConnectionTree connectionTree, final IAtomContainer structure) {
        // close rings
        IBond bond;
        final int maxSphereTree = connectionTree.getMaxSphere(false);
        for (int s = 0; s
                <= maxSphereTree; s++) {
            for (final ConnectionTreeNode nodeInSphere1 : connectionTree.getNodesInSphere(s, false)) {
                // set connections (parent nodes) in sphere nodes which have to be connected -> ring closures
                for (int s2 = s; s2
                        <= maxSphereTree; s2++) {
                    for (final ConnectionTreeNode nodeInSphere2 : connectionTree.getNodesInSphere(s2, false)) {
                        if ((structure.getBond(nodeInSphere1.getAtom(), nodeInSphere2.getAtom())
                                != null)
                                && !ConnectionTree.nodesFormRingClosure(nodeInSphere1, nodeInSphere2)) {
                            bond = structure.getBond(nodeInSphere1.getAtom(), nodeInSphere2.getAtom());
                            connectionTree.addRingClosureNode(nodeInSphere1.getKey(), nodeInSphere2.getKey(), bond);
                            connectionTree.addRingClosureNode(nodeInSphere2.getKey(), nodeInSphere1.getKey(), bond);
                        }
                    }
                }
            }
        }
    }

    public static void attachPseudoAtoms(final ConnectionTree connectionTree, final IAtomContainer structure) {
        int atomIndexInStructure;
        for (final ConnectionTreeNode node : connectionTree.getNodes(false)) {
            for (final IAtom connectedAtom : structure.getConnectedAtomsList(node.getAtom())) {
                atomIndexInStructure = structure.indexOf(connectedAtom);
                if (connectionTree.getBond(node.getKey(), atomIndexInStructure)
                        == null
                        && connectionTree.getBond(atomIndexInStructure, node.getKey())
                        == null) {
                    addPseudoNode(connectionTree, structure.getAtomCount()
                                          + connectionTree.getNodesCount(false), node.getKey(),
                                  structure.getBond(node.getAtom(), connectedAtom));
                }
            }
        }
    }

    /**
     * Function for extending a given connection tree only containing
     * its root node (0th sphere) by means of Breadth-First-Search (BFS).
     * Until a certain maximum sphere, each reachable next neighbor atom
     * is stored in a parent-child-relationship.
     * In addition, bonds within rings or between hetero atoms will be kept.
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

    private static boolean addPseudoNode(final ConnectionTree connectionTree, final int pseudoNodeKey,
                                         final int parentNodeKey, final IBond bondToParent) {
        if (!connectionTree.addNode(new PseudoAtom("R"), pseudoNodeKey, parentNodeKey, bondToParent)) {
            return false;
        }
        final ConnectionTreeNode pseudoNode = connectionTree.getNode(pseudoNodeKey);
        pseudoNode.setIsPseudoNode(true);
        //                pseudoNode.getAtom()
        //                          .setImplicitHydrogenCount(connectedAtom.getImplicitHydrogenCount());

        return true;
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

    /**
     * Reconstructs a structure from a given connection tree,
     * including ring closures.
     *
     * @param connectionTree connection tree
     *
     * @return IAtomContainer
     */
    public static IAtomContainer toAtomContainer(final ConnectionTree connectionTree) {
        // create new atom container and add the connection trees structure, beginning at the root atom
        final IAtomContainer ac = SilentChemObjectBuilder.getInstance()
                                                         .newAtomContainer();
        addToAtomContainer(connectionTree, ac, null, null);

        return ac;
    }

    /**
     * Adds a subtree to a node in another connection tree.
     *
     * @param connectionTree connection tree
     * @param parentNodeKey  parent node key in connection tree
     * @param subtree        subtree to add
     * @param bondToLink     bond
     */
    public static boolean addToConnectionTree(final ConnectionTree connectionTree, final int parentNodeKey,
                                              final ConnectionTree subtree, final IBond bondToLink) {
        return ConnectionTree.addSubtree(connectionTree, parentNodeKey, subtree, bondToLink);
    }

    /**
     * Adds the substructure of a connection tree to an atom container. <br>
     * The substructure can be linked via a bond and an atom index in the container, but this is optional.
     * If both, the bond and atom index to link, are not given (null) then the substructure will just be added
     * to the atom container without linkage.
     *
     * @param connectionTree
     * @param ac
     * @param atomIndexInStructureToLink
     * @param bondToLink
     */
    public static void addToAtomContainer(final ConnectionTree connectionTree, final IAtomContainer ac,
                                          final Integer atomIndexInStructureToLink, final IBond bondToLink) {
        List<ConnectionTreeNode> nodesInSphere;
        ConnectionTreeNode nodeInSphere, parentNode, partnerNode;
        IBond bond, bondToParent;
        // add root atom to given atom container and link it via a given linking bond
        ac.addAtom(connectionTree.getRootNode()
                                 .getAtom());
        if ((atomIndexInStructureToLink
                != null)
                && (bondToLink
                != null)) {
            final IBond bondToAdd = new Bond(ac.getAtom(atomIndexInStructureToLink), ac.getAtom(ac.getAtomCount()
                                                                                                        - 1));
            bondToAdd.setOrder(bondToLink.getOrder());
            bondToAdd.setIsInRing(bondToLink.isInRing());
            bondToAdd.setIsAromatic(bondToLink.isAromatic());
            bondToAdd.setAtom(ac.getAtom(atomIndexInStructureToLink), 0);
            bondToAdd.setAtom(ac.getAtom(ac.getAtomCount()
                                                 - 1), 1);
            ac.addBond(bondToAdd);
        }
        // for each sphere: add the atom which is stored as node to atom container and set bonds between parent nodes
        for (int s = 1; s
                <= connectionTree.getMaxSphere(false); s++) {
            // first add all atoms and its parents (previous sphere only, incl. bonds) to structure
            nodesInSphere = connectionTree.getNodesInSphere(s, false);
            for (int i = 0; i
                    < nodesInSphere.size(); i++) {
                nodeInSphere = nodesInSphere.get(i);
                if (nodeInSphere.isRingClosureNode()) {
                    continue;
                }
                ac.addAtom(nodeInSphere.getAtom());
                parentNode = nodeInSphere.getParent();
                bondToParent = nodeInSphere.getBondToParent();
                bond = new Bond(nodeInSphere.getAtom(), parentNode.getAtom(), bondToParent.getOrder());
                bond.setIsInRing(bondToParent.isInRing());
                bond.setIsAromatic(bondToParent.isAromatic());
                ac.addBond(bond);
            }
        }
        for (int s = 1; s
                <= connectionTree.getMaxSphere(true); s++) {
            // and as second add the remaining bonds (ring closures) to structure
            nodesInSphere = connectionTree.getNodesInSphere(s, true);
            for (int i = 0; i
                    < nodesInSphere.size(); i++) {
                nodeInSphere = nodesInSphere.get(i);
                if (!nodeInSphere.isRingClosureNode()) {
                    continue;
                }
                parentNode = nodeInSphere.getParent();
                partnerNode = nodeInSphere.getRingClosureParent();
                if (ac.getBond(ac.getAtom(ac.indexOf(partnerNode.getAtom())),
                               ac.getAtom(ac.indexOf(parentNode.getAtom())))
                        == null) {
                    bondToParent = nodeInSphere.getBondToParent();
                    bond = new Bond(parentNode.getAtom(), partnerNode.getAtom(), bondToParent.getOrder());
                    bond.setIsInRing(bondToParent.isInRing());
                    bond.setIsAromatic(bondToParent.isAromatic());
                    ac.addBond(bond);
                }
            }
        }
    }

    public static List<Integer> buildFragmentAtomIndicesList(final IAtomContainer structure, final int rootAtomIndex,
                                                             final Integer maxSphere, final Set<Integer> exclude,
                                                             final boolean withPseudoAtoms) {
        return buildFragmentTree(structure, rootAtomIndex, maxSphere, exclude, withPseudoAtoms).getKeys();
    }
}
