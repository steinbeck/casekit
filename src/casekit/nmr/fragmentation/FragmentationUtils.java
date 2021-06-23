package casekit.nmr.fragmentation;

import casekit.nmr.fragmentation.model.ConnectionTree;
import casekit.nmr.fragmentation.model.ConnectionTreeNode;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.Bond;
import org.openscience.cdk.silent.PseudoAtom;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class FragmentationUtils {

    public static boolean adjustNodeKeys(final ConnectionTree fragmentTree, final IAtomContainer structure) {
        int atomIndex;
        for (final ConnectionTreeNode node : fragmentTree.getNodes(false)) {
            atomIndex = structure.indexOf(node.getAtom());
            if (atomIndex
                    < 0) {
                return false;
            }
            node.setKey(atomIndex);
        }
        fragmentTree.initKeySet();

        return true;
    }

    public static void removeDuplicates(final List<ConnectionTree> fragmentTrees) {
        final List<Set<Integer>> keySets = new ArrayList<>();
        final List<ConnectionTree> fragmentsToRemove = new ArrayList<>();
        for (final ConnectionTree fragment : fragmentTrees) {
            // ignore pseudo nodes
            final Set<Integer> keySet = fragment.getNodes(false)
                                                .stream()
                                                .filter(node -> !node.isPseudoNode())
                                                .map(ConnectionTreeNode::getKey)
                                                .collect(Collectors.toSet());
            if (keySets.stream()
                       .noneMatch(keySetTemp -> keySetTemp.size()
                               == keySet.size()
                               && keySetTemp.containsAll(keySet))) {
                keySets.add(keySet);
            } else {
                fragmentsToRemove.add(fragment);
            }
        }
        fragmentTrees.removeAll(fragmentsToRemove);
    }

    public static IAtomContainer closeRings(final IAtomContainer substructure, final IAtomContainer structure) {
        final ConnectionTree fragmentTree = Fragmentation.buildFragmentTree(substructure, 0, null, new HashSet<>(),
                                                                            false);
        closeRings(fragmentTree, structure);

        return toAtomContainer(fragmentTree);
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

    private static boolean addPseudoNode(final ConnectionTree connectionTree, final int pseudoNodeKey,
                                         final int parentNodeKey, final IBond bondToParent) {
        if (!connectionTree.addNode(new PseudoAtom("R"), pseudoNodeKey, parentNodeKey, bondToParent)) {
            return false;
        }
        final ConnectionTreeNode pseudoNode = connectionTree.getNode(pseudoNodeKey);
        pseudoNode.setIsPseudoNode(true);

        return true;
    }
}
