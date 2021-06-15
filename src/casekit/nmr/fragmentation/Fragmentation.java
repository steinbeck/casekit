package casekit.nmr.fragmentation;

import casekit.nmr.fragmentation.model.ConnectionTree;
import casekit.nmr.fragmentation.model.ConnectionTreeNode;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.Bond;
import org.openscience.cdk.silent.PseudoAtom;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

import java.util.*;

public class Fragmentation {

    /**
     * Function for extending a given connection tree only containing
     * its root node (0th sphere) by means of Breadth-First-Search (BFS).
     * Until a certain maximum sphere, each reachable next neighbor atom
     * is stored in a parent-child-relationship.
     * In addition, bonds within rings or between hetero atoms will be kept.
     *
     * @param ac              atom container to go through
     * @param rootAtomIndex   root atom index to start from
     * @param maxSphere       spherical limit
     * @param exclude         atom indices which to exclude from search
     * @param withPseudoAtoms places pseudo atoms in the "outer" sphere
     *
     * @return connection tree
     */
    public static ConnectionTree BFS(final IAtomContainer ac, final int rootAtomIndex, final int maxSphere,
                                     final Set<Integer> exclude, final boolean withPseudoAtoms) {
        // create queue and connection tree for BFS
        final Queue<int[]> queue = new LinkedList<>();
        queue.add(new int[]{rootAtomIndex, 0});
        final ConnectionTree connectionTree = new ConnectionTree(ac.getAtom(rootAtomIndex), rootAtomIndex);

        BFS(ac, connectionTree, queue, new HashSet<>(), exclude, maxSphere, withPseudoAtoms);

        return connectionTree;
    }

    /**
     * Function for extending a given connection tree only containing
     * its root node (0th sphere) by means of Breadth-First-Search (BFS).
     * Until a certain maximum sphere, each reachable next neighbor atom
     * is stored in a parent-child-relationship.
     * In addition, bonds within rings or between hetero atoms will be kept.
     *
     * @param ac              atom container to go through
     * @param connectionTree  connection tree to expand, incl. the root node
     * @param queue           queue to use containing the atom index of the root node and start sphere
     * @param visited         atom indices which are already "visited" and
     *                        should be ignored
     * @param exclude         atom indices which to exclude from search
     * @param maxSphere       spherical limit
     * @param withPseudoAtoms places pseudo atoms in the "outer" sphere
     */
    private static void BFS(final IAtomContainer ac, final ConnectionTree connectionTree, final Queue<int[]> queue,
                            final Set<Integer> visited, final Set<Integer> exclude, final int maxSphere,
                            final boolean withPseudoAtoms) {
        // all nodes visited?
        if (queue.isEmpty()) {
            return;
        }
        final int[] queueValue = queue.remove();
        final int atomIndex = queueValue[0];
        final int sphere = queueValue[1];
        final IAtom atom = ac.getAtom(atomIndex);
        final ConnectionTreeNode node = connectionTree.getNode(atomIndex);
        // mark atom as visited
        visited.add(atomIndex);

        IBond bond;
        ConnectionTreeNode connectedAtomNode;
        // add nodes and bonds in lower spheres
        // go to all child nodes
        int connectedAtomIndex;
        for (final IAtom connectedAtom : ac.getConnectedAtomsList(atom)) {
            connectedAtomIndex = ac.indexOf(connectedAtom);
            bond = ac.getBond(atom, connectedAtom);
            // add children to queue if not already visited and connection is allowed or maxSphere is not reached yet
            if ((keepConnection(node.getAtom(), connectedAtom, bond)
                    || sphere
                    < maxSphere)
                    && !exclude.contains(connectedAtomIndex)) {
                // and not already waiting in queue
                if (!visited.contains(connectedAtomIndex)
                        && !queue.contains(connectedAtomIndex)) {
                    queue.add(new int[]{connectedAtomIndex, sphere
                            + 1});
                    connectionTree.addNode(connectedAtom, connectedAtomIndex, node.getKey(), bond);
                } else {
                    // node already exists in tree; add a further parent to connected atom (for ring closures)
                    connectedAtomNode = connectionTree.getNode(connectedAtomIndex);
                    if (connectedAtomNode
                            != null
                            && !ConnectionTree.hasRingClosureParent(node, connectedAtomNode)
                            && !ConnectionTree.hasRingClosureParent(connectedAtomNode, node)) {
                        connectionTree.addRingClosureNode(connectedAtomIndex, node.getKey(), bond);
                        connectionTree.addRingClosureNode(node.getKey(), connectedAtomIndex, bond);
                    }
                }
            } else if (withPseudoAtoms) {
                connectionTree.addNode(new PseudoAtom(connectedAtom), connectedAtomIndex, node.getKey(), bond);
            }
        }

        // further extension of connection tree
        BFS(ac, connectionTree, queue, visited, exclude, maxSphere, withPseudoAtoms);
    }

    public static boolean keepConnection(final IAtom atom1, final IAtom atom2, final IBond bond) {
        // hetero-hetero or carbon-hetero
        if ((isHeteroAtom(atom1)
                && isHeteroAtom(atom2))
                || (isCarbonAtom(atom1)
                && isHeteroAtom(atom2))
                || (isHeteroAtom(atom1)
                && isCarbonAtom(atom2))) {
            return true;
        }
        // do not cut ring bonds
        if (bond.isInRing()) {
            return true;
        }
        // carbon-carbon or carbon-hetero with higher bond order
        return ((isCarbonAtom(atom1)
                && isHeteroAtom(atom2))
                || (isHeteroAtom(atom1)
                && isCarbonAtom(atom2))
                || (isCarbonAtom(atom1)
                && isCarbonAtom(atom2)))
                && bond.getOrder()
                       .numeric()
                >= 2
                && !bond.isAromatic();

        //        // one carbon has bonds to multiple hetero atoms
        //        if (isCarbonAtom(atom1)
        //                && isHeteroAtom(atom2)) {
        //            int heteroAtomCount = 0;
        //            for (final IAtom atom3 : atom1.getContainer()
        //                                          .getConnectedAtomsList(atom1)) {
        //                if (isHeteroAtom(atom3)) {
        //                    heteroAtomCount++;
        //                }
        //            }
        //            if (heteroAtomCount
        //                    >= 2) {
        //                return true;
        //            }
        //        } else if (isHeteroAtom(atom1)
        //                && isCarbonAtom(atom2)) {
        //            int heteroAtomCount = 0;
        //            for (final IAtom atom3 : atom2.getContainer()
        //                                          .getConnectedAtomsList(atom2)) {
        //                if (isHeteroAtom(atom3)) {
        //                    heteroAtomCount++;
        //                }
        //            }
        //            if (heteroAtomCount
        //                    >= 2) {
        //                return true;
        //            }
        //        }
    }

    public static boolean isHeteroAtom(final IAtom atom) {
        return !atom.getSymbol()
                    .equals("H")
                && !isCarbonAtom(atom);
    }

    public static boolean isCarbonAtom(final IAtom atom) {
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
    public static IAtomContainer buildAtomContainer(final ConnectionTree connectionTree) {
        // create new atom container and add the connection trees structure, beginning at the root atom
        final IAtomContainer ac = SilentChemObjectBuilder.getInstance()
                                                         .newAtomContainer();
        addToAtomContainer(connectionTree, ac, null, null);

        return ac;
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
}
