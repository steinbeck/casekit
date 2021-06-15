/*
 * The MIT License
 *
 * Copyright (c) 2019 Michael Wenk [https://github.com/michaelwenk]
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
package casekit.nmr.hose;


import casekit.nmr.fragmentation.Fragmentation;
import casekit.nmr.fragmentation.model.ConnectionTree;
import casekit.nmr.fragmentation.model.ConnectionTreeNode;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.Atom;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

import java.util.*;

/**
 * Class to build HOSE code strings from molecules and vice versa
 * by using connection trees as intermediate forms.
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class HOSECodeBuilder {

    /**
     * Creates a partial sphere string content from the children of a given parent node.
     *
     * @param nodeInPrevSphere          parent node to create a partial sphere string content from
     * @param useBremserElementNotation whether to use Bremser notation
     *
     * @return
     *
     * @throws CDKException
     */
    private static List<String> buildPositionsInSphere(final ConnectionTreeNode nodeInPrevSphere,
                                                       final boolean useBremserElementNotation) throws CDKException {
        final List<ConnectionTreeNode> nodesInSphere = nodeInPrevSphere.getChildNodes();
        final List<String> positionsInSphere = new ArrayList<>();
        ConnectionTreeNode nodeInSphere;
        IBond bond;
        String position;
        for (int j = 0; j
                < nodesInSphere.size(); j++) {
            nodeInSphere = nodesInSphere.get(j);
            bond = nodeInPrevSphere.getBondsToChildren()
                                   .get(j);
            position = "";
            if (Utils.getSymbolForBond(bond)
                    == null) {
                throw new CDKException(Thread.currentThread()
                                             .getStackTrace()[1].getMethodName()
                                               + ": no bond information");
            }
            position += Utils.getSymbolForBond(bond);
            if (nodeInSphere.isRingClosureNode()) {
                position += "&";
            } else {
                if (useBremserElementNotation) {
                    position += Utils.toHOSECodeSymbol(nodeInSphere.getAtom()
                                                                   .getSymbol());
                } else {
                    position += nodeInSphere.getAtom()
                                            .getSymbol();
                }
                //                if(nodeInSphere.getAtom().getImplicitHydrogenCount() != null){
                //                    position += "[" + nodeInSphere.getAtom().getImplicitHydrogenCount() + "]";
                //                }
                position += buildFormalChargeCode(nodeInSphere.getAtom());
            }
            positionsInSphere.add(position);
        }

        return positionsInSphere;
    }

    /**
     * Builds the content of a sphere of the HOSE code which is to generate.
     *
     * @param connectionTree            connection tree to use
     * @param sphere                    sphere to selected from connection tree
     * @param delimiter                 sphere's delimiter
     * @param useBremserElementNotation whether to use Bremser notation
     *
     * @return
     *
     * @throws CDKException
     */
    private static String buildSphereString(final ConnectionTree connectionTree, final int sphere,
                                            final String delimiter,
                                            final boolean useBremserElementNotation) throws CDKException {
        StringBuilder sphereString = new StringBuilder();
        final List<ConnectionTreeNode> nodesInPrevSphere = connectionTree.getNodesInSphere(sphere
                                                                                                   - 1, true);
        ConnectionTreeNode nodeInPrevSphere;
        // for all nodes in previous sphere
        for (int i = 0; i
                < nodesInPrevSphere.size(); i++) {
            nodeInPrevSphere = nodesInPrevSphere.get(i);
            // skip ring closure nodes
            if (nodeInPrevSphere.isRingClosureNode()) {
                if ((i
                        == nodesInPrevSphere.size()
                        - 1)
                        && sphereString.toString()
                                       .endsWith(",")) {
                    sphereString = new StringBuilder(sphereString.substring(0, sphereString.length()
                            - 1));
                }
                continue;
            }
            // for all child nodes in the requested sphere
            if (nodeInPrevSphere.hasChildren()) {
                for (final String position : buildPositionsInSphere(nodeInPrevSphere, useBremserElementNotation)) {
                    sphereString.append(position);
                }
            }
            // add delimiter
            if (i
                    < nodesInPrevSphere.size()
                    - 1) {
                sphereString.append(delimiter);
            }
        }

        return sphereString.toString();
    }

    private static String buildFormalChargeCode(final IAtom atom) {
        if ((atom
                == null)
                || (atom.getFormalCharge()
                == null)
                || (atom.getFormalCharge()
                == 0)) {
            return "";
        }
        final String sign = atom.getFormalCharge()
                                    < 0
                            ? "-"
                            : "+";

        return Math.abs(atom.getFormalCharge())
                       == 1
               ? sign
               : "'"
                       + sign
                       + Math.abs(atom.getFormalCharge())
                       + "'";
    }

    /**
     * Actual function to build a HOSE code.
     *
     * @param connectionTree            connection tree to use
     * @param useBremserElementNotation whether to use Bremser notation
     *
     * @return
     *
     * @throws CDKException
     */
    private static String buildHOSECodeString(final ConnectionTree connectionTree,
                                              final boolean useBremserElementNotation) throws CDKException {
        final IAtom rootAtom = connectionTree.getRootNode()
                                             .getAtom();
        final int maxSphere = connectionTree.getMaxSphere(true);
        // zeroth sphere
        final StringBuilder HOSECode = new StringBuilder(rootAtom.getSymbol()
                                                                 + "-"
                                                                 + (rootAtom.getBondCount()
                + (rootAtom.getImplicitHydrogenCount()
                           == null
                   ? 0
                   : rootAtom.getImplicitHydrogenCount()))
                                                                 + buildFormalChargeCode(rootAtom));
        HOSECode.append(";");
        String delimiter;
        // go through each sphere of the connection tree
        for (int s = 1; s
                <= maxSphere; s++) {
            if (s
                    == 1) {
                delimiter = "";
            } else {
                delimiter = ",";
            }
            // create sphere string and add it to HOSE code string
            HOSECode.append(buildSphereString(connectionTree, s, delimiter, useBremserElementNotation));
            if (s
                    == 1) {
                HOSECode.append("(");
            }
            if (s
                    > 1
                    && s
                    < maxSphere) {
                HOSECode.append("/");
            }
        }
        if (maxSphere
                == 0) {
            HOSECode.append("(");
        }
        HOSECode.append(")");

        return HOSECode.toString();
    }

    public static String buildHOSECode(final ConnectionTree connectionTree,
                                       final boolean useBremserElementNotation) throws CDKException {
        return buildHOSECodeString(connectionTree, useBremserElementNotation);
    }

    public static String buildHOSECode(final IAtomContainer ac, final int rootAtomIndex, final Integer maxSphere,
                                       final boolean useBremserElementNotation) throws CDKException {
        return HOSECodeBuilder.buildHOSECode(HOSECodeBuilder.buildConnectionTree(ac, rootAtomIndex, maxSphere),
                                             useBremserElementNotation);
    }

    /**
     * Builds a connection tree of an atom container with specific start atom
     * and maximum number of spheres.
     * If the atoms in the atom container are not fully connected, then the
     * connection tree will be built until the last atom of all connected atoms
     * to the start atom is reached.
     *
     * @param ac            atom container
     * @param rootAtomIndex starting atom
     * @param maxSphere     if this is set to null, then the connection tree of whole
     *                      structure will be created
     *
     * @return
     *
     * @see ConnectionTree
     */
    public static ConnectionTree buildConnectionTree(final IAtomContainer ac, final int rootAtomIndex,
                                                     final Integer maxSphere) {
        return HOSECodeBuilder.buildConnectionTree(ac, rootAtomIndex, maxSphere, new HashSet<>());
    }

    /**
     * Builds a connection tree of an atom container with specific start atom
     * and maximum number of spheres.
     * If the atoms in the atom container are not fully connected, then the
     * connection tree will be built until the last atom of all connected atoms
     * to the start atom is reached.
     *
     * @param ac            atom container
     * @param rootAtomIndex starting atom
     * @param maxSphere     if this is set to null, then the connection tree of whole
     *                      structure will be created
     * @param visited       certain atom indices can be given here to ignore atoms
     *                      in BFS; they are then seen as already visited and not included in
     *                      the connection tree
     *
     * @return
     *
     * @see ConnectionTree
     */
    public static ConnectionTree buildConnectionTree(final IAtomContainer ac, final int rootAtomIndex,
                                                     final Integer maxSphere, final Set<Integer> visited) {
        // create queue for BFS and add root atom index
        final Queue<Integer> queue = new LinkedList<>();
        queue.add(rootAtomIndex);
        final ConnectionTree connectionTree = new ConnectionTree(ac.getAtom(rootAtomIndex), rootAtomIndex);
        BFS(ac, connectionTree, queue, new HashSet<>(visited), maxSphere);

        Utils.rankChildNodes(connectionTree);

        return connectionTree;
    }

    /**
     * Builds a connection tree from a given HOSE code. <br>
     * IMPORTANT: At the moment, ring closures can not be restored
     * from an HOSE code because of ambiguities.
     * So only a structural skeleton will be generated.
     *
     * @param HOSECode                  HOSE code
     * @param useBremserElementNotation whether given HOSE code contains Bremser notation
     *
     * @return
     *
     * @throws CDKException
     */
    public static ConnectionTree buildConnectionTree(final String HOSECode,
                                                     final boolean useBremserElementNotation) throws CDKException {
        final Map<Integer, ArrayList<Object[]>> ringClosures = new HashMap<>();
        final List<String> sphereStrings = Utils.splitHOSECodeIntoSpheres(HOSECode);
        IAtom atom;
        IBond bond;
        final int maxSphere;
        String bondTypeString, atomTypeString, childElementCore;
        Map<Integer, ArrayList<String>> positionsInSphere;
        ConnectionTreeNode parentNodeInPrevSphere;
        // set maxSphere
        maxSphere = sphereStrings.size()
                - 1;
        // zeroth sphere
        positionsInSphere = Utils.splitHOSECodeSphereIntoPositions(sphereStrings.get(0), true);
        // create root atom
        atom = new Atom(positionsInSphere.get(0)
                                         .get(0));
        // add charge to root atom
        atom.setFormalCharge(Integer.parseInt(positionsInSphere.get(0)
                                                               .get(1)));
        final ConnectionTree connectionTree = new ConnectionTree(atom, 0);
        // higher spheres
        for (int sphere = 1; sphere
                <= maxSphere; sphere++) {
            // get positions (sections separated by comma) of current sphere
            positionsInSphere = Utils.splitHOSECodeSphereIntoPositions(sphereStrings.get(sphere), false);
            // for all positions
            for (final int positionIndex : positionsInSphere.keySet()) {
                // for each child elements (symbols) in position
                for (final String childElement : positionsInSphere.get(positionIndex)) {
                    // ignore children containing null value from previous nodes; previous node has no further (unvisited in BFS) connected atoms
                    if (childElement
                            == null) {
                        continue;
                    }
                    bondTypeString = "";
                    atomTypeString = "";
                    childElementCore = childElement;
                    childElementCore = childElementCore.replace("+", "");
                    childElementCore = childElementCore.replace("-", "");
                    childElementCore = childElementCore.replaceAll("\\d", "");
                    // add new node and set bond to parent node or set a ring closure
                    if (childElementCore.contains("&")) { // ring closure
                        if (childElementCore.length()
                                == 2) {
                            bondTypeString = String.valueOf(childElementCore.charAt(0));
                        }
                        // the parent node/atom in previous sphere and its key we already have of a ring closure;
                        // the bond information we already have too (see below)
                        parentNodeInPrevSphere = connectionTree.getNodesInSphere(sphere
                                                                                         - 1, true)
                                                               .get(positionIndex);

                        if (!ringClosures.containsKey(sphere)) {
                            ringClosures.put(sphere, new ArrayList<>());
                        }
                        bond = SilentChemObjectBuilder.getInstance()
                                                      .newBond();
                        bond.setOrder(Utils.getBondOrderForSymbol(bondTypeString));
                        if (bondTypeString.equals("*")) {
                            bond.setIsInRing(true);
                            bond.setIsAromatic(true);
                        } else {
                            bond.setIsAromatic(false);
                        }
                        // store the ring closures and use them after looking at the HOSE code string
                        ringClosures.get(sphere)
                                    .add(new Object[]{parentNodeInPrevSphere, bond});


                        //                        // check whether the node in previous sphere is already involved in a ring closure; that should be not valid
                        //                        if(ConnectionTree.isAtRingClosure(parentNodeInPrevSphere)){
                        //                            continue;
                        //                        }
                        //
                        //                        // TODO: what we still not can detect for sure is the correct second node/atom of a ring closure
                        //                        ConnectionTreeNode parentNodeInSphere = null; // null is just a dummy value and should be replaced by the correct ConnectionTreeNode object
                        //
                        //                        // TODO: check that the detected node in sphere is not null; could be removed after the implementation of detection of that node
                        //                        if(parentNodeInSphere == null){
                        //                            continue;
                        //                        }
                        //                        // after that both node detections, check if that ring closure was already set beforehand by the reversed node order case
                        //                        if(ConnectionTree.nodesFormRingClosure(parentNodeInPrevSphere, parentNodeInSphere)){
                        //                            continue;
                        //                        }
                        //                        // otherwise build a new bond and fill it with
                        //                        bond = SilentChemObjectBuilder.getInstance().newBond();
                        //                        bond.setAtom(parentNodeInPrevSphere.getAtom(), 0);
                        //                        bond.setAtom(parentNodeInSphere.getAtom(), 1);
                        //                        bond.setOrder(Utils.getBondOrderForSymbol(bondTypeString));
                        //                        if (bondTypeString.equals("*")) {
                        //                            bond.setIsAromatic(true);
                        //                        } else {
                        //                            bond.setIsAromatic(false);
                        //                        }
                        //                        // set parent nodes as parents to each other, that one can detect them as ring closure afterwards
                        //                        parentNodeInPrevSphere.addParentNode(parentNodeInSphere, bond);
                        //                        connectionTree.addNode(null, -1 * parentNodeInPrevSphere.getKey(), parentNodeInPrevSphere.getKey(), bond, sphere + 1, true);
                        //                        parentNodeInSphere.addParentNode(parentNodeInPrevSphere, bond);
                        //                        connectionTree.addNode(null, -1 * parentNodeInSphere.getKey(), parentNodeInSphere.getKey(), bond, sphere + 1, true);

                    } else if (Utils.countAtoms(childElementCore)
                            == 1) { // each position contains either ring closures (&) or one element (e.g. C, Br), plus the bond information
                        if (childElementCore.length()
                                == 3) { // in case of bond type and an element with two letters, e.g. *Cl or =Br
                            bondTypeString = String.valueOf(childElementCore.charAt(0));
                            atomTypeString = String.valueOf(childElementCore.charAt(1));
                            atomTypeString += String.valueOf(childElementCore.charAt(2));
                        } else if (childElementCore.length()
                                == 2) { // in case of bond type and an element with one letter or an element with two letters, e.g. Cl or =N
                            if (Character.isLetter(childElementCore.charAt(0))) {
                                atomTypeString = String.valueOf(childElementCore.charAt(0));
                            } else {
                                bondTypeString = String.valueOf(childElementCore.charAt(0));
                            }
                            atomTypeString += String.valueOf(childElementCore.charAt(1));
                        } else if (childElementCore.length()
                                == 1) { // in case of an element with only one letter
                            atomTypeString = String.valueOf(childElementCore.charAt(0));
                        }
                        // there has to be some information (at least an element)
                        if (atomTypeString.isEmpty()) {
                            throw new CDKException(Thread.currentThread()
                                                         .getStackTrace()[1].getMethodName()
                                                           + ": no atom information in child element");
                        }
                        // otherwise set a new bond
                        bond = SilentChemObjectBuilder.getInstance()
                                                      .newBond();
                        if (useBremserElementNotation) {
                            atomTypeString = Utils.toElementSymbol(atomTypeString);
                        }
                        atom = new Atom(atomTypeString);
                        bond.setAtom(atom, 0);
                        bond.setAtom(connectionTree.getNodesInSphere(sphere
                                                                             - 1, true)
                                                   .get(positionIndex)
                                                   .getAtom(), 1);
                        bond.setOrder(Utils.getBondOrderForSymbol(bondTypeString));
                        bond.setIsAromatic(bondTypeString.equals("*"));
                        // set formal charge to atom
                        if (childElement.contains("-")) {
                            final String[] splitAtSign = childElement.split("-");
                            if (splitAtSign.length
                                    == 1) {
                                atom.setFormalCharge(-1);
                            } else {
                                atom.setFormalCharge(-1
                                                             * Integer.parseInt(splitAtSign[1]));
                            }
                        }
                        if (childElement.contains("+")) {
                            final String[] splitAtSign = childElement.split("\\+");
                            if (splitAtSign.length
                                    == 1) {
                                atom.setFormalCharge(1);
                            } else {
                                atom.setFormalCharge(Integer.parseInt(splitAtSign[1]));
                            }
                        }
                        // create a new node with the new build bond information to its parent node in the connection tree
                        connectionTree.addNode(atom, connectionTree.getNodesCount(false),
                                               connectionTree.getNodesInSphere(sphere
                                                                                       - 1, true)
                                                             .get(positionIndex)
                                                             .getKey(), bond);
                    } else {
                        throw new CDKException(Thread.currentThread()
                                                     .getStackTrace()[1].getMethodName()
                                                       + ": no valid components in child element");
                    }
                }
            }
        }

        // @TODO after storing the ring closures, try to close the rings
        //        System.out.println(" -> spheres count with ring closures: " + ringClosures.size());
        //        for (final int sphere : ringClosures.keySet()){
        //            System.out.println(" -> number of ring closures in sphere: " + sphere + " -> " + ringClosures.get(sphere).size());
        //        }

        Utils.rankChildNodes(connectionTree);

        return connectionTree;
    }

    /**
     * Function for extending a given connection tree only containing
     * its root node (0th sphere) by means of Breadth-First-Search (BFS).
     * Until a certain maximum sphere, each reachable next neighbor atom
     * is stored in a parent-child-relationship.
     *
     * @param ac             atom container to go through
     * @param connectionTree connection tree to expand, incl. the root node
     * @param queue          queue to use containing the atom index of the root node
     * @param visited        optional: atom indices which are already "visited" and
     *                       should be ignored
     * @param maxSphere      maximum number of spheres for connection tree extension
     */
    private static void BFS(final IAtomContainer ac, final ConnectionTree connectionTree, final Queue<Integer> queue,
                            final Set<Integer> visited, final Integer maxSphere) {
        // all nodes visited?
        if (queue.isEmpty()) {
            return;
        }
        final int atomIndex = queue.remove();
        final IAtom atom = ac.getAtom(atomIndex);
        final ConnectionTreeNode node = connectionTree.getNode(atomIndex);
        final int sphere = node.getSphere();
        // check whether the current sphere is to high, if maxSphere parameter is set
        if ((maxSphere
                != null)
                && (sphere
                > maxSphere)) {
            return;
        }
        // mark atom as visited
        visited.add(atomIndex);

        IBond bond;
        ConnectionTreeNode connectedAtomNode;
        if ((maxSphere
                != null)
                && (sphere
                == maxSphere)) {
            // set connections (parent nodes) in last sphere nodes which have to be connected -> ring closures
            // only parent nodes will be set to detect those ring closures again
            for (final ConnectionTreeNode nodeInLastSphere : connectionTree.getNodesInSphere(maxSphere, false)) {
                if ((ac.getBond(atom, nodeInLastSphere.getAtom())
                        != null)
                        && !ConnectionTree.hasRingClosureParent(node, nodeInLastSphere)
                        && !ConnectionTree.hasRingClosureParent(nodeInLastSphere, node)) {
                    bond = ac.getBond(node.getAtom(), nodeInLastSphere.getAtom());
                    connectionTree.addRingClosureNode(node.getKey(), nodeInLastSphere.getKey(), bond);
                    connectionTree.addRingClosureNode(nodeInLastSphere.getKey(), node.getKey(), bond);
                }
            }
        } else {
            // add nodes and bonds in lower spheres
            // go to all child nodes
            int connectedAtomIndex;
            for (final IAtom connectedAtom : ac.getConnectedAtomsList(atom)) {
                connectedAtomIndex = ac.indexOf(connectedAtom);
                bond = ac.getBond(atom, connectedAtom);
                // add children to queue if not already visited
                if (!visited.contains(connectedAtomIndex)) {
                    // and not already waiting in queue
                    if (!queue.contains(connectedAtomIndex)) {
                        queue.add(connectedAtomIndex);
                        connectionTree.addNode(connectedAtom, connectedAtomIndex, node.getKey(), bond);
                    } else {
                        // node already exists in tree; add a further parent to connected atom (for ring closures)
                        connectedAtomNode = connectionTree.getNode(connectedAtomIndex);
                        if (!ConnectionTree.hasRingClosureParent(node, connectedAtomNode)
                                && !ConnectionTree.hasRingClosureParent(connectedAtomNode, node)) {
                            connectionTree.addRingClosureNode(connectedAtomIndex, node.getKey(), bond);
                            connectionTree.addRingClosureNode(node.getKey(), connectedAtomIndex, bond);
                        }

                    }
                }
            }
        }
        // further extension of connectivity tree
        BFS(ac, connectionTree, queue, visited, maxSphere);
    }

    /**
     * Reconstructs a structure from a given HOSE code string. <br>
     * IMPORTANT: Ring closures are not restored, see
     * {@link #buildConnectionTree(String, boolean)}.
     *
     * @param HOSECode                  HOSE code
     * @param useBremserElementNotation whether the HOSE code includes Bremser notation
     *
     * @return IAtomContainer
     *
     * @see #buildConnectionTree(String, boolean)
     * @see Fragmentation#toAtomContainer(ConnectionTree)
     */
    public static IAtomContainer buildAtomContainer(final String HOSECode,
                                                    final boolean useBremserElementNotation) throws CDKException {
        return Fragmentation.toAtomContainer(HOSECodeBuilder.buildConnectionTree(HOSECode, useBremserElementNotation));
    }
}