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

import casekit.nmr.hose.model.ConnectionTree;
import casekit.nmr.hose.model.ConnectionTreeNode;
import org.openscience.cdk.interfaces.IBond;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Utils {

    /**
     * Returns the summed subtree weight starting at a specific node in a connection
     * tree. The weight of starting node is included here.
     *
     * @param node node to start und calculate the weight for
     *
     * @return
     */
    public static int calculateSubtreeWeight(final ConnectionTreeNode node) {
        return getSubtreeWeight(node, null);
    }

    /**
     * Returns the summed subtree weight starting at a specific node in a connection
     * tree. The weight of starting node and the weight of bond to its given parent
     * node is included here.
     *
     * @param node
     * @param parentNode
     *
     * @return
     */
    private static int getSubtreeWeight(final ConnectionTreeNode node, final ConnectionTreeNode parentNode) {
        int weight = getNodeWeight(node, parentNode);
        for (final ConnectionTreeNode childNode : node.getChildNodes()) {
            weight += getSubtreeWeight(childNode, node);
        }

        return weight;
    }

    /**
     * Returns the weight for a node and its connection to a parent node
     * (optional).
     *
     * @param node       node to get the weight from
     * @param parentNode parent node of node or null
     *
     * @return the priority weight for node; plus the weight of
     * the bond to its parent node if the parent node is not null
     *
     * @see Utils#getSymbolPriorityWeight(String)
     */
    public static Integer getNodeWeight(final ConnectionTreeNode node, final ConnectionTreeNode parentNode) {
        int weight = 0;
        if (parentNode
                != null) {
            final String bondSymbol = getSymbolForBond(node.getBondToParent());
            if (bondSymbol
                    == null) {
                return null;
            }
            // add weight for bond type priority
            if (!bondSymbol.isEmpty()) {
                weight += getSymbolPriorityWeight(bondSymbol);
            }

        }
        // add weight for further symbol priority
        if (node.isRingClosureNode()) {
            weight += getSymbolPriorityWeight("&");
        } else {
            weight += getSymbolPriorityWeight(node.getAtom()
                                                  .getSymbol());
            //            weight -= node.getAtom()
            //                          .getImplicitHydrogenCount()
            //                    * getSymbolPriorityWeight("H");
        }

        return weight;
    }

    /**
     * Returns an ArrayList of ranked child node indices for a tree node.
     *
     * @param node node to rank the children
     *
     * @return
     *
     * @see #getNodeWeight(ConnectionTreeNode, ConnectionTreeNode)
     */
    private static List<Integer> getRankedChildNodesIndices(final ConnectionTreeNode node) {
        final List<ConnectionTreeNode> childNodes = node.getChildNodes();
        final List<Integer> rankedChildNodesIndices = new ArrayList<>();
        for (int i = 0; i
                < childNodes.size(); i++) {
            rankedChildNodesIndices.add(i);
        }
        rankedChildNodesIndices.sort((childNodeIndex1, childNodeIndex2) -> {
            final int nodeWeightsComp = -1
                    * Integer.compare(getNodeWeight(childNodes.get(childNodeIndex1), node),
                                      getNodeWeight(childNodes.get(childNodeIndex2), node));
            if (nodeWeightsComp
                    != 0) {
                return nodeWeightsComp;
            }
            return -1
                    * Integer.compare(calculateSubtreeWeight(childNodes.get(childNodeIndex1)),
                                      calculateSubtreeWeight(childNodes.get(childNodeIndex2)));
        });

        return rankedChildNodesIndices;
    }

    /**
     * Sorts the child nodes of a node by HOSE code priority and weight.
     *
     * @param node node with child nodes to rank
     *
     * @see #getNodeWeight(ConnectionTreeNode, ConnectionTreeNode)
     */
    private static void rankChildNodes(final ConnectionTreeNode node) {
        final List<Integer> rankedChildNodesIndices = getRankedChildNodesIndices(node);
        final List<ConnectionTreeNode> rankedChildNodes = new ArrayList<>();
        final List<IBond> rankedChildNodeBonds = new ArrayList<>();
        for (int i = 0; i
                < rankedChildNodesIndices.size(); i++) {
            rankedChildNodes.add(node.getChildNodes()
                                     .get(rankedChildNodesIndices.get(i)));
            rankedChildNodeBonds.add(node.getBondsToChildren()
                                         .get(rankedChildNodesIndices.get(i)));
        }
        node.getChildNodes()
            .clear();
        node.getBondsToChildren()
            .clear();
        node.getChildNodes()
            .addAll(rankedChildNodes);
        node.getBondsToChildren()
            .addAll(rankedChildNodeBonds);
    }

    /**
     * Sorts the child nodes of each node in the connection tree by HOSE code
     * priority and weight.
     *
     * @param connectionTree connection tree where to rank the child nodes of
     *                       each node.
     *
     * @see #rankChildNodes(ConnectionTreeNode)
     */
    public static void rankChildNodes(final ConnectionTree connectionTree) {
        List<ConnectionTreeNode> nodesInSphere;
        for (int sphere = 0; sphere
                < connectionTree.getMaxSphere(); sphere++) {
            nodesInSphere = connectionTree.getNodesInSphere(sphere, true);
            // for all nodes in sphere
            for (int i = 0; i
                    < nodesInSphere.size(); i++) {
                // findHits all child nodes of that node
                if (nodesInSphere.get(i)
                                 .hasChildren()) {
                    rankChildNodes(nodesInSphere.get(i));
                }
            }
        }
    }

    /**
     * Returns the number of non-empty spheres.
     * For example: C-3;() -> 1, C-3;=N(C/) -> 3
     *
     * @param HOSECode
     *
     * @return
     */
    public static int getSpheresCount(final String HOSECode) {
        int spheresCount = 0;
        for (final String sphere : splitHOSECodeIntoSpheres(HOSECode)) {
            if (!sphere.trim()
                       .isEmpty()) {
                spheresCount++;
            }
        }
        return spheresCount;
    }

    /**
     * Splits a HOSE code into a list of spheres as strings.
     *
     * @param HOSECode HOSE code
     *
     * @return ArrayList of all sphere strings
     */
    public static List<String> splitHOSECodeIntoSpheres(final String HOSECode) {
        final List<String> HOSECodeSpheres = new ArrayList<>();
        final String[] splitSpheres_0_1 = HOSECode.split(";");
        final String[] splitSpheres_1_2 = splitSpheres_0_1[1].split("\\(");
        final String[] splitSpheres_2_n = splitSpheres_1_2[1].substring(0, splitSpheres_1_2[1].length()
                - 1)
                                                             .split("/");
        HOSECodeSpheres.add(splitSpheres_0_1[0]);
        HOSECodeSpheres.add(splitSpheres_1_2[0]);
        for (int s = 0; s
                < splitSpheres_2_n.length; s++) {
            HOSECodeSpheres.add(splitSpheres_2_n[s]);
        }

        return HOSECodeSpheres;
    }

    /**
     * Splits a HOSE code sphere into its positions. Each position includes all
     * its elements.
     * Example: {@code /CC,*N&,C/} results in: {@code {0: [C,C], 1: [*N,&], 2: [C]}}
     *
     * @param HOSECodeSphere HOSE code sphere
     * @param isCenterSphere whether center (zeroth) sphere is given
     *
     * @return HashMap of ArrayLists containing elements for each position of that HOSE code sphere
     */
    public static Map<Integer, ArrayList<String>> splitHOSECodeSphereIntoPositions(final String HOSECodeSphere,
                                                                                   final boolean isCenterSphere) {
        final Map<Integer, ArrayList<String>> positions = new HashMap<>();
        // zeroth sphere
        if (isCenterSphere) {
            positions.put(0, new ArrayList<>());
            // add element
            positions.get(0)
                     .add(HOSECodeSphere.split("-")[0]);
            // zeroth sphere contains charges
            if (HOSECodeSphere.endsWith("-")) { // add negative formal charge with value 1
                positions.get(0)
                         .add("-1");
            } else if (HOSECodeSphere.endsWith("+")) {// add positive formal charge with value 1
                positions.get(0)
                         .add("1");
            } else if (HOSECodeSphere.endsWith("'")) { // add formal charge with a higher value
                positions.get(0)
                         .add(HOSECodeSphere.split("'")[1]);
            } else {
                positions.get(0)
                         .add("0");
            }

            return positions;
        }
        // higher spheres
        char c;
        StringBuilder elem = new StringBuilder();
        int positionCounter = 0;
        boolean formalChargeDetected = false;
        for (int i = 0; i
                < HOSECodeSphere.length(); i++) {
            c = HOSECodeSphere.charAt(i);
            if ((c
                    == '=')
                    || (c
                    == '%')
                    || (c
                    == '*')) {
                if (elem.length()
                        > 0) {
                    if (!positions.containsKey(positionCounter)) {
                        positions.put(positionCounter, new ArrayList<>());
                    }
                    positions.get(positionCounter)
                             .add(elem.toString());
                    elem = new StringBuilder();
                }
                elem.append(c);
            } else if (Character.isUpperCase(c)
                    || (c
                    == '&')) {
                if ((elem.length()
                        > 0)
                        && (Character.isLetter(elem.charAt(elem.length()
                                                                   - 1))
                        || (elem.charAt(elem.length()
                                                - 1)
                        == '&'))) {
                    if (!positions.containsKey(positionCounter)) {
                        positions.put(positionCounter, new ArrayList<>());
                    }
                    positions.get(positionCounter)
                             .add(elem.toString());
                    elem = new StringBuilder();
                    elem.append(c);
                } else if (formalChargeDetected) {
                    if (!positions.containsKey(positionCounter)) {
                        positions.put(positionCounter, new ArrayList<>());
                    }
                    positions.get(positionCounter)
                             .add(elem.toString());
                    elem = new StringBuilder();
                    elem.append(c);
                    formalChargeDetected = false;
                } else {
                    elem.append(c);
                }
            } else if (Character.isLowerCase(c)) {
                elem.append(c);
            } else if ((c
                    == '-')
                    || (c
                    == '+')
                    || Character.isDigit(c)) {
                elem.append(c);
                formalChargeDetected = true;
            } else if (c
                    == ',') {
                if (!positions.containsKey(positionCounter)) {
                    positions.put(positionCounter, new ArrayList<>());
                }
                if (elem.length()
                        == 0) {
                    positions.get(positionCounter)
                             .add(null);
                } else {
                    positions.get(positionCounter)
                             .add(elem.toString());
                    elem = new StringBuilder();
                }
                positionCounter++;
                formalChargeDetected = false;
            }
        }
        // add last element
        if (elem.length()
                > 0) {
            if (!positions.containsKey(positionCounter)) {
                positions.put(positionCounter, new ArrayList<>());
            }
            positions.get(positionCounter)
                     .add(elem.toString());
        } else if (HOSECodeSphere.endsWith(",")) {
            if (!positions.containsKey(positionCounter)) {
                positions.put(positionCounter, new ArrayList<>());
            }
            positions.get(positionCounter)
                     .add(null);
        }

        return positions;
    }

    /**
     * Counts the number of occurring atoms within a given HOSE code.
     *
     * @param HOSECode HOSE code to analyse
     *
     * @return number of atoms within HOSE code
     */
    public static int countAtoms(final String HOSECode) {
        int counter = 0;
        for (int k = 0; k
                < HOSECode.length(); k++) {
            // Check for uppercase letters
            if (Character.isLetter(HOSECode.charAt(k))
                    && Character.isUpperCase(HOSECode.charAt(k))) {
                counter++;
            }
        }

        return counter;
    }

    /**
     * Returns the weight/cost for an HOSE code symbol regarding its priority.
     *
     * @param symbol HOSE code symbol
     *
     * @return weight/cost for the symbol
     */
    public static int getSymbolPriorityWeight(final String symbol) {
        switch (symbol) {
            case "%":
                return 15;
            case "=":
                return 14;
            case "*":
                return 13;
            case "C":
                return 12;
            case "O":
                return 11;
            case "N":
                return 10;
            case "S":
                return 9;
            case "P":
                return 8;
            case "Si":
            case "Q":
                return 7;
            case "B":
                return 6;
            case "F":
                return 5;
            case "Cl":
            case "X":
                return 4;
            case "Br":
            case "Y":
                return 3;
            case "I":
                return 2;
            case "&":
            case "H":
                return 1;
        }

        return 0;
    }

    /**
     * Converts an element symbol into notation as shown in origin article by
     * Bremser.
     * That includes: Si -> Q, Cl -> X, Br -> Y
     *
     * @param element
     *
     * @return HOSE code symbol as in origin article by Bremser
     */
    public static String toHOSECodeSymbol(final String element) {
        if (element.equals("Si")) {
            return "Q";
        }
        if (element.equals("Cl")) {
            return "X";
        }
        if (element.equals("Br")) {
            return "Y";
        }

        return element;
    }

    /**
     * Converts an HOSE code symbol as shown in origin article by
     * Bremser into default element notation.
     * That includes: Q -> Si, X -> Cl, Y -> Br
     *
     * @param symbol
     *
     * @return default element notation
     */
    public static String toElementSymbol(final String symbol) {
        if (symbol.equals("Q")) {
            return "Si";
        }
        if (symbol.equals("X")) {
            return "Cl";
        }
        if (symbol.equals("Y")) {
            return "Br";
        }

        return symbol;
    }

    /**
     * Returns the notation of bond information used in HOSE code.
     * The bond has to contain its bond order and aromaticity information.
     *
     * @param bond bond containing bond order and aromaticity information
     *
     * @return HOSE code symbol for a bond
     */
    public static String getSymbolForBond(final IBond bond) {
        if (bond
                != null) {
            if (bond.isAromatic()) {
                return "*";
            }
            switch (bond.getOrder()) {
                case SINGLE:
                    return "";
                case DOUBLE:
                    return "=";
                case TRIPLE:
                    return "%";
            }
        }

        return null;
    }

    /**
     * Returns the bond order from an HOSE code bond symbol.
     * One has to consider that in this direction the aromatic HOSE code symbol
     * (*) is ambiguous. It means either a single or a double aromatic bond.
     * For that case, a single bond will always be returned.
     *
     * @param symbol HOSE code bond symbol
     *
     * @return bond order for a bond symbol or if symbol is unknown
     */
    public static IBond.Order getBondOrderForSymbol(final String symbol) {

        switch (symbol) {
            case "":
            case "*":
                return IBond.Order.SINGLE;
            case "=":
                return IBond.Order.DOUBLE;
            case "%":
                return IBond.Order.TRIPLE;
        }

        return null;
    }
}
