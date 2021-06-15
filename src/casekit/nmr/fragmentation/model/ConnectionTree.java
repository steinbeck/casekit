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
package casekit.nmr.fragmentation.model;

import casekit.nmr.hose.Utils;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;

import java.util.*;

/**
 * Represents a tree of connected atoms (nodes) of a molecule
 * in a parent-child-relationship.
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 * @see ConnectionTreeNode
 */
public class ConnectionTree {
    private final ConnectionTreeNode root;
    private final Set<Integer> keySet;
    private int maxSphere;

    public ConnectionTree(final IAtom rootAtom, final int key) {
        this.root = new ConnectionTreeNode(rootAtom, key, 0, null, null);
        this.keySet = new HashSet<>();
        this.keySet.add(this.root.getKey());
        this.maxSphere = 0;
    }

    /**
     * Checks whether two nodes form a ring closures in one connection tree.
     *
     * @param node1 first node
     * @param node2 second node
     *
     * @return
     */
    public static boolean nodesFormRingClosure(final ConnectionTreeNode node1, final ConnectionTreeNode node2) {
        return ConnectionTree.hasRingClosureParent(node1, node2)
                && ConnectionTree.hasRingClosureParent(node2, node1);
    }

    public static boolean hasRingClosureParent(final ConnectionTreeNode node,
                                               final ConnectionTreeNode ringClosureParentNode) {
        for (final ConnectionTreeNode childNode : node.getChildNodes()) {
            if (childNode.isRingClosureNode()
                    && (childNode.getRingClosureParent()
                                 .getKey()
                    == ringClosureParentNode.getKey())) {
                return true;
            }
        }

        return false;
    }

    /**
     * Checks whether a node is at a ring closures.
     *
     * @param node node to check
     *
     * @return
     */
    public static boolean isAtRingClosure(final ConnectionTreeNode node) {
        for (final ConnectionTreeNode childNode : node.getChildNodes()) {
            if (childNode.isRingClosureNode()) {
                return true;
            }
        }
        return false;
    }

    /**
     * Returns a subtree of a given connection tree by using a given subtree root node.
     *
     * @param connectionTree connection tree
     * @param rootNodeKey    root node key for subtree to create
     *
     * @return
     */
    public static ConnectionTree buildSubtree(final ConnectionTree connectionTree, final int rootNodeKey) {
        if (!connectionTree.containsKey(rootNodeKey)) {
            return null;
        }

        final ConnectionTreeNode rootNode = connectionTree.getNode(rootNodeKey);
        final ConnectionTree subtree = new ConnectionTree(rootNode.getAtom(), rootNode.getKey());

        buildSubtree(subtree, rootNode, 1);

        return subtree;
    }

    private static void buildSubtree(final ConnectionTree subtree, final ConnectionTreeNode parentNode,
                                     final int sphere) {
        int childNodeIndex = 0;
        for (final ConnectionTreeNode childNode : parentNode.getChildNodes()) {
            if (childNode.isRingClosureNode()) {
                subtree.addRingClosureNode(parentNode.getKey(), childNode.getRingClosureParent()
                                                                         .getKey(), childNode.getBondToParent());
            } else {
                subtree.addNode(childNode.getAtom(), childNode.getKey(), parentNode.getKey(),
                                parentNode.getBondsToChildren()
                                          .get(childNodeIndex));
            }
            buildSubtree(subtree, childNode, sphere
                    + 1);
            childNodeIndex++;
        }
    }

    public ConnectionTreeNode getRootNode() {
        return this.root;
    }

    public boolean addNode(final IAtom newNodeAtomData, final int newNodeKey, final int parentNodeKey,
                           final IBond bondToParent) {
        if (this.containsKey(newNodeKey)) {
            return false;
        }
        final ConnectionTreeNode parentNode = this.getNode(parentNodeKey);
        this.addNode(new ConnectionTreeNode(newNodeAtomData, newNodeKey, parentNode.getSphere()
                + 1, parentNode, bondToParent), parentNode);

        return true;
    }

    public boolean addRingClosureNode(final int parentNodeKey, final int ringClosureParentNodeKey,
                                      final IBond bondToParent) {
        if (!this.containsKey(parentNodeKey)
                || !this.containsKey(ringClosureParentNodeKey)) {
            return false;
        }
        final ConnectionTreeNode parentNode = this.getNode(parentNodeKey);
        this.addNode(new ConnectionTreeNode(this.getNode(ringClosureParentNodeKey), parentNode.getSphere()
                + 1, parentNode, bondToParent), parentNode);

        return true;
    }

    private void addNode(final ConnectionTreeNode newNode, final ConnectionTreeNode parentNode) {
        parentNode.addChildNode(newNode, newNode.getBondToParent());

        if (!newNode.isRingClosureNode()) {
            this.keySet.add(newNode.getKey());
        }
        if (newNode.getSphere()
                > this.maxSphere) {
            this.maxSphere = newNode.getSphere();
        }
    }

    public int getMaxSphere(final boolean withRingClosureNodes) {
        if (!withRingClosureNodes
                && this.getNodesInSphere(this.maxSphere, false)
                       .isEmpty()) {
            return this.maxSphere
                    - 1;
        }
        return this.maxSphere;
    }

    public int getNodesCount(final boolean withRingClosureNodes) {
        return this.getNodes(withRingClosureNodes)
                   .size();
    }

    public int getNodesCountInSphere(final int sphere, final boolean withRingClosureNodes) {
        return this.getNodesInSphere(sphere, withRingClosureNodes)
                   .size();
    }

    public List<Integer> getKeys() {
        final List<Integer> keys = new ArrayList<>();
        for (int s = 0; s
                <= this.getMaxSphere(false); s++) {
            for (final ConnectionTreeNode nodeInSphere : this.getNodesInSphere(s, false)) {
                keys.add(nodeInSphere.getKey());
            }
        }

        return keys;
    }

    public List<ConnectionTreeNode> getNodes(final boolean withRingClosureNodes) {
        final List<ConnectionTreeNode> nodes = new ArrayList<>();
        for (int s = 0; s
                <= this.getMaxSphere(withRingClosureNodes); s++) {
            nodes.addAll(this.getNodesInSphere(s, withRingClosureNodes));
        }

        return nodes;
    }

    public boolean containsKey(final int key) {
        return this.keySet.contains(key);
    }

    public ConnectionTreeNode getNode(final int key) {
        if (!this.containsKey(key)) {
            return null;
        }

        return this.findNode(key, this.root);
    }

    private ConnectionTreeNode findNode(final int key, final ConnectionTreeNode currentNode) {
        if (currentNode.isRingClosureNode()) {
            return null;
        }
        if (currentNode.getKey()
                == key) {
            return currentNode;
        }
        ConnectionTreeNode result = null;
        for (final ConnectionTreeNode childNode : currentNode.getChildNodes()) {
            result = this.findNode(key, childNode);
            if ((result
                    != null)
                    && (result.getKey()
                    == key)) {
                break;
            }
        }

        return result;
    }

    public int getNodeIndexInSphere(final ConnectionTreeNode node, final int sphere) {

        return this.getNodesInSphere(sphere, true)
                   .indexOf(node);
    }

    public List<Integer> getNodeKeysInSphere(final int sphere) {
        final List<Integer> keys = new ArrayList<>();
        for (final ConnectionTreeNode treeNode : this.getNodesInSphere(sphere, false)) {
            if (!treeNode.isRingClosureNode()) {
                keys.add(treeNode.getKey());
            }
        }

        return keys;
    }

    public List<ConnectionTreeNode> getNodesInSphere(final int sphere, final boolean withRingClosureNodes) {
        final List<ConnectionTreeNode> nodesInSphere = this.findNodesInSphere(sphere, this.root, new ArrayList<>());
        if (withRingClosureNodes) {
            return nodesInSphere;
        }
        // remove ring closure nodes
        final List<ConnectionTreeNode> nodesInSphereToRemove = new ArrayList<>();
        for (final ConnectionTreeNode nodeInSphere : nodesInSphere) {
            if (nodeInSphere.isRingClosureNode()) {
                nodesInSphereToRemove.add(nodeInSphere);
            }
        }
        nodesInSphere.removeAll(nodesInSphereToRemove);

        return nodesInSphere;
    }

    private List<ConnectionTreeNode> findNodesInSphere(final int sphere, final ConnectionTreeNode currentNode,
                                                       final List<ConnectionTreeNode> indicesInSphere) {
        if (currentNode.getSphere()
                == sphere) {
            indicesInSphere.add(currentNode);
            return indicesInSphere;
        }
        for (final ConnectionTreeNode childNode : currentNode.getChildNodes()) {
            this.findNodesInSphere(sphere, childNode, indicesInSphere);
        }

        return indicesInSphere;
    }

    public IBond getBond(final int nodeKey1, final int nodeKey2) {
        if (!this.containsKey(nodeKey1)
                || !this.containsKey(nodeKey2)) {
            return null;
        }
        // node1 and node2 have parent-child-relationship
        final ConnectionTreeNode node1 = this.getNode(nodeKey1);
        for (final ConnectionTreeNode childNode : node1.getChildNodes()) {
            if (!childNode.isRingClosureNode()
                    && (childNode.getKey()
                    == nodeKey2)) {
                return childNode.getBondToParent();
            }
        }
        // if nodes form a ring closure
        if (ConnectionTree.nodesFormRingClosure(node1, this.getNode(nodeKey2))) {
            for (final ConnectionTreeNode childNode : node1.getChildNodes()) {
                if (childNode.isRingClosureNode()
                        && (childNode.getRingClosureParent()
                                     .getKey()
                        == nodeKey2)) {
                    return childNode.getBondToParent();
                }
            }
        }

        return null;
    }

    public boolean hasParent(final int key, final int parentKey) {
        if (!this.containsKey(key)
                || !this.containsKey(parentKey)) {
            return false;
        }

        return this.getNode(key)
                   .getParent()
                   .getKey()
                == parentKey;
    }

    public boolean hasChild(final int key, final int childKey) {
        if (!this.containsKey(key)
                || !this.containsKey(childKey)) {
            return false;
        }

        return this.getNode(key)
                   .hasChild(childKey);
    }

    public void removeNode(final int key) {
        final ConnectionTreeNode node = this.getNode(key);
        if (node
                != null) {
            final List<ConnectionTreeNode> children = new ArrayList<>(node.getChildNodes());
            for (final ConnectionTreeNode childNode : children) {
                if (childNode.isRingClosureNode()) {
                    childNode.getRingClosureParent()
                             .removeChildNode(childNode);
                } else {
                    this.removeNode(childNode.getKey());
                }
                node.removeChildNode(childNode);
                this.keySet.remove(childNode.getKey());
            }
            final ConnectionTreeNode parent = node.getParent();
            if (parent
                    != null) {
                parent.removeChildNode(node);
            }
            this.keySet.remove(node.getKey());
        }
    }

    /**
     * @param parentKey
     * @param childKey1
     * @param childKey2
     *
     * @return
     *
     * @deprecated
     */
    public boolean swapChildNodes(final int parentKey, final int childKey1, final int childKey2) {
        if (!this.containsKey(parentKey)
                || !this.containsKey(childKey1)
                || !this.containsKey(childKey2)
                || !this.hasChild(parentKey, childKey1)
                || !this.hasChild(parentKey, childKey2)
                || !this.hasParent(childKey1, parentKey)
                || !this.hasParent(childKey2, parentKey)) {
            return false;
        }
        final ConnectionTreeNode parentNode = this.getNode(parentKey);
        final ConnectionTreeNode childNode1 = this.getNode(childKey1);
        final ConnectionTreeNode childNode2 = this.getNode(childKey2);

        final int indexChildNode1 = parentNode.getChildNodes()
                                              .indexOf(childNode1);
        final int indexChildNode2 = parentNode.getChildNodes()
                                              .indexOf(childNode2);
        Collections.swap(parentNode.getChildNodes(), indexChildNode1, indexChildNode2);
        Collections.swap(parentNode.getBondsToChildren(), indexChildNode1, indexChildNode2);


        return (parentNode.getChildNodes()
                          .indexOf(childNode1)
                == indexChildNode2)
                && (parentNode.getChildNodes()
                              .indexOf(childNode2)
                == indexChildNode1);
    }

    @Override
    public String toString() {
        final StringBuilder treeStringBuilder = new StringBuilder();
        for (int s = 0; s
                <= this.maxSphere; s++) {
            treeStringBuilder.append("[")
                             .append(s)
                             .append("]");
            for (final ConnectionTreeNode nodeInSphere : this.getNodesInSphere(s, true)) {
                treeStringBuilder.append(" ");
                if (nodeInSphere.isRingClosureNode()) {
                    treeStringBuilder.append("-")
                                     .append(nodeInSphere.getParent()
                                                         .getKey())
                                     .append(": ");
                } else {
                    treeStringBuilder.append(nodeInSphere.getKey())
                                     .append(": ");
                }
                if (nodeInSphere.hasAParent()) {
                    treeStringBuilder.append(Utils.getSymbolForBond(nodeInSphere.getBondToParent()));
                }
                if (nodeInSphere.isRingClosureNode()) {
                    treeStringBuilder.append("&");
                } else {
                    treeStringBuilder.append(nodeInSphere.getAtom()
                                                         .getSymbol());
                }
                treeStringBuilder.append(" {");
                if (nodeInSphere.isRingClosureNode()) {
                    treeStringBuilder.append(nodeInSphere.getRingClosureParent()
                                                         .getKey());
                } else {
                    if (nodeInSphere.getChildNodes()
                                    .size()
                            > 1) {
                        for (final ConnectionTreeNode childNode : nodeInSphere.getChildNodes()) {
                            if (childNode.isRingClosureNode()) {
                                treeStringBuilder.append("-")
                                                 .append(childNode.getRingClosureParent()
                                                                  .getKey())
                                                 .append(" ");
                            } else {
                                treeStringBuilder.append(childNode.getKey())
                                                 .append(" ");
                            }
                        }
                    } else if (nodeInSphere.getChildNodes()
                                           .size()
                            == 1) {
                        if (nodeInSphere.getChildNodes()
                                        .get(0)
                                        .isRingClosureNode()) {
                            treeStringBuilder.append("-")
                                             .append(nodeInSphere.getChildNodes()
                                                                 .get(0)
                                                                 .getRingClosureParent()
                                                                 .getKey());
                        } else {
                            treeStringBuilder.append(nodeInSphere.getChildNodes()
                                                                 .get(0)
                                                                 .getKey());
                        }

                    }
                }

                treeStringBuilder.append("} ");
            }
        }

        return treeStringBuilder.toString();
    }
}
