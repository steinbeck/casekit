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
package casekit.nmr.fragments.model;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;

import java.util.ArrayList;
import java.util.List;

/**
 * Represents a node in a connection tree {@link ConnectionTree}.
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class ConnectionTreeNode {

    private final List<ConnectionTreeNode> children;
    private final List<IBond> bondsToChildren;
    private ConnectionTreeNode parent;
    private IBond bondToParent;
    private IAtom atom;
    private Integer key;
    private int sphere;
    private boolean isRingClosure;
    private ConnectionTreeNode ringClosureParent;
    private boolean isPseudoNode;

    /**
     * Pre-defined constructor for creating a non-ring closure node.
     *
     * @param atom
     * @param key
     * @param sphere
     */
    public ConnectionTreeNode(final IAtom atom, final int key, final int sphere, final ConnectionTreeNode parent,
                              final IBond bondToParent) {
        this.atom = atom;
        this.key = key;
        this.sphere = sphere;
        this.parent = parent;
        this.bondToParent = bondToParent;
        this.children = new ArrayList<>();
        this.bondsToChildren = new ArrayList<>();
        this.isRingClosure = false;
        this.isPseudoNode = false;
    }

    /**
     * Pre-defined constructor for creating a ring closure node.
     *
     * @param sphere
     * @param ringClosurePartner
     */
    public ConnectionTreeNode(final ConnectionTreeNode ringClosurePartner, final int sphere,
                              final ConnectionTreeNode parent, final IBond bondToParent) {
        this.sphere = sphere;
        this.parent = parent;
        this.bondToParent = bondToParent;
        this.children = new ArrayList<>();
        this.bondsToChildren = new ArrayList<>();
        this.isRingClosure = true;
        this.ringClosureParent = ringClosurePartner;
        this.isPseudoNode = false;
    }

    public IAtom getAtom() {
        return this.atom;
    }

    public void setAtom(final IAtom atom) {
        this.atom = atom;
    }

    public ConnectionTreeNode getParent() {
        return this.parent;
    }

    public void setParent(final ConnectionTreeNode parent) {
        this.parent = parent;
    }

    public IBond getBondToParent() {
        return this.bondToParent;
    }

    public void setBondToParent(final IBond bondToParent) {
        this.bondToParent = bondToParent;
    }

    public List<ConnectionTreeNode> getChildNodes() {
        return this.children;
    }

    public List<IBond> getBondsToChildren() {
        return this.bondsToChildren;
    }

    public Integer getKey() {
        return this.key;
    }

    public void setKey(final int key) {
        this.key = key;
    }

    public int getSphere() {
        return this.sphere;
    }

    public void setSphere(final int sphere) {
        this.sphere = sphere;
    }

    public void setIsRingClosureNode(final boolean isRingClosureNode) {
        this.isRingClosure = isRingClosureNode;
    }

    public boolean isRingClosureNode() {
        return this.isRingClosure;
    }

    public ConnectionTreeNode getRingClosureParent() {
        return this.ringClosureParent;
    }

    public void setRingClosureParent(final ConnectionTreeNode ringClosureParent) {
        this.ringClosureParent = ringClosureParent;
    }

    public void setIsPseudoNode(final boolean isPseudoNode) {
        this.isPseudoNode = isPseudoNode;
    }

    public boolean isPseudoNode() {
        return this.isPseudoNode;
    }

    public boolean addChildNode(final ConnectionTreeNode childNode, final IBond bondToChild) {
        return this.addChildNode(childNode, bondToChild, this.getChildNodes()
                                                             .size());
    }

    public boolean addChildNode(final ConnectionTreeNode childNode, final IBond bondToChild, final int pos) {
        if (!this.checkListIndex(this.getChildNodes()
                                     .size(), pos)) {
            return false;
        }
        this.getChildNodes()
            .add(pos, childNode);
        this.getBondsToChildren()
            .add(pos, bondToChild);

        return true;
    }

    private boolean checkListIndex(final int listLength, final int pos) {
        return (pos
                >= 0)
                && (pos
                <= listLength);
    }

    public boolean removeChildNode(final ConnectionTreeNode childNode) {
        final int indexOfChildNode = this.getChildNodes()
                                         .indexOf(childNode);
        if (indexOfChildNode
                == -1) {
            return false;
        }
        this.getChildNodes()
            .remove(indexOfChildNode);
        this.getBondsToChildren()
            .remove(indexOfChildNode);

        return true;
    }

    public boolean hasAParent() {
        return this.parent
                != null;
    }

    public boolean hasChild(final int childKey) {
        for (final ConnectionTreeNode childNode : this.children) {
            if (!childNode.isRingClosureNode()
                    && (childNode.getKey()
                    == childKey)) {
                return true;
            }
        }

        return false;
    }

    public boolean hasChildren() {
        return !this.children.isEmpty();
    }

    @Override
    public String toString() {
        return "ConnectionTreeNode{"
                + "key="
                + this.key
                + ", sphere="
                + this.sphere
                + ", isRingClosure="
                + this.isRingClosure
                + ", isPseudoNode="
                + this.isPseudoNode
                + ", ..."
                //                + ", children="
                //                + this.children
                //                + ", bondsToChildren="
                //                + this.bondsToChildren
                //                + ", parent="
                //                + this.parent
                //                + ", bondToParent="
                //                + this.bondToParent
                //                + ", atom="
                //                + this.atom
                //                + ", ringClosureParent="
                //                + this.ringClosureParent
                + '}';
    }
}
