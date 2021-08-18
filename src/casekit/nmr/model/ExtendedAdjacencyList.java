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
package casekit.nmr.model;

import casekit.nmr.utils.Utils;
import lombok.AllArgsConstructor;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.openscience.cdk.graph.matrix.ConnectionMatrix;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType.Hybridization;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.Atom;
import org.openscience.cdk.silent.Bond;
import org.openscience.cdk.silent.PseudoAtom;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
@NoArgsConstructor
@AllArgsConstructor
@Getter
@Setter
public class ExtendedAdjacencyList {

    private int[][][] adjacencyList; // connected atom index, bond order, bond is in ring, bond is aromatic
    private String[] atomTypes;
    private Integer[][] atomProperties;// hydrogenCounts, valencies, formalCharges, isInRingAtoms, isAromaticAtoms
    private Hybridization[] hybridizations;
    private int bondCount;


    public ExtendedAdjacencyList(final IAtomContainer ac) {
        final double[][] connectionMatrix = ConnectionMatrix.getMatrix(ac);
        this.adjacencyList = new int[connectionMatrix.length][][];
        List<int[]> connectedAtomsList;
        int[][] temp;
        IBond bond;
        for (int i = 0; i
                < connectionMatrix.length; i++) {
            connectedAtomsList = new ArrayList<>();
            for (int j = 0; j
                    < connectionMatrix[i].length; j++) {
                if (connectionMatrix[i][j]
                        >= 1) {
                    bond = ac.getBond(ac.getAtom(i), ac.getAtom(j));
                    connectedAtomsList.add(new int[]{j, (int) connectionMatrix[i][j], bond.isInRing()
                                                                                      ? 1
                                                                                      : 0, bond.isAromatic()
                                                                                           ? 1
                                                                                           : 0});
                }
            }
            temp = new int[connectedAtomsList.size()][];
            for (int k = 0; k
                    < connectedAtomsList.size(); k++) {
                temp[k] = connectedAtomsList.get(k);
            }
            this.adjacencyList[i] = temp;
        }
        this.atomTypes = new String[this.adjacencyList.length];
        this.hybridizations = new Hybridization[this.adjacencyList.length];
        this.atomProperties = new Integer[this.adjacencyList.length][];

        IAtom atom;
        for (int i = 0; i
                < this.adjacencyList.length; i++) {
            atom = ac.getAtom(i);
            this.setAtomProperties(i, atom.getSymbol(), atom.getImplicitHydrogenCount(), atom.getValency(),
                                   atom.getFormalCharge(), atom.isInRing(), atom.isAromatic(), atom.getHybridization());
        }
        this.updateBondCount();
    }

    private void setAtomProperties(final int atomIndex, final String atomType, final Integer implicitHydrogenCount,
                                   final Integer valency, final Integer formalCharge, final Boolean isInRing,
                                   final Boolean isAromatic, final Hybridization hybridization) {
        this.atomTypes[atomIndex] = atomType;
        this.atomProperties[atomIndex] = new Integer[5];
        this.atomProperties[atomIndex][0] = implicitHydrogenCount;
        this.atomProperties[atomIndex][1] = valency;
        this.atomProperties[atomIndex][2] = formalCharge;
        this.atomProperties[atomIndex][3] = isInRing
                                            ? 1
                                            : 0;
        this.atomProperties[atomIndex][4] = isAromatic
                                            ? 1
                                            : 0;
        this.hybridizations[atomIndex] = hybridization;
    }

    private void updateBondCount() {
        int bondCounter = 0;
        for (int i = 0; i
                < this.adjacencyList.length; i++) {
            bondCounter += this.adjacencyList[i].length;
        }
        this.bondCount = bondCounter
                / 2;
    }

    public int getAtomCount() {
        return this.adjacencyList.length;
    }


    public IAtomContainer toAtomContainer() {
        final IAtomContainer ac = SilentChemObjectBuilder.getInstance()
                                                         .newAtomContainer();
        IAtom atom;
        for (int i = 0; i
                < this.adjacencyList.length; i++) {
            if (this.atomTypes[i].equals("R")) {
                atom = new PseudoAtom("R");
            } else {
                atom = new Atom(this.atomTypes[i]);
            }
            atom.setImplicitHydrogenCount(this.atomProperties[i][0]);
            atom.setValency(this.atomProperties[i][1]);
            atom.setFormalCharge(this.atomProperties[i][2]);
            atom.setIsInRing(this.atomProperties[i][3]
                                     == 1);
            atom.setIsAromatic(this.atomProperties[i][4]
                                       == 1);
            atom.setHybridization(this.hybridizations[i]);

            ac.addAtom(atom);
        }
        IBond bond;
        for (int i = 0; i
                < this.adjacencyList.length; i++) {
            for (int k = 0; k
                    < this.adjacencyList[i].length; k++) {
                if (ac.getBond(ac.getAtom(i), ac.getAtom(this.adjacencyList[i][k][0]))
                        == null) {
                    bond = new Bond(ac.getAtom(i), ac.getAtom(this.adjacencyList[i][k][0]),
                                    Utils.getBondOrder(this.adjacencyList[i][k][1]));
                    bond.setIsInRing(this.adjacencyList[i][k][2]
                                             == 1);
                    bond.setIsAromatic(this.adjacencyList[i][k][3]
                                               == 1);
                    ac.addBond(bond);
                }
            }
        }

        return ac;
    }

    public ExtendedAdjacencyList buildClone() {
        return new ExtendedAdjacencyList(this.toAtomContainer());
    }

    @Override
    public String toString() {
        return "ExtendedAdjacencyList{"
                + "atomCount="
                + this.getAtomCount()
                + ", bondCount="
                + this.bondCount
                + ", adjacencyList="
                + Arrays.deepToString(this.adjacencyList)
                + ", atomTypes="
                + Arrays.toString(this.atomTypes)
                + ", atomProperties="
                + Arrays.deepToString(this.atomProperties)
                + ", hybridizations="
                + Arrays.toString(this.hybridizations)
                + '}';
    }
}
