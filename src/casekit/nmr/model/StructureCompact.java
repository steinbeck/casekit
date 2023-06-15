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
import org.openscience.cdk.interfaces.IAtomType;
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
public class StructureCompact {

    private int[][][] bondProperties; // connected atom index, bond order, bond is in ring, bond is aromatic, bond stereo configuration
    private Integer[][] atomProperties; // element symbol, hybridization, implicitHydrogenCount, valency, formalCharge, isInRingAtom, isAromaticAtom

    public StructureCompact(final IAtomContainer ac) {
        final double[][] connectionMatrix = ConnectionMatrix.getMatrix(ac);
        this.bondProperties = new int[connectionMatrix.length][][];
        List<int[]> connectedAtomsList;
        int[][] temp;
        IBond bond;
        for (int i = 0; i
                < connectionMatrix.length; i++) {
            connectedAtomsList = new ArrayList<>();
            for (int j = i
                    + 1; j
                         < connectionMatrix[i].length; j++) {
                if (connectionMatrix[i][j]
                        >= 1) {
                    bond = ac.getBond(ac.getAtom(i), ac.getAtom(j));
                    connectedAtomsList.add(new int[]{j, (int) connectionMatrix[i][j], bond.isInRing()
                                                                                      ? 1
                                                                                      : 0, bond.isAromatic()
                                                                                           ? 1
                                                                                           : 0,
                                                     bond.getStereo().ordinal()});
                }
            }
            temp = new int[connectedAtomsList.size()][];
            for (int k = 0; k
                    < connectedAtomsList.size(); k++) {
                temp[k] = connectedAtomsList.get(k);
            }
            this.bondProperties[i] = temp;
        }
        this.atomProperties = new Integer[connectionMatrix.length][];

        IAtom atom;
        for (int i = 0; i
                < connectionMatrix.length; i++) {
            atom = ac.getAtom(i);
            this.atomProperties[i] = new Integer[7];
            this.atomProperties[i][0] = atom.getSymbol()
                                            .equals("R")
                                        ? -1
                                        : atom.getAtomicNumber();
            this.atomProperties[i][1] = atom.getHybridization()
                                                == null
                                        ? -1
                                        : atom.getHybridization()
                                              .ordinal();
            this.atomProperties[i][2] = atom.getImplicitHydrogenCount();
            this.atomProperties[i][3] = atom.getValency();
            this.atomProperties[i][4] = atom.getFormalCharge();
            this.atomProperties[i][5] = atom.isInRing()
                                        ? 1
                                        : 0;
            this.atomProperties[i][6] = atom.isAromatic()
                                        ? 1
                                        : 0;
        }
    }

    public int atomCount() {
        return this.atomProperties.length;
    }

    public int bondCount() {
        int bondCount = 0;
        for (int i = 0; i
                < this.bondProperties.length; i++) {
            bondCount += this.bondProperties[i].length;
        }
        return bondCount;
    }

    public IAtomContainer toAtomContainer() {
        final IAtomContainer ac = SilentChemObjectBuilder.getInstance()
                                                         .newAtomContainer();
        IAtom atom;
        for (int i = 0; i
                < this.atomProperties.length; i++) {
            atom = this.atomProperties[i][0]
                           == -1
                   ? new PseudoAtom("R")
                   : new Atom(this.atomProperties[i][0]);
            atom.setHybridization(this.atomProperties[i][1]
                                          == -1
                                  ? null
                                  : IAtomType.Hybridization.values()[this.atomProperties[i][1]]);
            atom.setImplicitHydrogenCount(this.atomProperties[i][2]);
            atom.setValency(this.atomProperties[i][3]);
            atom.setFormalCharge(this.atomProperties[i][4]);
            atom.setIsInRing(this.atomProperties[i][5]
                                     == 1);
            atom.setIsAromatic(this.atomProperties[i][6]
                                       == 1);

            ac.addAtom(atom);
        }
        IBond bond;
        for (int i = 0; i
                < this.bondProperties.length; i++) {
            for (int k = 0; k
                    < this.bondProperties[i].length; k++) {
                bond = new Bond(ac.getAtom(i), ac.getAtom(this.bondProperties[i][k][0]),
                                Utils.getBondOrder(this.bondProperties[i][k][1]));
                bond.setIsInRing(this.bondProperties[i][k][2]
                                         == 1);
                bond.setIsAromatic(this.bondProperties[i][k][3]
                                           == 1);
                if (this.bondProperties[i][k].length
                        == 5) {
                    // with stereo information
                    bond.setStereo(IBond.Stereo.values()[this.bondProperties[i][k][4]]);
                } else {
                    // without stereo information
                    bond.setStereo(IBond.Stereo.NONE);
                }
                ac.addBond(bond);
            }
        }

        return ac;
    }

    public StructureCompact buildClone() {
        return new StructureCompact(Arrays.copyOf(this.bondProperties, this.bondProperties.length),
                                    Arrays.copyOf(this.atomProperties, this.atomProperties.length));

    }

    @Override
    public String toString() {
        return "StructureCompact{"
                + "bondProperties="
                + Arrays.deepToString(this.bondProperties)
                + ", atomProperties="
                + Arrays.deepToString(this.atomProperties)
                + '}';
    }
}
