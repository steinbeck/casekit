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

import casekit.nmr.Utils;
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
import org.openscience.cdk.silent.SilentChemObjectBuilder;

import java.util.Arrays;

/**
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
@NoArgsConstructor
@AllArgsConstructor
@Getter
@Setter
public class ExtendedConnectionMatrix {

    private double[][] connectionMatrix;
    private String[] atomTypes;
    private Integer[][] atomPropertiesNumeric;// hydrogenCounts, valencies, formalCharges;
    private Hybridization[] hybridizations;
    private Boolean[][] atomPropertiesBoolean;// isInRingAtoms, isAromaticAtoms;
    private Boolean[][][] bondProperties;
    private int bondCount;


    public ExtendedConnectionMatrix(final IAtomContainer ac) {
        this.connectionMatrix = ConnectionMatrix.getMatrix(ac);
        this.atomTypes = new String[this.connectionMatrix.length];
        this.hybridizations = new Hybridization[this.connectionMatrix.length];
        this.atomPropertiesNumeric = new Integer[this.connectionMatrix.length][];
        this.atomPropertiesBoolean = new Boolean[this.connectionMatrix.length][];
        this.bondProperties = new Boolean[this.connectionMatrix.length][][];

        this.init(ac);
    }

    private void init(final IAtomContainer ac) {
        IAtom atom1, atom2;
        IBond bond;
        for (int i = 0; i
                < this.connectionMatrix.length; i++) {
            atom1 = ac.getAtom(i);
            this.setAtomProperties(i, atom1.getSymbol(), atom1.getImplicitHydrogenCount(), atom1.getValency(),
                                   atom1.getFormalCharge(), atom1.isInRing(), atom1.isAromatic(),
                                   atom1.getHybridization());

            this.bondProperties[i] = new Boolean[this.connectionMatrix.length][2];
            for (int k = 0; k
                    < this.connectionMatrix.length; k++) {
                atom2 = ac.getAtom(k);
                bond = ac.getBond(atom1, atom2);
                if (bond
                        != null) {
                    this.setBondProperty(i, k, bond.isInRing(), bond.isAromatic());
                }
            }
        }
        this.updateBondCount();
    }

    private void init(final ExtendedConnectionMatrix extendedConnectionMatrix) {
        for (int i = 0; i
                < this.getAtomCount(); i++) {
            if (i
                    < extendedConnectionMatrix.getAtomCount()) {
                this.setAtomProperties(i, extendedConnectionMatrix.getAtomType(i),
                                       extendedConnectionMatrix.getHydrogenCount(i),
                                       extendedConnectionMatrix.getValency(i),
                                       extendedConnectionMatrix.getFormalCharge(i),
                                       extendedConnectionMatrix.isInRing(i), extendedConnectionMatrix.isAromatic(i),
                                       extendedConnectionMatrix.getHybridization(i));


            }
            this.bondProperties[i] = new Boolean[this.getAtomCount()][2];
            if (i
                    < extendedConnectionMatrix.getAtomCount()) {
                for (int k = 0; k
                        < extendedConnectionMatrix.getAtomCount(); k++) {
                    this.connectionMatrix[i][k] = extendedConnectionMatrix.getBondOrder(i, k);
                    this.setBondProperty(i, k, extendedConnectionMatrix.isInRing(i, k),
                                         extendedConnectionMatrix.isAromatic(i, k));

                }
            } else {
                for (int k = 0; k
                        < this.getAtomCount(); k++) {
                    this.connectionMatrix[i][k] = 0.0;
                    //                    this.setBondProperty(i, k, null, null);
                }
            }
        }
        this.updateBondCount();
    }

    private void extendConnectionMatrix() {
        this.extendConnectionMatrix(1);
    }

    private void extendConnectionMatrix(final int extensionSize) {
        this.connectionMatrix = new double[this.getAtomCount()
                + extensionSize][this.getAtomCount()
                + extensionSize];
        this.atomTypes = new String[this.connectionMatrix.length];
        this.hybridizations = new Hybridization[this.connectionMatrix.length];
        this.atomPropertiesNumeric = new Integer[this.connectionMatrix.length][];
        this.atomPropertiesBoolean = new Boolean[this.connectionMatrix.length][];
        this.bondProperties = new Boolean[this.connectionMatrix.length][][];
    }

    public void addAtom(final String atomType, final Integer implicitHydrogenCount, final Integer valency,
                        final Integer formalCharge, final Boolean isInRing, final Boolean isAromatic,
                        final Hybridization hybridization) {
        // create backup object
        final ExtendedConnectionMatrix extendedConnectionMatrixBackup = this.buildClone();
        // extend the sizes of all matrices by one
        this.extendConnectionMatrix();
        // fill all information in again from backup object
        this.init(extendedConnectionMatrixBackup);
        // set information for new atom
        this.setAtomProperties(this.getAtomCount()
                                       - 1, atomType, implicitHydrogenCount, valency, formalCharge, isInRing,
                               isAromatic, hybridization);
    }

    public boolean addBond(final int atomIndex1, final int atomIndex2, final double order, final Boolean isInRing,
                           final Boolean isAromatic) {
        if (!this.hasAtom(atomIndex1)
                || !this.hasAtom(atomIndex2)) {
            return false;
        }
        if (!this.isValidBondAddition(atomIndex1, atomIndex2, order, isAromatic)) {
            return false;
        }
        this.connectionMatrix[atomIndex1][atomIndex2] = order;
        this.connectionMatrix[atomIndex2][atomIndex1] = order;
        this.setBondProperty(atomIndex1, atomIndex2, isInRing, isAromatic);
        this.setBondProperty(atomIndex2, atomIndex1, isInRing, isAromatic);

        this.updateBondCount();

        return true;
    }

    public boolean isValidBondAddition(final int atomIndex1, final int atomIndex2, final double order,
                                       final boolean isAromatic) {
        if (!this.hasAtom(atomIndex1)
                || !this.hasAtom(atomIndex2)) {
            return false;
        }

        return this.isValidBondAddition(atomIndex1, order, isAromatic)
                && this.isValidBondAddition(atomIndex2, order, isAromatic);
    }

    public boolean isValidBondAddition(final int atomIndex, final double order, final boolean isAromatic) {
        float bondOrderSum = this.getBondOrderSum(atomIndex, true);
        if (isAromatic) {
            bondOrderSum += 1.5;
        } else {
            bondOrderSum += order;
        }
        // -1 for cases with heterocyclic aromatics, like the N in the small aromatic ring in coffein if we want to add the bond to the CH3 group
        if (this.isAromatic(atomIndex)
                && (!this.getAtomType(atomIndex)
                         .equals("C"))) {
            bondOrderSum -= 1;
        }

        return bondOrderSum
                <= this.getValency(atomIndex);
    }

    private void setAtomProperties(final int atomIndex, final String atomType, final Integer implicitHydrogenCount,
                                   final Integer valency, final Integer formalCharge, final Boolean isInRing,
                                   final Boolean isAromatic, final Hybridization hybridization) {
        this.atomTypes[atomIndex] = atomType;
        this.atomPropertiesNumeric[atomIndex] = new Integer[3];
        this.atomPropertiesNumeric[atomIndex][0] = implicitHydrogenCount;
        this.atomPropertiesNumeric[atomIndex][1] = valency;
        this.atomPropertiesNumeric[atomIndex][2] = formalCharge;
        this.atomPropertiesBoolean[atomIndex] = new Boolean[2];
        this.atomPropertiesBoolean[atomIndex][0] = isInRing;
        this.atomPropertiesBoolean[atomIndex][1] = isAromatic;
        this.hybridizations[atomIndex] = hybridization;
    }

    private void setBondProperty(final int atomIndex1, final int atomIndex2, final Boolean isInRing,
                                 final Boolean isAromatic) {
        this.bondProperties[atomIndex1][atomIndex2][0] = isInRing;
        this.bondProperties[atomIndex1][atomIndex2][1] = isAromatic;
    }

    private void updateBondCount() {
        int bondCounter = 0;
        for (int i = 0; i
                < this.getAtomCount(); i++) {
            for (int j = i
                    + 1; j
                         < this.getAtomCount(); j++) {
                if (this.connectionMatrix[i][j]
                        > 0.0) {
                    bondCounter++;
                }
            }
        }
        this.bondCount = bondCounter;
    }

    public Boolean hasBond(final int atomIndex1, final int atomIndex2) {
        if (!this.hasAtom(atomIndex1)
                || !this.hasAtom(atomIndex2)) {
            return null;
        }

        return this.getBondOrder(atomIndex1, atomIndex2)
                > 0.0;
    }

    public Double getBondOrder(final int atomIndex1, final int atomIndex2) {
        if (!this.hasAtom(atomIndex1)
                || !this.hasAtom(atomIndex2)) {
            return null;
        }

        return this.connectionMatrix[atomIndex1][atomIndex2];
    }

    public Float getBondOrderSum(final int atomIndex, final boolean includeHydrogens) {
        if (!this.hasAtom(atomIndex)) {
            return null;
        }
        float bondOrderSum = (float) 0.0;
        for (int j = 0; j
                < this.connectionMatrix[atomIndex].length; j++) {
            if ((this.isAromatic(atomIndex, j)
                    != null)
                    && this.isAromatic(atomIndex, j)) {
                bondOrderSum += 1.5;
            } else {
                bondOrderSum += this.getBondOrder(atomIndex, j);
            }
        }
        if (includeHydrogens) {
            bondOrderSum += this.getHydrogenCount(atomIndex);
        }

        return bondOrderSum;
    }

    public String getAtomType(final int atomIndex) {
        if (!this.hasAtom(atomIndex)) {
            return null;
        }

        return this.atomTypes[atomIndex];
    }

    public Integer getHydrogenCount(final int atomIndex) {
        if (!this.hasAtom(atomIndex)) {
            return null;
        }

        return this.atomPropertiesNumeric[atomIndex][0];
    }

    public Integer getValency(final int atomIndex) {
        if (!this.hasAtom(atomIndex)) {
            return null;
        }

        return this.atomPropertiesNumeric[atomIndex][1];
    }

    public Integer getFormalCharge(final int atomIndex) {
        if (!this.hasAtom(atomIndex)) {
            return null;
        }

        return this.atomPropertiesNumeric[atomIndex][2];
    }

    public Boolean isInRing(final int atomIndex) {
        if (!this.hasAtom(atomIndex)) {
            return null;
        }

        return this.atomPropertiesBoolean[atomIndex][0];
    }

    public Boolean isAromatic(final int atomIndex) {
        if (!this.hasAtom(atomIndex)) {
            return null;
        }

        return this.atomPropertiesBoolean[atomIndex][1];
    }

    public Hybridization getHybridization(final int atomIndex) {
        if (!this.hasAtom(atomIndex)) {
            return null;
        }

        return this.hybridizations[atomIndex];
    }

    public Boolean isInRing(final int atomIndex1, final int atomIndex2) {
        if (!this.hasAtom(atomIndex1)
                || !this.hasAtom(atomIndex2)) {
            return null;
        }

        return this.bondProperties[atomIndex1][atomIndex2][0];
    }

    public Boolean isAromatic(final int atomIndex1, final int atomIndex2) {
        if (!this.hasAtom(atomIndex1)
                || !this.hasAtom(atomIndex2)) {
            return null;
        }

        return this.bondProperties[atomIndex1][atomIndex2][1];
    }

    public int getAtomCount() {
        return this.connectionMatrix.length;
    }

    public Boolean isUnsaturated(final int atomIndex) {
        if (!this.hasAtom(atomIndex)) {
            return null;
        }

        return this.getBondOrderSum(atomIndex, true)
                < this.getValency(atomIndex);
    }

    public boolean hasAtom(final int atomIndex) {
        return (atomIndex
                >= 0)
                && (atomIndex
                < this.getAtomCount());
    }

    public IAtomContainer toAtomContainer() {
        final IAtomContainer ac = SilentChemObjectBuilder.getInstance()
                                                         .newAtomContainer();
        IAtom atom;
        for (int i = 0; i
                < this.connectionMatrix.length; i++) {
            atom = new Atom(this.atomTypes[i]);
            atom.setImplicitHydrogenCount(this.atomPropertiesNumeric[i][0]);
            atom.setValency(this.atomPropertiesNumeric[i][1]);
            atom.setFormalCharge(this.atomPropertiesNumeric[i][2]);
            atom.setHybridization(this.hybridizations[i]);
            atom.setIsInRing(this.atomPropertiesBoolean[i][0]);
            atom.setIsAromatic(this.atomPropertiesBoolean[i][1]);

            ac.addAtom(atom);
        }
        IBond bond;
        for (int i = 0; i
                < this.bondProperties.length; i++) {
            for (int k = i
                    + 1; k
                         < this.bondProperties.length; k++) {
                if (this.connectionMatrix[i][k]
                        > 0.0) {
                    bond = new Bond(ac.getAtom(i), ac.getAtom(k),
                                    Utils.getBondOrder((int) this.connectionMatrix[i][k]));
                    bond.setIsInRing(this.bondProperties[i][k][0]);
                    bond.setIsAromatic(this.bondProperties[i][k][1]);
                    ac.addBond(bond);
                }
            }
        }

        return ac;
    }

    public ExtendedConnectionMatrix buildClone() {
        return new ExtendedConnectionMatrix(this.toAtomContainer());
    }

    @Override
    public String toString() {
        return "ExtendedConnectionMatrix{"
                + "connectionMatrix="
                + Arrays.toString(this.connectionMatrix)
                + ", atomTypes="
                + Arrays.toString(this.atomTypes)
                + ", atomPropertiesNumeric="
                + Arrays.toString(this.atomPropertiesNumeric)
                + ", hybridizations="
                + Arrays.toString(this.hybridizations)
                + ", atomPropertiesBoolean="
                + Arrays.toString(this.atomPropertiesBoolean)
                + ", bondProperties="
                + Arrays.toString(this.bondProperties)
                + ", bondCount="
                + this.bondCount
                + '}';
    }
}
