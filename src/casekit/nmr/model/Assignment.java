/*
 * The MIT License
 *
 * Copyright 2018 Michael Wenk [https://github.com/michaelwenk].
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package casekit.nmr.model;

import lombok.AllArgsConstructor;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import java.util.Arrays;

/**
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
@NoArgsConstructor
@AllArgsConstructor
@Getter
@Setter
public class Assignment
        implements Cloneable {

    private String[] nuclei;
    private int[][][] assignments;


    public void initAssignments(final int length) {
        final int[][][] temp = new int[this.getNDim()][length][0];
        for (int i = 0; i
                < this.getNDim(); i++) {
            for (int j = 0; j
                    < length; j++) {
                temp[i][j] = new int[]{};
            }
        }
        this.assignments = temp;
    }

    public int getNDim() {
        return this.getNuclei().length;
    }

    public boolean containsDim(final int dim) {
        return dim
                >= 0
                && dim
                < this.getNDim();
    }

    public boolean compareNuclei(final String[] nuclei) {
        return Arrays.equals(this.getNuclei(), nuclei);
    }

    /**
     * Sets an assignment array with value for an index position.
     *
     * @param dim
     * @param index
     * @param assignment
     *
     * @return
     */
    public boolean setAssignment(final int dim, final int index, final int[] assignment) {
        if (!this.containsDim(dim)
                || !this.checkIndex(dim, index)) {
            return false;
        }
        this.assignments[dim][index] = assignment;

        return true;
    }

    public boolean setAssignments(final int dim, final int[][] assignments) {
        if (!this.containsDim(dim)
                || this.getSize()
                != assignments.length) {
            return false;
        }
        for (int i = 0; i
                < this.getSize(); i++) {
            this.setAssignment(dim, i, assignments[i]);
        }

        return true;
    }

    public int[] getAssignment(final int dim, final int index) {
        if (!this.containsDim(dim)
                || !this.checkIndex(dim, index)) {
            return null;
        }

        return this.assignments[dim][index];
    }

    public int getAssignment(final int dim, final int index, final int equivalenceIndex) {
        if (!this.containsDim(dim)
                || !this.checkIndex(dim, index)) {
            return -1;
        }

        return this.assignments[dim][index][equivalenceIndex];
    }

    public void addAssignmentEquivalence(final int dim, final int index, final int equivalenceIndex) {
        final int[] temp = this.getAssignment(dim, index);
        final int[] equivalenceIndices = new int[temp.length
                + 1];
        for (int j = 0; j
                < temp.length; j++) {
            equivalenceIndices[j] = temp[j];
        }
        equivalenceIndices[equivalenceIndices.length
                - 1] = equivalenceIndex;

        this.setAssignment(dim, index, equivalenceIndices);
    }

    public Integer getIndex(final int dim, final int assignment) {
        if (!this.containsDim(dim)) {
            return null;
        }
        for (int index = 0; index
                < this.assignments[dim].length; index++) {
            if (Arrays.stream(this.getAssignment(dim, index))
                      .anyMatch(value -> value
                              == assignment)) {
                return index;
            }
        }

        return -1;
    }

    public int[][] getAssignments(final int dim) {
        if (!this.containsDim(dim)) {
            return null;
        }

        return this.assignments[dim];
    }

    public int getSize() {
        if (this.getNDim()
                > 0) {
            return this.assignments[0].length;
        }
        return 0;
    }

    public int getSetAssignmentsCount(final int dim) {
        int setAssignmentsCounter = 0;
        if (this.containsDim(dim)) {
            for (int j = 0; j
                    < this.assignments[dim].length; j++) {
                if (this.assignments[dim][j].length
                        > 0) {
                    setAssignmentsCounter++;
                }
            }
        }
        return setAssignmentsCounter;
    }

    public int getSetAssignmentsCountWithEquivalences(final int dim) {
        int setAssignmentsCounter = 0;
        if (this.containsDim(dim)) {
            for (int j = 0; j
                    < this.assignments[dim].length; j++) {
                setAssignmentsCounter += this.assignments[dim][j].length;
            }
        }
        return setAssignmentsCounter;
    }

    public boolean addAssignment(final int dim, final int[] assignment) {
        if (!this.containsDim(dim)) {
            return false;
        }

        final int[][][] newAssignments = new int[this.getNDim()][][];
        for (int d = 0; d
                < this.getNDim(); d++) {
            newAssignments[d] = new int[this.assignments[d].length
                    + 1][];
            for (int i = 0; i
                    < this.assignments[d].length; i++) {
                newAssignments[d][i] = this.assignments[d][i];
            }
        }
        newAssignments[dim][this.assignments[dim].length] = assignment;
        this.assignments = newAssignments;

        return true;
    }

    private boolean checkIndex(final int dim, final int index) {
        return (index
                >= 0)
                && (index
                < this.assignments[dim].length);
    }

    @Override
    public Assignment clone() throws CloneNotSupportedException {
        return (Assignment) super.clone();
    }

    @Override
    public String toString() {
        final StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append("Assignments:\n");

        for (int i = 0; i
                < this.getNDim(); i++) {
            stringBuilder.append(Arrays.toString(this.assignments[i]))
                         .append("\n");
        }

        return stringBuilder.toString();
    }
}
