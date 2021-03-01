/*
 * This class was adopted and modified from an earlier version by Christoph Steinbeck
 */

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

import casekit.nmr.model.dimensional.Dimensional;

import java.util.Arrays;

/**
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Signal extends Dimensional {

    private Double[] shifts;
    private String multiplicity;
    private Double intensity;
    private String kind;
    private int equivalencesCount;
    private int phase;


    public Signal() {
    }

    public Signal(final String[] nuclei) {
        this(nuclei, null, null, null, null, 1, 0);
    }

    public Signal(final String[] nuclei, final Double[] shifts, final String multiplicity, final String kind, final Double intensity, final int equivalencesCount, final int phase) {
        super(nuclei);
        this.shifts = shifts; // this.initShifts(shifts, this.getNDim());
        this.multiplicity = multiplicity;
        this.kind = kind;
        this.intensity = intensity;
        this.equivalencesCount = equivalencesCount;
        this.phase = phase;
    }

    //    private Double[] initShifts(final Double[] shifts, final int nDim) {
    ////        if((shifts == null) || (shifts.length != nDim)){
    ////            throw new Exception("Number of given nuclei (" + nDim + ") and shifts (" + shifts.length + ") is not the same!!!");
    ////        }
    //        final Double[] tempShifts = new Double[nDim];
    //        for (int d = 0; d < nDim; d++) {
    //            tempShifts[d] = shifts[d];
    //        }
    //
    //        return tempShifts;
    //    }

    public boolean setShift(final Double shift, final int dim) {
        if (!this.containsDim(dim)) {
            return false;
        }
        this.shifts[dim] = shift;

        return true;
    }

    public Double getShift(final int dim) {
        if (!this.containsDim(dim)) {
            return null;
        }
        return this.shifts[dim];
    }

    public void setIntensity(final Double intensity) {
        this.intensity = intensity;
    }

    public Double getIntensity() {
        return this.intensity;
    }

    public void setMultiplicity(final String multiplicity) {
        this.multiplicity = multiplicity;
    }

    public String getMultiplicity() {
        return this.multiplicity;
    }

    public String getKind() {
        return kind;
    }

    public void setKind(final String kind) {
        this.kind = kind;
    }

    //    public void setPhase(final Integer phase) {
    //        this.phase = phase;
    //    }
    //
    //    public Integer getPhase() {
    //        return this.phase;
    //    }


    public Signal buildClone() {
        return new Signal(this.getNuclei(), this.shifts, this.multiplicity, this.kind, this.intensity, this.equivalencesCount, this.phase);
    }

    @Override
    public String toString() {
        return "Signal{" + "shifts=" + Arrays.toString(shifts) + ", multiplicity='" + multiplicity + '\'' + ", intensity=" + intensity + ", kind='" + kind + '\'' + ", equivalencesCount=" + equivalencesCount + ", phase=" + phase + '}';
    }

    public Double[] getShifts() {
        return shifts;
    }

    public int getEquivalencesCount() {
        return equivalencesCount;
    }

    public void setEquivalencesCount(final int equivalencesCount) {
        this.equivalencesCount = equivalencesCount;
    }

    public int getPhase() {
        return phase;
    }

    public void setPhase(final int phase) {
        this.phase = phase;
    }
}
