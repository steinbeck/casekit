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
package casekit.NMR.model;

import casekit.NMR.model.dimensional.DimensionalNMR;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Signal extends DimensionalNMR {

    /**
     * Am array of doubles to store the chemical shift of
     */
    private Double[] shifts;

    private Double intensity;
    private String multiplicity;

//    private Integer phase;
//    public final static int PHASE_NONE = 0, PHASE_POSITIVE = 1, PHASE_NEGATIVE = 2;
//    public final static String[] PHASENAMES = {"NONE", "POSITIVE", "NEGATIVE"};

    
    public Signal(final String[] nuclei) throws Exception {
        this(nuclei, null);
    }
    
    public Signal(final String[] nuclei, final Double[] shifts) throws Exception {
        this(nuclei, shifts, null, null);
    }

    public Signal(final String[] nuclei, final Double[] shifts, final String multiplicity, final Double intensity) throws Exception {
        super(nuclei);
        this.shifts = this.initShifts(shifts, this.getNDim());
        this.multiplicity = multiplicity;
        this.intensity = intensity;
    }
    
    private Double[] initShifts(final Double[] shifts, final int nDim) throws Exception {
        if((shifts == null) || (shifts.length != nDim)){
            throw new Exception("Number of given nuclei (" + nDim + ") and shifts (" + shifts.length + ") is not the same!!!");
        }
        final Double[] tempShifts = new Double[nDim];
        for (int d = 0; d < nDim; d++) {
            tempShifts[d] = shifts[d];
        }
        
        return tempShifts;
    }

    public boolean setShift(final Double shift, final int dim) {
        if(!this.containsDim(dim)){
            return false;
        }
        this.shifts[dim] = shift;
        
        return true;
    }

    public Double getShift(final int dim) {
        if(!this.containsDim(dim)){
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

//    public void setPhase(final Integer phase) {
//        this.phase = phase;
//    }
//
//    public Integer getPhase() {
//        return this.phase;
//    }


    public Signal getClone() throws Exception {
//        final Signal clone = new Signal(this.getDimNames(), this.shifts, this.multiplicity, this.intensity);
//        clone.setPhase(this.phase);
//
//        return clone;

        return new Signal(this.getNuclei(), this.shifts, this.multiplicity, this.intensity);
    }
    
}
