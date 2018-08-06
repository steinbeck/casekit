/*
* This class was copied and modified from NMRSignal class in casekit.model package (by Christoph Steinbeck)
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

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Signal {
    
    private final int ndim;

    /**
     * Am array of doubles to store the chemical shift of
     */
    private final Double[] shifts;
    private final String[] nuclei;
    private final int[] assignedAtomIndices;

    /* Signal intensity in arbitrary values */
    private Double intensity;
    private String multiplicity;

    private Integer phase;
    public final static int PHASE_NONE = 0, PHASE_POSITIVE = 1, PHASE_NEGATIVE = 2;
    public final static String[] PHASENAMES = {"NONE", "POSITIVE", "NEGATIVE"};

    public Signal(final String[] nuclei, final Double[] shifts) {
        this.nuclei = nuclei;
        this.ndim = this.nuclei.length;
        this.shifts = shifts;
        this.assignedAtomIndices = this.initAssignedAtomIndices(this.ndim);
    }
    
    public Signal(final String[] nuclei, final Double[] shifts, final Double intensity) {
        this.nuclei = nuclei;
        this.ndim = this.nuclei.length;
        this.shifts = this.initShifts(shifts, this.ndim);
        this.assignedAtomIndices = this.initAssignedAtomIndices(this.ndim);
        this.intensity = intensity;
    }
    
    private Double[] initShifts(final Double[] shifts, final int ndim){
        final Double[] tempShifts = new Double[ndim];
        for (int d = 0; d < ndim; d++) {
            tempShifts[d] = shifts[d];
        }
        
        return tempShifts;
    }
    
    private int[] initAssignedAtomIndices(final int ndim){
        final int[] tempAssignedAtomIndices = new int[ndim];
        for (int d = 0; d < this.ndim; d++) {
            tempAssignedAtomIndices[d] = -1;
        }
        
        return tempAssignedAtomIndices;
    }
    
    public int getDim(){
        return this.ndim;
    }
    
    public String[] getNuclei(){
        return this.nuclei;
    }
    
    public void setShift(final Double shift, final int dim) {
        this.shifts[dim] = shift;
    }

    public Double getShift(final int dim) {
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
    
    public boolean setAssignedAtomIndices(final int[] indices){
        if(indices.length != this.ndim){
            return false;
        }
        for (int d = 0; d < this.ndim; d++) {
            this.assignedAtomIndices[d] = indices[d];
        }
        
        return true;
    }
    
    public int[] getAssignedIndices(){
        return this.assignedAtomIndices;
    }
    
    public boolean setAssignedAtomIndex(final int index, final int dim){
        if(dim < 0 || dim >= this.ndim){
            return false;
        }
        this.assignedAtomIndices[dim] = index;
        
        return true;
    }
    
    public int getAssignedAtomIndex(final int dim){
        return this.assignedAtomIndices[dim];
    }

    public void setPhase(final int phase) {
        this.phase = phase;
    }
    
    public Integer getPhase() {
        return this.phase;
    }

    @Override
    public String toString() {
        String s = "";
        s += ndim + " -dimensional NMRSignal for nuclei ";
        for (int f = 0; f < this.nuclei.length; f++) {
            s += this.nuclei[f] + "; ";
        }
        s += "\nShiftlist: ";
        for (int f = 0; f < this.shifts.length; f++) {
            s += this.shifts[f] + "; ";
        }
        s += "\n\n";
        return s;
    }
    
    /**
     *
     * @return
     * @deprecated 
     */
    public Signal getClone(){
        final Signal signalClone = new Signal(this.nuclei, this.shifts);
        signalClone.setIntensity(this.intensity);
        signalClone.setMultiplicity(this.multiplicity);
        signalClone.setPhase(this.phase);
        signalClone.setAssignedAtomIndices(this.assignedAtomIndices);
        
        return signalClone;
    }
    
}
