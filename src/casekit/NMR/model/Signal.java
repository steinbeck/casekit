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
    
    private final int nDim;

    /**
     * Am array of doubles to store the chemical shift of
     */
    private Double[] shifts;
    private final String[] nuclei;

    /* Signal intensity in arbitrary values */
    private Double intensity;
    private String multiplicity;

    private Integer phase;
    public final static int PHASE_NONE = 0, PHASE_POSITIVE = 1, PHASE_NEGATIVE = 2;
    public final static String[] PHASENAMES = {"NONE", "POSITIVE", "NEGATIVE"};

    
    public Signal(final String[] nuclei) {
        this.nuclei = nuclei;
        this.nDim = this.nuclei.length;
        this.shifts = this.initShifts(null, this.nDim);
    }
    
    public Signal(final String[] nuclei, final Double[] shifts) {
        this.nuclei = nuclei;
        this.nDim = this.nuclei.length;
        this.shifts = this.initShifts(shifts, this.nDim);
    }
    
    public Signal(final String[] nuclei, final Double[] shifts, final Double intensity) {
        this(nuclei, shifts);
        this.intensity = intensity;
    }
    
    private Double[] initShifts(final Double[] shifts, final int nDim){
        final Double[] tempShifts = new Double[nDim];
        for (int d = 0; d < nDim; d++) {
            if((shifts != null) && (shifts.length == nDim)){
                tempShifts[d] = shifts[d];
            } else {
                tempShifts[d] = null;
            }
        }
        
        return tempShifts;
    }
    
    public int getDimCount(){
        return this.nDim;
    }
    
    public String[] getNuclei(){
        return this.nuclei;
    }
    
    public boolean setShift(final Double shift, final int dim) {
        if(!this.checkDimension(dim)){
            return false;
        }
        this.shifts[dim] = shift;
        
        return true;
    }

    public Double getShift(final int dim) {
        if(!this.checkDimension(dim)){
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

    public void setPhase(final Integer phase) {
        this.phase = phase;
    }
    
    public Integer getPhase() {
        return this.phase;
    }

    public boolean checkDimension(final int dim){
       return (dim >= 0) && (dim < this.nDim);
   }
    
    /**
     *
     * @return
     */
    public Signal getClone(){
        final Signal signalClone = new Signal(this.nuclei, this.shifts);
        signalClone.setIntensity(this.intensity);
        signalClone.setMultiplicity(this.multiplicity);
        signalClone.setPhase(this.phase);
        
        return signalClone;
    }
    
}
