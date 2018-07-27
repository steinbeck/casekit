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
package NMR;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Signal {
    
    private final int ndim;

    /**
     * Am array of doubles to store the chemical shift of
     */
    private final Double shift[];
    private final String[] nuclei;

    /* Signal intensity in arbitrary values */
    private Double intensity;
    private String multiplicity;

    private Integer phase;
    public final static int DIM_ONE = 1, DIM_TWO = 2, DIM_THREE = 3, DIM_FOUR = 4;
    public final static int SHIFT_PROTON = 0, SHIFT_HETERO = 1;
    public final static int PHASE_NONE = 0, PHASE_POSITIVE = 1, PHASE_NEGATIVE = 2;
    public final static String[] PHASENAMES = {"NONE", "POSITIVE", "NEGATIVE"};

    public Signal(final String[] nuclei, final Double[] shift, final Double intensity, final Integer phase, final String multiplicity) {
        this.nuclei = nuclei;
        this.ndim = this.nuclei.length;
        this.shift = shift;
        this.intensity = intensity;
        this.phase = phase;
        this.multiplicity = multiplicity;
    }
    
    public int getDim(){
        return this.ndim;
    }
    
    public void setShift(final Double shift, final int dim) {
        this.shift[dim] = shift;
    }

    public Double getShift(final int dim) {
        return this.shift[dim];
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
        for (int f = 0; f < this.shift.length; f++) {
            s += this.shift[f] + "; ";
        }
        s += "\n\n";
        return s;
    }
    
    public String[] getNuclei(){
        return this.nuclei;
    }
}
