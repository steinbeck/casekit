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
public class Signal {

    private String[] nuclei;
    private Double[] shifts;
    private String multiplicity;
    private String kind;
    private Double intensity;
    private int equivalencesCount;
    private Integer phase;


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

    public Signal buildClone() {
        return new Signal(this.getNuclei()
                              .clone(), this.shifts.clone(), this.multiplicity, this.kind, this.intensity,
                          this.equivalencesCount, this.phase);
    }

    @Override
    public String toString() {
        return "Signal{"
                + "nuclei="
                + Arrays.toString(this.nuclei)
                + ", shifts="
                + Arrays.toString(this.shifts)
                + ", multiplicity='"
                + this.multiplicity
                + '\''
                + ", kind='"
                + this.kind
                + '\''
                + ", intensity="
                + this.intensity
                + ", equivalencesCount="
                + this.equivalencesCount
                + ", phase="
                + this.phase
                + '}';
    }
}
