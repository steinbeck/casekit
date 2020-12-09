/*
 * The MIT License
 *
 * Copyright 2019 Michael Wenk [https://github.com/michaelwenk]
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package casekit.nmr.model.dimensional;

import java.util.Arrays;

public class Dimensional {

    private String[] nuclei;
    private int nDim;

    protected Dimensional() {
    }

    protected Dimensional(final String[] nuclei) {
        this.nuclei = nuclei;
        this.nDim = nuclei.length;
    }

    public final String[] getNuclei(){
        return this.nuclei;
    }

    public void setNuclei(String[] nuclei) {
        this.nuclei = nuclei;
    }

    public boolean compareNuclei(final String[] nuclei){
        return Arrays.equals(this.getNuclei(), nuclei);
    }

    /**
     * Returns the number of dimensions.
     *
     * @return
     */
    public final int getNDim() {
        return this.nDim;
    }

    public void setNDim(int nDim) {
        this.nDim = nDim;
    }

    /**
     * Checks whether the input dimension exists by dimension number. The dimension
     * indexing starts at 0.
     *
     * @param dim input dimension number
     * @return
     */
    public final boolean containsDim(final int dim){
        return (dim >= 0) && (dim < this.getNDim());
    }

    /**
     * Checks whether the input dimension count is equal to the number of dimensions of this object.
     *
     * @param nDim number of input dimensions
     * @return
     */
    public final boolean compareNDim(final int nDim){
        return nDim == this.getNDim();
    }
}
