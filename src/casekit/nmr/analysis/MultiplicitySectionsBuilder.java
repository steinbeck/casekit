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
package casekit.nmr.analysis;

import casekit.nmr.model.Signal;
import casekit.nmr.model.Spectrum;
import org.openscience.cdk.exception.CDKException;

import java.util.*;

/**
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class MultiplicitySectionsBuilder {

    private final HashSet<String> multiplicities;
    private int minLimit, maxLimit, stepSize, steps;

    public MultiplicitySectionsBuilder() {
        this.multiplicities = new HashSet<>();
        this.init();
    }

    private void init() {
        this.multiplicities.clear();
        this.multiplicities.add("S");
        this.multiplicities.add("D");
        this.multiplicities.add("T");
        this.multiplicities.add("Q");
        this.minLimit = -20;
        this.maxLimit = 260;
        this.stepSize = 5;
        this.steps = (this.maxLimit - this.minLimit) / this.stepSize; // ppm range from -20 to 260 in 5 ppm steps
    }

    /**
     * Resets to following default values: <br> <br>
     * multiplicties: S, D, T. Q <br>
     * min. ppm limit: -20 <br>
     * max. ppm limit: 260 <br>
     * step size: 5
     */
    public void reset() {
        this.init();
    }

    public Map<String, ArrayList<Integer>> buildMultiplicitySections(final Spectrum spectrum) throws CDKException {
        final HashMap<String, ArrayList<Integer>> multSections = new HashMap<>();
        // init
        for (final String mult : this.multiplicities) {
            multSections.put(mult, new ArrayList<>());
        }
        // set the mult. sections 
        Signal signal;
        int shiftSection;
        for (int i = 0; i < spectrum.getSignalCount(); i++) {
            signal = spectrum.getSignal(i);
            if ((signal == null) || (signal.getShift(0) == null) || (signal.getMultiplicity() == null) || (signal.getIntensity() == null) || (!this.multiplicities.contains(signal.getMultiplicity()))) {
                throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": signal, shift or multiplicity is missing");
            }
            shiftSection = (int) ((signal.getShift(0) - this.minLimit) / this.stepSize);
            multSections.get(signal.getMultiplicity()).add(shiftSection);
        }

        return multSections;
    }

    public Set<String> getMultiplicities() {
        return this.multiplicities;
    }

    public boolean addMultiplicity(final String mult) {
        return this.multiplicities.add(mult);
    }

    public boolean removeMultiplicity(final String mult) {
        return this.multiplicities.remove(mult);
    }

    public boolean containsMultiplicity(final String mult) {
        return this.multiplicities.contains(mult);
    }

    public int getMinLimit() {
        return this.minLimit;
    }

    public void setMinLimit(final int minLimit) {
        this.minLimit = minLimit;
    }

    public int getMaxLimit() {
        return this.maxLimit;
    }

    public void setMaxLimit(final int maxLimit) {
        this.maxLimit = maxLimit;
    }

    public int getStepSize() {
        return this.stepSize;
    }

    public void setStepSize(final int stepSize) {
        this.stepSize = stepSize;
    }
}
