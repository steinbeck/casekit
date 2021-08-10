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
        this.multiplicities.add("s");
        this.multiplicities.add("d");
        this.multiplicities.add("t");
        this.multiplicities.add("q");
        this.multiplicities.add("unknown");
        this.minLimit = -50;
        this.maxLimit = 300;
        this.stepSize = 5;
        this.updateSteps(); // ppm range from -20 to 260 in 5 ppm steps
    }

    /**
     * Resets to following default values: <br> <br>
     * multiplicties: S, D, T. Q <br>
     * min. ppm limit: -50 <br>
     * max. ppm limit: 300 <br>
     * step size: 5
     */
    public void reset() {
        this.init();
    }

    public Map<String, List<Integer>> buildMultiplicitySections(final Spectrum spectrum, final int dim) {
        final Map<String, List<Integer>> multiplicitySections = new HashMap<>();
        // init
        for (final String multiplicity : this.multiplicities) {
            multiplicitySections.put(multiplicity, new ArrayList<>());
        }
        // set the mult. sections
        Signal signal;
        Integer shiftSection;
        String multiplicity;
        for (int i = 0; i
                < spectrum.getSignalCount(); i++) {
            signal = spectrum.getSignal(i);
            shiftSection = this.calculateShiftSection(signal, dim);
            if (shiftSection
                    == null) {
                System.err.println("MultiplicitySectionsBuilder: signal or its chemical shift is missing: "
                                           + signal);
                continue;
            }
            multiplicity = this.checkMultiplicity(signal);
            if (multiplicity
                    == null) {
                System.err.println("MultiplicitySectionsBuilder: signal multiplicity is not in list: "
                                           + signal);
                continue;
            }
            multiplicitySections.get(multiplicity)
                                .add(shiftSection);
        }

        return multiplicitySections;
    }

    public Integer calculateShiftSection(final Signal signal, final int dim) {
        if (signal
                == null
                || signal.getShift(dim)
                == null) {
            return null;
        }
        return (int) ((signal.getShift(dim)
                - this.minLimit)
                / this.stepSize);
    }

    public String checkMultiplicity(final Signal signal) {
        final String multiplicity = signal.getMultiplicity()
                                            != null
                                    ? signal.getMultiplicity()
                                    : "unknown";
        if (!this.multiplicities.contains(multiplicity)) {
            return null;
        }

        return multiplicity;
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
        this.updateSteps();
    }

    public int getMaxLimit() {
        return this.maxLimit;
    }

    public void setMaxLimit(final int maxLimit) {
        this.maxLimit = maxLimit;
        this.updateSteps();
    }

    public int getStepSize() {
        return this.stepSize;
    }

    public void setStepSize(final int stepSize) {
        this.stepSize = stepSize;
        this.updateSteps();
    }

    public int getSteps() {
        return this.steps;
    }

    private void updateSteps() {
        this.steps = this.calculateSteps(this.minLimit, this.maxLimit, this.stepSize);
    }

    public int calculateSteps(final int minLimit, final int maxLimit, final int stepSize) {
        return (maxLimit
                - minLimit)
                / stepSize;
    }
}
