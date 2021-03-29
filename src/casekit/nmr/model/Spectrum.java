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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Spectrum {

    private String[] nuclei;
    /**
     * An arbitrary name or description that can be assigned to this spectrum for identification purposes.
     */
    private String description;
    /**
     * An arbitrary name to identify the type of this spectrum, like COSY, NOESY, HSQC, etc. I
     * decided not to provide static Strings with given experiment type since the there are
     * numerous experiments yielding basically identical information having different names
     */
    private String specType;
    /**
     * The proton frequency of the spectrometer used to record this spectrum.
     */
    private Double spectrometerFrequency;
    private String solvent;
    private String standard;
    private List<Signal> signals;
    private int signalCount;

    public Spectrum() {
    }

    public Spectrum(final String[] nuclei) {
        this.nuclei = nuclei;
        this.signals = new ArrayList<>();
        this.signalCount = 0;
    }

    public Spectrum(final String[] nuclei, final String description, final String specType,
                    final Double spectrometerFrequency, final String solvent, final String standard,
                    final List<Signal> signals, final int signalCount) {
        this.nuclei = nuclei;
        this.description = description;
        this.specType = specType;
        this.spectrometerFrequency = spectrometerFrequency;
        this.solvent = solvent;
        this.standard = standard;
        this.signals = new ArrayList<>(signals);
        this.signalCount = signalCount;
    }

    public String[] getNuclei() {
        return this.nuclei;
    }

    public void setNuclei(final String[] nuclei) {
        this.nuclei = nuclei;
    }

    public int getNDim() {
        return this.getNuclei().length;
    }

    public boolean containsDim(final int dim) {
        return dim
                >= 0
                && dim
                <= this.getNDim();
    }

    public boolean compareNuclei(final String[] nuclei) {
        return Arrays.equals(this.getNuclei(), nuclei);
    }

    public void setSpecType(final String specType) {
        this.specType = specType;
    }

    public String getSpecType() {
        return this.specType;
    }

    public void setSpecDescription(final String description) {
        this.description = description;
    }

    public String getSpecDescription() {
        return this.description;
    }

    public int getSignalCount() {
        return this.signalCount;
    }

    public int getSignalCountWithEquivalences() {
        int sum = 0;
        for (final Signal signal : this.getSignals()) {
            sum += signal.getEquivalencesCount();
        }
        return sum;
    }

    /**
     * Adds a signal to this spectrum.
     *
     * @param signal signal to add
     *
     * @return
     */
    public boolean addSignal(final Signal signal) {
        return this.addSignal(signal, 0.0);
    }

    /**
     * Adds a signal to this spectrum and stores an equivalent signal index.
     *
     * @param signal        signal to add
     * @param pickPrecision precision to find equivalent signals to store in
     *
     * @return
     */
    public boolean addSignal(final Signal signal, final double pickPrecision) {
        if ((signal
                == null)
                || !this.compareNuclei(signal.getNuclei())) {
            return false;
        }

        // check for equivalent signals in all dimensions
        final List<Integer> closestSignalList = this.pickClosestSignal(signal.getShift(0), 0, pickPrecision);
        for (int dim = 1; dim
                < this.getNDim(); dim++) {
            closestSignalList.retainAll(this.pickClosestSignal(signal.getShift(dim), dim, pickPrecision));
        }

        if (closestSignalList.isEmpty()) {
            // add signal at the end of signal list
            this.signals.add(signal);
            this.signalCount++;
        } else {
            Signal closestSignal;
            for (final Integer closestSignalIndex : closestSignalList) {
                closestSignal = this.getSignal(closestSignalIndex);
                if (closestSignal.getMultiplicity()
                                 .equals(signal.getMultiplicity())) {
                    closestSignal.setEquivalencesCount(closestSignal.getEquivalencesCount()
                                                               + 1);
                }
            }
        }

        return true;

    }

    public boolean removeSignal(final Signal signal) {
        return this.removeSignal(this.getSignalIndex(signal));
    }

    public boolean removeSignal(final int signalIndex) {
        if (!this.checkSignalIndex(signalIndex)) {
            return false;
        }
        if (this.signals.remove(signalIndex)
                != null) {
            this.signalCount--;

            return true;
        }

        return false;
    }

    private boolean checkSignalIndex(final Integer signalIndex) {
        return (signalIndex
                != null)
                && (signalIndex
                >= 0)
                && (signalIndex
                < this.getSignalCount());
    }

    /**
     * Returns a NMR Signal at position number in the List
     *
     * @param signalIndex
     *
     * @return
     */
    public Signal getSignal(final int signalIndex) {
        if (!this.checkSignalIndex(signalIndex)) {
            return null;
        }

        return this.signals.get(signalIndex);
    }

    public List<Signal> getSignals() {
        return this.signals;
    }

    public Double getShift(final int signalIndex, final int dim) {
        if (!this.checkSignalIndex(signalIndex)) {
            return null;
        }

        return this.getSignal(signalIndex)
                   .getShift(dim);
    }

    public List<Double> getShifts(final int dim) {
        return this.getSignals()
                   .stream()
                   .map(signal -> signal.getShift(dim))
                   .collect(Collectors.toList());
    }

    public String getMultiplicity(final int signalIndex) {
        if (!this.checkSignalIndex(signalIndex)) {
            return null;
        }

        return this.getSignal(signalIndex)
                   .getMultiplicity();
    }

    public Boolean hasEquivalences(final int signalIndex) {
        if (!this.checkSignalIndex(signalIndex)) {
            return null;
        }

        return this.getEquivalencesCount(signalIndex)
                > 1;
    }

    public Integer getEquivalencesCount(final int signalIndex) {
        if (!this.checkSignalIndex(signalIndex)) {
            return null;
        }

        return this.getSignal(signalIndex)
                   .getEquivalencesCount();
    }

    public List<Integer> getEquivalencesCounts() {
        return this.getSignals()
                   .stream()
                   .map(Signal::getEquivalencesCount)
                   .collect(Collectors.toList());
    }

    /**
     * Returns the position of an NMRSignal the List
     *
     * @param signal
     *
     * @return
     */
    public int getSignalIndex(final Signal signal) {
        for (int s = 0; s
                < this.signals.size(); s++) {
            if (this.signals.get(s)
                    == signal) {
                return s;
            }
        }
        return -1;
    }

    public void setSpectrometerFrequency(final Double sf) {
        this.spectrometerFrequency = sf;
    }

    public Double getSpectrometerFrequency() {
        return spectrometerFrequency;
    }

    public void setSolvent(final String solvent) {
        this.solvent = solvent;
    }

    public String getSolvent() {
        return solvent;
    }

    public void setStandard(final String standard) {
        this.standard = standard;
    }

    public String getStandard() {
        return standard;
    }


    /**
     * Returns the signal index (or indices) closest to the given shift. If no signal is found within the interval
     * defined by {@code pickPrecision}, an empty list is returned.
     *
     * @param shift         query shift
     * @param dim           dimension in spectrum to look in
     * @param pickPrecision tolerance value for search window
     *
     * @return
     */
    public List<Integer> pickClosestSignal(final double shift, final int dim, final double pickPrecision) {
        final List<Integer> matchIndices = new ArrayList<>();
        if (!this.containsDim(dim)) {
            return matchIndices;
        }
        double minDiff = pickPrecision;
        // detect the minimal difference between a signal shift to the given query shift
        for (int s = 0; s
                < this.getSignalCount(); s++) {
            if (Math.abs(this.getShift(s, dim)
                                 - shift)
                    < minDiff) {
                minDiff = Math.abs(this.getShift(s, dim)
                                           - shift);
            }
        }
        for (int s = 0; s
                < this.getSignalCount(); s++) {
            if (Math.abs(this.getShift(s, dim)
                                 - shift)
                    == minDiff) {
                matchIndices.add(s);
            }
        }

        return matchIndices;
    }

    /**
     * Returns a list of signal indices within the interval defined by
     * pickPrecision. That list is sorted by the distances to the query shift.
     * If none is found an empty ArrayList is returned.
     *
     * @param shift         query shift
     * @param dim           dimension in spectrum to look in
     * @param pickPrecision tolerance value for search window
     *
     * @return
     */
    public List<Integer> pickSignals(final Double shift, final int dim, final double pickPrecision) {
        final List<Integer> pickedSignals = new ArrayList<>();
        if (!this.containsDim(dim)) {
            return pickedSignals;
        }
        for (int s = 0; s
                < this.getSignalCount(); s++) {
            if (Math.abs(this.getShift(s, dim)
                                 - shift)
                    <= pickPrecision) {
                pickedSignals.add(s);
            }
        }
        // sort signal indices by distance to query shift
        pickedSignals.sort(Comparator.comparingDouble(pickedSignalIndex -> Math.abs(shift
                                                                                            - this.getShift(
                pickedSignalIndex, dim))));

        return pickedSignals;
    }

    public Spectrum buildClone() {
        final Spectrum clone = new Spectrum(this.getNuclei());
        for (int i = 0; i
                < this.getSignalCount(); i++) {
            clone.addSignal(this.getSignal(i)
                                .buildClone());
        }
        clone.setSpecDescription(this.description);
        clone.setSolvent(this.solvent);
        clone.setSpecType(this.specType);
        clone.setSpectrometerFrequency(this.spectrometerFrequency);
        clone.setStandard(this.standard);

        return clone;
    }

    @Override
    public String toString() {
        return "Spectrum{"
                + "description='"
                + description
                + '\''
                + ", specType='"
                + specType
                + '\''
                + ", spectrometerFrequency="
                + spectrometerFrequency
                + ", solvent='"
                + solvent
                + '\''
                + ", standard='"
                + standard
                + '\''
                + ", signals="
                + signals
                + ", signalCount="
                + signalCount
                + '}';
    }

    public String getDescription() {
        return description;
    }
}
