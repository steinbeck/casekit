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
public class SignalCompact {

    private String[] strings; // nucleus dim 1, nucleus dim 2, ... , multiplicity, signal kind
    private Double[] doubles; // shift dim 1, shift dim 2, ... , intensity
    private Integer[] integers; // dimensions, equivalence count, phase

    public SignalCompact(final Signal signal) {
        this.strings = new String[signal.getNDim()
                + 3];
        this.doubles = new Double[signal.getNDim()
                + 1];
        for (int dim = 0; dim
                < signal.getNDim(); dim++) {
            this.strings[dim] = signal.getNuclei()[dim];
            this.doubles[dim] = signal.getShift(dim);
        }
        this.strings[signal.getNDim()] = signal.getMultiplicity();
        this.strings[signal.getNDim()
                + 1] = signal.getKind();
        this.strings[signal.getNDim()
                + 2] = signal.getId();
        this.doubles[signal.getNDim()] = signal.getIntensity();
        this.integers = new Integer[]{signal.getNDim(), signal.getEquivalencesCount(), signal.getPhase()};
    }

    public int dimensions() {
        return this.integers[0];
    }

    public String[] nuclei() {
        final int nDim = this.integers[0];
        final String[] nuclei = new String[nDim];
        for (int dim = 0; dim
                < nDim; dim++) {
            nuclei[dim] = this.strings[dim];
        }
        return nuclei;
    }

    public Signal toSignal() {
        final Signal signal = new Signal();
        signal.setNuclei(this.nuclei());
        signal.setMultiplicity(this.strings[this.dimensions()]);
        signal.setKind(this.strings[this.dimensions()
                + 1]);
        // @TODO remove following condition if not needed anymore (sherlock dataset storage is currently without signal ID)
        signal.setId(this.strings.length
                             - this.dimensions()
                             >= 3
                     ? this.strings[this.dimensions()
                + 2]
                     : null);
        signal.setShifts(new Double[this.dimensions()]);
        for (int dim = 0; dim
                < this.dimensions(); dim++) {
            signal.setShift(this.doubles[dim], dim);
        }
        signal.setIntensity(this.doubles[this.dimensions()]);
        signal.setEquivalencesCount(this.integers[1]);
        signal.setPhase(this.integers[2]);

        return signal;
    }


    public SignalCompact buildClone() {
        return new SignalCompact(this.strings.clone(), this.doubles.clone(), this.integers.clone());
    }

    @Override
    public String toString() {
        return "SignalCompact{"
                + "strings="
                + Arrays.toString(this.strings)
                + ", doubles="
                + Arrays.toString(this.doubles)
                + ", integers="
                + Arrays.toString(this.integers)
                + '}';
    }
}

