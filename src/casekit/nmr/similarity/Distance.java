package casekit.nmr.similarity;

public class Distance {

    private final int signalIndexSpectrum1;
    private final int signalIndexSpectrum2;
    private final double value;

    public Distance(final int signalIndexSpectrum1, final int signalIndexSpectrum2, final double value) {
        this.signalIndexSpectrum1 = signalIndexSpectrum1;
        this.signalIndexSpectrum2 = signalIndexSpectrum2;
        this.value = value;
    }

    public int getSignalIndexSpectrum1() {
        return this.signalIndexSpectrum1;
    }

    public int getSignalIndexSpectrum2() {
        return this.signalIndexSpectrum2;
    }

    public double getValue() {
        return this.value;
    }

    @Override
    public String toString() {
        return "Distance{"
                + "signalIndexSpectrum1="
                + this.signalIndexSpectrum1
                + ", signalIndexSpectrum2="
                + this.signalIndexSpectrum2
                + ", value="
                + this.value
                + '}';
    }
}
