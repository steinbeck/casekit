package casekit.nmr.model;

import lombok.AllArgsConstructor;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
@NoArgsConstructor
@AllArgsConstructor
@Getter
@Setter
public class SpectrumCompact {

    private String[] nuclei;
    private Map<String, String> meta;
    private SignalCompact[] signals;

    public SpectrumCompact(final Spectrum spectrum) {
        this.nuclei = spectrum.getNuclei();
        this.meta = spectrum.getMeta();
        this.signals = spectrum.getSignals()
                               .stream()
                               .map(SignalCompact::new)
                               .toArray(SignalCompact[]::new);
    }

    public void addMetaInfo(final String key, final String value) {
        if (this.meta
                == null) {
            this.meta = new HashMap<>();
        }
        this.meta.put(key, value);
    }

    public void removeMetaInfo(final String key) {
        this.meta.remove(key);
    }

    public Spectrum toSpectrum() {
        final Spectrum spectrum = new Spectrum();
        spectrum.setNuclei(this.nuclei);
        spectrum.setMeta(this.meta);
        spectrum.setSignals(Arrays.stream(this.signals)
                                  .map(SignalCompact::toSignal)
                                  .collect(Collectors.toList()));
        spectrum.setSignalCount(this.signals.length);

        return spectrum;
    }

    public SpectrumCompact buildClone() {
        final Map<String, String> metaTemp = this.meta
                                                     == null
                                             ? new HashMap<>()
                                             : new HashMap<>(this.meta);
        return new SpectrumCompact(this.nuclei.clone(), metaTemp, Arrays.stream(this.signals)
                                                                        .map(SignalCompact::buildClone)
                                                                        .toArray(SignalCompact[]::new));
    }

    @Override
    public String toString() {
        return "SpectrumCompact{"
                + "nuclei="
                + Arrays.toString(this.nuclei)
                + ", meta="
                + this.meta
                + ", signals="
                + Arrays.toString(this.signals)
                + '}';
    }
}
