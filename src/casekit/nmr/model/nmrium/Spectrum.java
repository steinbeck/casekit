/*
 * MIT License
 *
 * Copyright (c) 2020 Michael Wenk (https://github.com/michaelwenk)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package casekit.nmr.model.nmrium;

import casekit.nmr.model.Signal;
import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import lombok.ToString;

import java.util.ArrayList;
import java.util.Map;

@NoArgsConstructor
@Getter
@Setter
@ToString
@JsonIgnoreProperties(ignoreUnknown = true)
public class Spectrum {

    private String id;
    private Default<Range> ranges;
    private Default<Zone> zones;
    private Map<String, Object> info;

    public casekit.nmr.model.Spectrum toSpectrum(final boolean considerSignalKind) {
        final int dimension = (int) this.info.get("dimension");
        final boolean isFid = (boolean) this.info.get("isFid");

        if (!isFid) {
            if (dimension
                    == 1) {
                final String nucleus = (String) this.info.get("nucleus");
                final casekit.nmr.model.Spectrum spectrum = new casekit.nmr.model.Spectrum();
                spectrum.setNuclei(new String[]{nucleus});
                this.ranges.getValues()
                           .forEach(range -> range.getSignal()
                                                  .forEach(signal1D -> {
                                                      if (considerSignalKind
                                                              && signal1D.getKind()
                                                                         .equals("signal")) {
                                                          spectrum.addSignal(new Signal(new String[]{nucleus},
                                                                                        new Double[]{
                                                                                                signal1D.getDelta()},
                                                                                        signal1D.getMultiplicity(),
                                                                                        signal1D.getKind(), null, 0, 0,
                                                                                        null));
                                                      }
                                                  }));
                spectrum.addMetaInfo("solvent", (String) this.info.get("solvent"));
                spectrum.addMetaInfo("spectrumType", (String) this.info.get("experiment"));

                return spectrum;

            } else if (dimension
                    == 2) {
                final String[] nuclei = ((ArrayList<String>) this.info.get("nucleus")).toArray(new String[]{});
                final casekit.nmr.model.Spectrum spectrum = new casekit.nmr.model.Spectrum();
                spectrum.setNuclei(nuclei);

                this.zones.getValues()
                          .forEach(zone -> zone.getSignal()
                                               .forEach(signal2D -> {
                                                   if (considerSignalKind
                                                           && signal2D.getKind()
                                                                      .equals("signal")) {
                                                       spectrum.addSignal(new Signal(nuclei, new Double[]{
                                                               (Double) signal2D.getX()
                                                                                .get("delta"), (Double) signal2D.getY()
                                                                                                                .get("delta")},
                                                                                     signal2D.getMultiplicity(),
                                                                                     signal2D.getKind(), null, 0, 0,
                                                                                     signal2D.getPathLength()));
                                                   }
                                               }));
                spectrum.addMetaInfo("solvent", (String) this.info.get("solvent"));
                spectrum.addMetaInfo("spectrumType", (String) this.info.get("experiment"));

                return spectrum;
            }
        }

        return null;
    }
}
