package casekit.nmr.similarity;

import casekit.nmr.model.Signal;
import casekit.nmr.model.Spectrum;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public class Utilities {

    public static List<Distance> buildDistanceList(final Spectrum spectrum1, final Spectrum spectrum2, final int dim1,
                                                   final int dim2, final double shiftTolerance,
                                                   final boolean checkMultiplicity,
                                                   final boolean checkEquivalencesCount,
                                                   final boolean allowLowerEquivalencesCount) {
        final List<Distance> distanceList = new ArrayList<>();
        Double distanceValue;
        for (int i = 0; i
                < spectrum1.getSignalCount(); i++) {
            for (int j = 0; j
                    < spectrum2.getSignalCount(); j++) {
                distanceValue = getDistanceValue(spectrum1.getSignal(i), spectrum2.getSignal(j), dim1, dim2,
                                                 checkMultiplicity, checkEquivalencesCount, allowLowerEquivalencesCount,
                                                 shiftTolerance);
                if (distanceValue
                        != null) {
                    distanceList.add(new Distance(i, j, distanceValue));
                }
            }
        }
        distanceList.sort(Comparator.comparingDouble(Distance::getValue));

        return distanceList;
    }

    public static Double getDistanceValue(final Signal signal1, final Signal signal2, final int dim1, final int dim2,
                                          final boolean checkMultiplicity, final boolean checkEquivalencesCount,
                                          final boolean allowLowerEquivalencesCount, final double shiftTolerance) {
        boolean passed = true;
        // @TODO maybe consider further parameters to check ? e.g. intensity
        if (checkMultiplicity) {
            passed = (signal1.getMultiplicity()
                    == null
                    && signal2.getMultiplicity()
                    == null)
                    || (signal1.getMultiplicity()
                    != null
                    && signal1.getMultiplicity()
                              .equals(signal2.getMultiplicity()));
        }
        if (passed
                && checkEquivalencesCount) {
            if (allowLowerEquivalencesCount) {
                passed = signal1.getEquivalencesCount()
                        <= signal2.getEquivalencesCount();
            } else {
                passed = signal1.getEquivalencesCount()
                        == signal2.getEquivalencesCount();
            }
        }
        if (!passed) {
            return null;
        }
        final double distanceValue = Math.abs(signal1.getShift(dim1)
                                                      - signal2.getShift(dim2));

        return distanceValue
                       > shiftTolerance
               ? null
               : distanceValue;
    }
}
