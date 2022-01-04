package casekit.nmr.lsd.inputfile;

import casekit.nmr.lsd.Constants;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class Utilities {
    public static String buildSSTR(final int sstrIndex, final String atomType, final List<Integer> hybridization,
                                   final List<Integer> protonsCount) {
        if (hybridization.isEmpty()) {
            hybridization.addAll(Arrays.stream(Constants.defaultHybridizationMap.get(atomType))
                                       .boxed()
                                       .collect(Collectors.toList()));
        }
        if (protonsCount.isEmpty()) {
            protonsCount.addAll(Arrays.stream(Constants.defaultProtonsCountPerValencyMap.get(atomType))
                                      .boxed()
                                      .collect(Collectors.toList()));
        }
        final StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append("SSTR S")
                     .append(sstrIndex)
                     .append(" ")
                     .append(atomType)
                     .append(" ");
        if (hybridization.size()
                == 1) {
            stringBuilder.append(hybridization.get(0))
                         .append(" ");
            if (protonsCount.size()
                    == 1) {
                stringBuilder.append(protonsCount.get(0));
            } else {
                stringBuilder.append(buildMultipleValuesString(protonsCount));
            }
        } else {
            stringBuilder.append(buildMultipleValuesString(hybridization));
            stringBuilder.append(" ");
            if (protonsCount.size()
                    == 1) {
                stringBuilder.append(protonsCount.get(0));
            } else {
                stringBuilder.append(buildMultipleValuesString(protonsCount));
            }
        }

        return stringBuilder.toString();
    }

    private static String buildMultipleValuesString(final List<Integer> values) {
        final StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append("(");
        for (int l = 0; l
                < values.size(); l++) {
            stringBuilder.append(values.get(l));
            if (l
                    < values.size()
                    - 1) {
                stringBuilder.append(" ");
            }
        }
        stringBuilder.append(")");

        return stringBuilder.toString();
    }
}
