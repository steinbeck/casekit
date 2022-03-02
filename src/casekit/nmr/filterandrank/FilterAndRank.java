package casekit.nmr.filterandrank;

import casekit.nmr.analysis.MultiplicitySectionsBuilder;
import casekit.nmr.model.Assignment;
import casekit.nmr.model.DataSet;
import casekit.nmr.model.Spectrum;
import casekit.nmr.similarity.Similarity;
import casekit.nmr.utils.Statistics;
import org.openscience.cdk.fingerprint.BitSetFingerprint;

import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

public class FilterAndRank {

    public static List<DataSet> filterAndRank(final List<DataSet> dataSetList, final Spectrum querySpectrum,
                                              final double shiftTolerance, final double maxAverageDeviation,
                                              final boolean checkMultiplicity, final boolean checkEquivalencesCount,
                                              final MultiplicitySectionsBuilder multiplicitySectionsBuilder,
                                              final boolean allowIncompleteMatch) {
        return rank(filter(dataSetList, querySpectrum, shiftTolerance, maxAverageDeviation, checkMultiplicity,
                           checkEquivalencesCount, multiplicitySectionsBuilder, allowIncompleteMatch));
    }

    public static List<DataSet> filter(final List<DataSet> dataSetList, final Spectrum querySpectrum,
                                       final double shiftTolerance, final double maxAverageDeviation,
                                       final boolean checkMultiplicity, final boolean checkEquivalencesCount,
                                       final MultiplicitySectionsBuilder multiplicitySectionsBuilder,
                                       final boolean allowIncompleteMatch) {
        if (querySpectrum.getNDim()
                == 1
                && querySpectrum.getNuclei()[0].equals("13C")) {
            return dataSetList.stream()
                              .filter(dataSet -> checkDataSet(dataSet, querySpectrum, shiftTolerance,
                                                              maxAverageDeviation, checkMultiplicity,
                                                              checkEquivalencesCount, multiplicitySectionsBuilder,
                                                              allowIncompleteMatch)
                                      != null)
                              .collect(Collectors.toList());
        }

        return dataSetList;
    }

    public static DataSet checkDataSet(final DataSet dataSet, final Spectrum querySpectrum, final double shiftTolerance,
                                       final double maxAverageDeviation, final boolean checkMultiplicity,
                                       final boolean checkEquivalencesCount,
                                       final MultiplicitySectionsBuilder multiplicitySectionsBuilder,
                                       final boolean allowIncompleteMatch) {
        final BitSetFingerprint bitSetFingerprintQuerySpectrum = Similarity.getBitSetFingerprint(querySpectrum, 0,
                                                                                                 multiplicitySectionsBuilder);
        final Spectrum spectrum = dataSet.getSpectrum()
                                         .toSpectrum();
        final Assignment spectralMatchAssignment = Similarity.matchSpectra(spectrum, querySpectrum, 0, 0,
                                                                           shiftTolerance, checkMultiplicity,
                                                                           checkEquivalencesCount, false);

        dataSet.addAttachment("querySpectrumSignalCount", querySpectrum.getSignalCount());
        dataSet.addAttachment("querySpectrumSignalCountWithEquivalences",
                              querySpectrum.getSignalCountWithEquivalences());
        dataSet.addAttachment("setAssignmentsCountWithEquivalences",
                              spectralMatchAssignment.getSetAssignmentsCountWithEquivalences(0));
        final boolean isCompleteSpectralMatch = querySpectrum.getSignalCount()
                == spectralMatchAssignment.getSetAssignmentsCount(0);
        final boolean isCompleteSpectralMatchWithEquivalences = querySpectrum.getSignalCountWithEquivalences()
                == spectralMatchAssignment.getSetAssignmentsCountWithEquivalences(0);
        dataSet.addAttachment("setAssignmentsCount", spectralMatchAssignment.getSetAssignmentsCount(0));
        dataSet.addAttachment("setAssignmentsCountWithEquivalences",
                              spectralMatchAssignment.getSetAssignmentsCountWithEquivalences(0));
        dataSet.addAttachment("isCompleteSpectralMatch", isCompleteSpectralMatch);
        dataSet.addAttachment("isCompleteSpectralMatchWithEquivalences", isCompleteSpectralMatchWithEquivalences);
        dataSet.addAttachment("spectralMatchAssignment", spectralMatchAssignment);

        Double[] deviations = Similarity.getDeviations(spectrum, querySpectrum, 0, 0, spectralMatchAssignment);
        if (allowIncompleteMatch) {
            deviations = Arrays.stream(deviations)
                               .filter(Objects::nonNull)
                               .toArray(Double[]::new);
        }
        final Double averageDeviation = Statistics.calculateAverageDeviation(deviations);
        if (averageDeviation
                != null
                && averageDeviation
                <= maxAverageDeviation) {
            dataSet.addAttachment("averageDeviation", averageDeviation);
            final Double rmsd = Statistics.calculateRMSD(deviations);
            dataSet.addAttachment("rmsd", rmsd);

            final BitSetFingerprint bitSetFingerprintDataSet = Similarity.getBitSetFingerprint(spectrum, 0,
                                                                                               multiplicitySectionsBuilder);
            final Double tanimotoCoefficient = Similarity.calculateTanimotoCoefficient(bitSetFingerprintQuerySpectrum,
                                                                                       bitSetFingerprintDataSet);
            dataSet.addAttachment("tanimoto", tanimotoCoefficient);

            return dataSet;
        }

        return null;
    }

    public static List<DataSet> rank(final List<DataSet> dataSetList) {
        dataSetList.sort((dataSet1, dataSet2) -> {
            final int avgDevComparison = compareNumericDataSetAttachmentKey(dataSet1, dataSet2, "averageDeviation");
            if (avgDevComparison
                    != 0) {
                return avgDevComparison;
            }

            final boolean isCompleteSpectralMatchDataSet1 = Boolean.parseBoolean(dataSet1.getMeta()
                                                                                         .get("isCompleteSpectralMatch"));
            final boolean isCompleteSpectralMatchDataSet2 = Boolean.parseBoolean(dataSet2.getMeta()
                                                                                         .get("isCompleteSpectralMatch"));
            if (isCompleteSpectralMatchDataSet1
                    && !isCompleteSpectralMatchDataSet2) {
                return -1;
            } else if (!isCompleteSpectralMatchDataSet1
                    && isCompleteSpectralMatchDataSet2) {
                return 1;
            }
            final int setAssignmentsCountComparison = compareNumericDataSetAttachmentKey(dataSet1, dataSet2,
                                                                                         "setAssignmentsCount");
            if (setAssignmentsCountComparison
                    != 0) {
                return -1
                        * setAssignmentsCountComparison;
            }

            return 0;
        });

        return dataSetList;
    }

    private static int compareNumericDataSetAttachmentKey(final DataSet dataSet1, final DataSet dataSet2,
                                                          final String attachmentKey) {
        Double valueDataSet1 = null;
        Double valueDataSet2 = null;
        try {
            valueDataSet1 = Double.parseDouble(String.valueOf(dataSet1.getAttachment()
                                                                      .get(attachmentKey)));
        } catch (final NullPointerException | NumberFormatException e) {
            //                e.printStackTrace();
        }
        try {
            valueDataSet2 = Double.parseDouble(String.valueOf(dataSet2.getAttachment()
                                                                      .get(attachmentKey)));
        } catch (final NullPointerException | NumberFormatException e) {
            //                e.printStackTrace();
        }

        if (valueDataSet1
                != null
                && valueDataSet2
                != null) {
            if (valueDataSet1
                    < valueDataSet2) {
                return -1;
            } else if (valueDataSet1
                    > valueDataSet2) {
                return 1;
            }
            return 0;
        }
        if (valueDataSet1
                != null) {
            return -1;
        } else if (valueDataSet2
                != null) {
            return 1;
        }

        return 0;
    }
}
