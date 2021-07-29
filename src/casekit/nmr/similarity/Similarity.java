/*
 * The MIT License
 *
 * Copyright 2019 Michael Wenk [https://github.com/michaelwenk]
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package casekit.nmr.similarity;

import casekit.nmr.analysis.MultiplicitySectionsBuilder;
import casekit.nmr.model.Assignment;
import casekit.nmr.model.Signal;
import casekit.nmr.model.Spectrum;
import casekit.nmr.utils.Statistics;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.BitSetFingerprint;
import org.openscience.cdk.similarity.Tanimoto;

import java.util.*;

public class Similarity {

    /**
     * Checks whether two spectra contain given dimensions.
     *
     * @param spectrum1 first spectrum
     * @param spectrum2 second spectrum
     * @param dim1      dimension to select in first spectrum
     * @param dim2      dimension to select in second spectrum
     *
     * @return true if both spectra contain the selected dimension
     */
    private static boolean checkDimensions(final Spectrum spectrum1, final Spectrum spectrum2, final int dim1,
                                           final int dim2) {
        return spectrum1.containsDim(dim1)
                && spectrum2.containsDim(dim2);
    }

    /**
     * Calculates the continuous Tanimoto coefficient between two spectra in given dimensions.
     *
     * @param spectrum1 first spectrum
     * @param spectrum2 second spectrum
     * @param dim1      dimension in first spectrum to take the shifts from
     * @param dim2      dimension in second spectrum to take the shifts from
     *
     * @return
     */
    public static Double calculateTanimotoCoefficient(final Spectrum spectrum1, final Spectrum spectrum2,
                                                      final int dim1, final int dim2,
                                                      final MultiplicitySectionsBuilder multiplicitySectionsBuilder) {
        if (!checkDimensions(spectrum1, spectrum2, dim1, dim2)) {
            return null;
        }

        try {
            return Tanimoto.calculate(getFingerprint(spectrum1, dim1, multiplicitySectionsBuilder),
                                      getFingerprint(spectrum2, dim2, multiplicitySectionsBuilder));
        } catch (final CDKException e) {
            e.printStackTrace();
        }

        return null;
    }

    public static BitSetFingerprint getFingerprint(final Spectrum spectrum, final int dim,
                                                   final MultiplicitySectionsBuilder multiplicitySectionsBuilder) throws CDKException {
        final BitSetFingerprint bitSetFingerprint = new BitSetFingerprint(multiplicitySectionsBuilder.getSteps());
        final Map<String, List<Integer>> multiplicitySections = multiplicitySectionsBuilder.buildMultiplicitySections(
                spectrum, dim);
        for (final Map.Entry<String, List<Integer>> entry : multiplicitySections.entrySet()) {
            for (final int section : entry.getValue()) {
                bitSetFingerprint.set(section, true);
            }
        }

        return bitSetFingerprint;
    }

    /**
     * Returns deviations between two already matched spectra.
     *
     * @param spectrum1   first spectrum
     * @param spectrum2   second spectrum
     * @param dim1        dimension in first spectrum to take the shifts from
     * @param dim2        dimension in second spectrum to take the shifts from
     * @param assignments assignments from previous matching
     *
     * @return
     */
    public static Double[] getDeviations(final Spectrum spectrum1, final Spectrum spectrum2, final int dim1,
                                         final int dim2, final Assignment assignments) {
        final Double[] deviations = new Double[spectrum1.getSignalCount()];
        Signal matchedSignalInSpectrum2;
        for (int i = 0; i
                < spectrum1.getSignalCount(); i++) {
            if (assignments.getAssignment(0, i).length
                    == 0) {
                deviations[i] = null;
            } else {
                matchedSignalInSpectrum2 = spectrum2.getSignal(assignments.getAssignment(0, i)[0]);
                deviations[i] = Math.abs(spectrum1.getSignal(i)
                                                  .getShift(dim1)
                                                 - matchedSignalInSpectrum2.getShift(dim2));
            }
        }

        return deviations;
    }

    /**
     * Returns deviations between matched shifts of two spectra.
     * The matching procedure is already included here.
     *
     * @param spectrum1                   first spectrum
     * @param spectrum2                   second spectrum
     * @param dim1                        dimension in first spectrum to take the shifts from
     * @param dim2                        dimension in second spectrum to take the shifts from
     * @param shiftTol                    shift tolerance
     * @param checkMultiplicity           indicates whether to compare the multiplicity of matched signals
     * @param checkEquivalencesCount      indicates whether to compare the equivalences counts of matched signals
     * @param allowLowerEquivalencesCount indicates to allow a lower equivalences counts spectrum 2
     *
     * @return
     *
     * @see #matchSpectra(Spectrum, Spectrum, int, int, double, boolean, boolean, boolean)
     */
    public static Double[] getDeviations(final Spectrum spectrum1, final Spectrum spectrum2, final int dim1,
                                         final int dim2, final double shiftTol, final boolean checkMultiplicity,
                                         final boolean checkEquivalencesCount,
                                         final boolean allowLowerEquivalencesCount) {
        final Assignment matchAssignments = matchSpectra(spectrum1, spectrum2, dim1, dim2, shiftTol, checkMultiplicity,
                                                         checkEquivalencesCount, allowLowerEquivalencesCount);

        return getDeviations(spectrum1, spectrum2, dim1, dim2, matchAssignments);
    }

    /**
     * Returns the average of all deviations of matched shifts between two
     * spectra.
     *
     * @param spectrum1                   first spectrum
     * @param spectrum2                   second spectrum
     * @param dim1                        dimension in first spectrum to take the shifts from
     * @param dim2                        dimension in second spectrum to take the shifts from
     * @param shiftTol                    Tolerance value [ppm] used during peak picking in
     *                                    shift comparison
     * @param checkMultiplicity           indicates whether to compare the multiplicity of matched signals
     * @param checkEquivalencesCount      indicates whether to compare the equivalences counts of matched signals
     * @param allowLowerEquivalencesCount indicates to allow a lower equivalences counts spectrum 2
     *
     * @return
     *
     * @see #getDeviations(Spectrum, Spectrum, int, int, double, boolean, boolean, boolean)
     * @see Statistics#calculateAverageDeviation(Double[])
     */
    public static Double calculateAverageDeviation(final Spectrum spectrum1, final Spectrum spectrum2, final int dim1,
                                                   final int dim2, final double shiftTol,
                                                   final boolean checkMultiplicity,
                                                   final boolean checkEquivalencesCount,
                                                   final boolean allowLowerEquivalencesCount) {
        return Statistics.calculateAverageDeviation(
                Similarity.getDeviations(spectrum1, spectrum2, dim1, dim2, shiftTol, checkMultiplicity,
                                         checkEquivalencesCount, allowLowerEquivalencesCount));
    }

    /**
     * Returns the average of all deviations of matched shifts between two
     * spectra.
     *
     * @param spectrum1   first spectrum
     * @param spectrum2   second spectrum
     * @param dim1        dimension in first spectrum to take the shifts from
     * @param dim2        dimension in second spectrum to take the shifts from
     * @param assignments assignments from previous matching
     *
     * @return
     *
     * @see #getDeviations(Spectrum, Spectrum, int, int, double, boolean, boolean, boolean)
     * @see Statistics#calculateAverageDeviation(Double[])
     */
    public static Double calculateAverageDeviation(final Spectrum spectrum1, final Spectrum spectrum2, final int dim1,
                                                   final int dim2, final Assignment assignments) {
        return Statistics.calculateAverageDeviation(
                Similarity.getDeviations(spectrum1, spectrum2, dim1, dim2, assignments));
    }

    /**
     * Returns the average of all deviations of matched shifts between two
     * spectra.
     *
     * @param spectrum1                   first spectrum
     * @param spectrum2                   second spectrum
     * @param dim1                        dimension in first spectrum to take the shifts from
     * @param dim2                        dimension in second spectrum to take the shifts from
     * @param shiftTol                    Tolerance value [ppm] used during peak picking in
     *                                    shift comparison
     * @param checkMultiplicity           indicates whether to compare the multiplicity of matched signals
     * @param checkEquivalencesCount      indicates whether to compare the equivalences counts of matched signals
     * @param allowLowerEquivalencesCount indicates to allow a lower equivalences counts spectrum 2
     *
     * @return
     *
     * @see #getDeviations(Spectrum, Spectrum, int, int, double, boolean, boolean, boolean)
     * @see Statistics#calculateAverageDeviation(Double[])
     */
    public static Double calculateRMSD(final Spectrum spectrum1, final Spectrum spectrum2, final int dim1,
                                       final int dim2, final double shiftTol, final boolean checkMultiplicity,
                                       final boolean checkEquivalencesCount,
                                       final boolean allowLowerEquivalencesCount) {
        return Statistics.calculateRMSD(
                Similarity.getDeviations(spectrum1, spectrum2, dim1, dim2, shiftTol, checkMultiplicity,
                                         checkEquivalencesCount, allowLowerEquivalencesCount));
    }

    /**
     * Returns the average of all deviations of matched shifts between two
     * spectra.
     *
     * @param spectrum1   first spectrum
     * @param spectrum2   second spectrum
     * @param dim1        dimension in first spectrum to take the shifts from
     * @param dim2        dimension in second spectrum to take the shifts from
     * @param assignments assignments from previous matching
     *
     * @return
     *
     * @see #getDeviations(Spectrum, Spectrum, int, int, double, boolean, boolean, boolean)
     * @see Statistics#calculateAverageDeviation(Double[])
     */
    public static Double calculateRMSD(final Spectrum spectrum1, final Spectrum spectrum2, final int dim1,
                                       final int dim2, final Assignment assignments) {
        return Statistics.calculateRMSD(Similarity.getDeviations(spectrum1, spectrum2, dim1, dim2, assignments));
    }

    /**
     * Returns the closest shift matches between two spectra in selected dimensions
     * as an Assignment object with one set dimension only. <br>
     *
     * @param spectrum1                   first spectrum (possible subspectrum)
     * @param spectrum2                   second spectrum
     * @param dim1                        dimension in first spectrum to take the shifts from
     * @param dim2                        dimension in second spectrum to take the shifts from
     * @param shiftTol                    Tolerance value [ppm] used during spectra shift
     *                                    comparison
     * @param checkMultiplicity           indicates whether to compare the multiplicity of matched signals
     * @param checkEquivalencesCount      indicates whether to compare the equivalences counts of matched signals
     * @param allowLowerEquivalencesCount indicates to allow a lower equivalences counts spectrum 2
     *
     * @return Assignments with signal indices of spectrum and matched indices
     * in query spectrum; null if one of the spectra does not
     * contain the selected dimension
     */
    public static Assignment matchSpectra(final Spectrum spectrum1, final Spectrum spectrum2, final int dim1,
                                          final int dim2, final double shiftTol, final boolean checkMultiplicity,
                                          final boolean checkEquivalencesCount,
                                          final boolean allowLowerEquivalencesCount) {
        if (!Similarity.checkDimensions(spectrum1, spectrum2, dim1, dim2)) {
            return null;
        }
        final Assignment matchAssignments = new Assignment();
        matchAssignments.setNuclei(new String[]{spectrum1.getNuclei()[dim1]});
        matchAssignments.initAssignments(spectrum1.getSignalCount());
        final Set<Integer> assigned = new HashSet<>();
        List<Integer> pickedSignalIndicesSpectrum2;
        boolean passed;

        for (int i = 0; i
                < spectrum1.getSignalCount(); i++) {
            if (spectrum1.getShift(i, dim1)
                    == null) {
                continue;
            }

            // @TODO add solvent deviation value for picking closest signal(s)
            pickedSignalIndicesSpectrum2 = new ArrayList<>();
            for (final int pickedSignalIndexSpectrum2 : spectrum2.pickSignals(spectrum1.getShift(i, dim1), dim2,
                                                                              shiftTol)) {
                passed = true;
                // @TODO maybe consider further parameters to check ? e.g. intensity
                if (checkMultiplicity) {
                    passed = (spectrum1.getMultiplicity(i)
                            == null
                            && spectrum2.getMultiplicity(pickedSignalIndexSpectrum2)
                            == null)
                            || (spectrum1.getMultiplicity(i)
                            != null
                            && spectrum1.getMultiplicity(i)
                                        .equals(spectrum2.getMultiplicity(pickedSignalIndexSpectrum2)));
                }
                if (passed
                        && checkEquivalencesCount) {
                    if (allowLowerEquivalencesCount) {
                        passed = spectrum1.getEquivalencesCount(i)
                                <= spectrum2.getEquivalencesCount(pickedSignalIndexSpectrum2);
                    } else {
                        passed = spectrum1.getEquivalencesCount(i)
                                == spectrum2.getEquivalencesCount(pickedSignalIndexSpectrum2);
                    }
                }

                if (passed) {
                    pickedSignalIndicesSpectrum2.add(pickedSignalIndexSpectrum2);
                }
            }
            for (final int pickedSignalIndexSpectrum2 : pickedSignalIndicesSpectrum2) {
                if (!assigned.contains(pickedSignalIndexSpectrum2)) {
                    // add signal to list of already assigned signals
                    assigned.add(pickedSignalIndexSpectrum2);
                    for (int k = 0; k
                            < spectrum1.getEquivalencesCount(i); k++) {
                        matchAssignments.addAssignmentEquivalence(0, i, pickedSignalIndexSpectrum2);
                    }
                    break;
                }
            }
        }

        return matchAssignments;
    }

    /**
     * Returns the closest shift matches between two spectra in all dimensions
     * as one Assignment object with N set dimensions.
     * N here means the number of dimensions in both spectra. <br>
     * Despite intensities are expected, they are still not considered here.
     *
     * @param spectrum1                   first spectrum (possible subspectrum)
     * @param spectrum2                   second spectrum
     * @param shiftTols                   tolerance values [ppm] per each dimension used during spectra shift
     *                                    comparisons
     * @param checkMultiplicity           indicates whether to compare the multiplicity of matched signals
     * @param checkEquivalencesCount      indicates whether to compare the equivalences counts of matched signals
     * @param allowLowerEquivalencesCount indicates to allow a lower equivalences counts spectrum 2
     *
     * @return Assignments with signal indices of spectrum1 and matched indices
     * in spectrum2 for each dimension; null if the number of
     * dimensions in both spectra is not the same or is different than the number of given
     * shift tolerances
     *
     * @see #matchSpectra(Spectrum, Spectrum, int, int, double, boolean, boolean, boolean)
     */
    public static Assignment matchSpectra(final Spectrum spectrum1, final Spectrum spectrum2, final double[] shiftTols,
                                          final boolean checkMultiplicity, final boolean checkEquivalencesCount,
                                          final boolean allowLowerEquivalencesCount) {
        if ((spectrum1.getNDim()
                != spectrum2.getNDim())
                || (spectrum1.getNDim()
                != shiftTols.length)) {
            return null;
        }
        final Assignment matchAssignment = new Assignment();
        matchAssignment.setNuclei(spectrum1.getNuclei());
        matchAssignment.initAssignments(spectrum1.getSignalCount());
        for (int dim = 0; dim
                < spectrum1.getNDim(); dim++) {
            matchAssignment.setAssignments(dim, matchSpectra(spectrum1, spectrum2, dim, dim, shiftTols[dim],
                                                             checkMultiplicity, checkEquivalencesCount,
                                                             allowLowerEquivalencesCount).getAssignments(0));
        }

        return matchAssignment;
    }

}
