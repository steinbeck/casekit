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

package casekit.nmr.utils;

import casekit.nmr.Utils;
import casekit.nmr.model.Assignment;
import casekit.nmr.model.Signal;
import casekit.nmr.model.Spectrum;
import org.apache.commons.lang3.ArrayUtils;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.similarity.Tanimoto;

import java.util.*;

public class Match {


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
     * Combines selected dimensions of two spectra while considering possible equivalent signals
     * via the {@code pickPrecision} parameter and multiplicity comparison.
     * In {@code spectrum1}, the equivalent signals have to be set.
     *
     * @param spectrum1     first spectrum, incl. equivalent signals
     * @param spectrum2     second spectrum
     * @param pickPrecision tolerance value used for signal shift matching to
     *                      find equivalent signals
     * @param dim1          dimension of first spectrum to combine
     * @param dim2          dimension of second spectrum to combine
     *
     * @return null if one spectrum does not contain the selected dimension
     */
    public static Spectrum combineSpectra(final Spectrum spectrum1, final Spectrum spectrum2, final int dim1,
                                          final int dim2, final double pickPrecision) throws Exception {
        if (!Match.checkDimensions(spectrum1, spectrum2, dim1, dim2)) {
            return null;
        }
        // create new spectra which is to fill with signals of both spectra
        final Spectrum combinedSpectrum = spectrum1.buildClone();
        // fill in signals from spectrum2
        // consider the possibility of potential equivalent signals here
        for (final Signal signalSpectrum2 : spectrum2.getSignals()) {
            combinedSpectrum.addSignal(signalSpectrum2.buildClone(), pickPrecision);
        }
        return combinedSpectrum;
    }

    /**
     * Calculates the Tanimoto coefficient between two spectra in given dimensions.
     *
     * @param spectrum1 first spectrum
     * @param spectrum2 second spectrum
     * @param dim1      dimension in first spectrum to take the shifts from
     * @param dim2      dimension in second spectrum to take the shifts from
     *
     * @return
     *
     * @throws CDKException
     */
    public static Float calculateTanimotoCoefficient(final Spectrum spectrum1, final Spectrum spectrum2, final int dim1,
                                                     final int dim2) throws CDKException {
        if (!Match.checkDimensions(spectrum1, spectrum2, dim1, dim2)) {
            return null;
        }
        final double[] shiftsSpectrum1 = ArrayUtils.toPrimitive(spectrum1.getShifts(dim1)
                                                                         .toArray(
                                                                                 new Double[spectrum1.getSignalCount()]));
        Arrays.parallelSort(shiftsSpectrum1);
        final double[] shiftsSpectrum2 = ArrayUtils.toPrimitive(spectrum2.getShifts(dim2)
                                                                         .toArray(
                                                                                 new Double[spectrum2.getSignalCount()]));
        Arrays.parallelSort(shiftsSpectrum2);

        return Tanimoto.calculate(shiftsSpectrum1, shiftsSpectrum2);
    }

    /**
     * Returns deviations between matched shifts of two spectra.
     * The matching procedure is already included here.
     *
     * @param spectrum1 first spectrum
     * @param spectrum2 second spectrum
     * @param dim1      dimension in first spectrum to take the shifts from
     * @param dim2      dimension in second spectrum to take the shifts from
     * @param shiftTol
     *
     * @return
     *
     * @see #matchSpectra(Spectrum, Spectrum, int, int, double, boolean, boolean)
     */
    public static Double[] getDeviations(final Spectrum spectrum1, final Spectrum spectrum2, final int dim1,
                                         final int dim2, final double shiftTol, final boolean checkMultiplicity,
                                         final boolean checkEquivalencesCount) {
        final Double[] deviations = new Double[spectrum1.getSignalCount()];
        final Assignment matchAssignments = matchSpectra(spectrum1, spectrum2, dim1, dim2, shiftTol, checkMultiplicity,
                                                         checkEquivalencesCount);
        Signal matchedSignalInSpectrum2;
        for (int i = 0; i
                < spectrum1.getSignalCount(); i++) {
            if (matchAssignments.getAssignment(0, i).length
                    == 0) {
                deviations[i] = null;
            } else {
                matchedSignalInSpectrum2 = spectrum2.getSignal(matchAssignments.getAssignment(0, i)[0]);
                deviations[i] = Math.abs(spectrum1.getSignal(i)
                                                  .getShift(dim1)
                                                 - matchedSignalInSpectrum2.getShift(dim2));
            }
        }
        return deviations;
    }

    /**
     * Returns the average of all deviations within a given input array.
     *
     * @param deviations array of deviations
     *
     * @return
     */
    public static Double calculateAverageDeviation(final Double[] deviations) {
        // every signal has to have a match
        for (final Double deviation : deviations) {
            if (deviation
                    == null) {
                return null;
            }
        }

        return Utils.getMean(deviations);
    }

    /**
     * Returns the average of all deviations of matched shifts between two
     * spectra.
     *
     * @param spectrum1 first spectrum
     * @param spectrum2 second spectrum
     * @param dim1      dimension in first spectrum to take the shifts from
     * @param dim2      dimension in second spectrum to take the shifts from
     * @param shiftTol  Tolerance value [ppm] used during peak picking in
     *                  shift comparison
     *
     * @return
     *
     * @see #getDeviations(Spectrum, Spectrum, int, int, double, boolean, boolean)
     * @see #calculateAverageDeviation(Double[])
     */
    public static Double calculateAverageDeviation(final Spectrum spectrum1, final Spectrum spectrum2, final int dim1,
                                                   final int dim2, final double shiftTol,
                                                   final boolean checkMultiplicity,
                                                   final boolean checkEquivalencesCount) {
        return Match.calculateAverageDeviation(
                Match.getDeviations(spectrum1, spectrum2, dim1, dim2, shiftTol, checkMultiplicity,
                                    checkEquivalencesCount));
    }

    /**
     * Returns the average of all deviations within a given input array.
     *
     * @param data array of deviations
     *
     * @return
     */
    public static Double calculateRMSD(final Double[] data) {
        // every signal has to have a match
        for (final Double value : data) {
            if (value
                    == null) {
                return null;
            }
        }

        return casekit.nmr.utils.Utils.getRMSD(data);
    }

    /**
     * Returns the average of all deviations of matched shifts between two
     * spectra.
     *
     * @param spectrum1 first spectrum
     * @param spectrum2 second spectrum
     * @param dim1      dimension in first spectrum to take the shifts from
     * @param dim2      dimension in second spectrum to take the shifts from
     * @param shiftTol  Tolerance value [ppm] used during peak picking in
     *                  shift comparison
     *
     * @return
     *
     * @see #getDeviations(Spectrum, Spectrum, int, int, double, boolean, boolean)
     * @see #calculateAverageDeviation(Double[])
     */
    public static Double calculateRMSD(final Spectrum spectrum1, final Spectrum spectrum2, final int dim1,
                                       final int dim2, final double shiftTol, final boolean checkMultiplicity,
                                       final boolean checkEquivalencesCount) {
        return Match.calculateRMSD(Match.getDeviations(spectrum1, spectrum2, dim1, dim2, shiftTol, checkMultiplicity,
                                                       checkEquivalencesCount));
    }

    /**
     * Returns the closest shift matches between two spectra in selected dimensions
     * as an Assignment object with one set dimension only. <br>
     * Despite intensities are expected, they are still not considered here.
     *
     * @param spectrum      first spectrum
     * @param querySpectrum query spectrum (Subspectrum)
     * @param dim1          dimension in first spectrum to take the shifts from
     * @param dim2          dimension in second spectrum to take the shifts from
     * @param shiftTol      Tolerance value [ppm] used during spectra shift
     *                      comparison
     *
     * @return Assignments with signal indices of spectrum and matched indices
     * in query spectrum; null if one of the spectra does not
     * contain the selected dimension
     */
    public static Assignment matchSpectra(final Spectrum spectrum, final Spectrum querySpectrum, final int dim1,
                                          final int dim2, final double shiftTol, final boolean checkMultiplicity,
                                          final boolean checkEquivalencesCount) {
        if (!Match.checkDimensions(spectrum, querySpectrum, dim1, dim2)) {
            return null;
        }
        final Assignment matchAssignments = new Assignment();
        matchAssignments.setNuclei(new String[]{spectrum.getNuclei()[dim1]});
        matchAssignments.initAssignments(spectrum.getSignalCount());
        final Set<Integer> assigned = new HashSet<>();
        List<Integer> pickedSignalIndicesSpectrum2;
        boolean passed;

        for (int i = 0; i
                < spectrum.getSignalCount(); i++) {
            if (spectrum.getShift(i, dim1)
                    == null) {
                continue;
            }

            // @TODO add solvent deviation value for picking closest signal(s)
            pickedSignalIndicesSpectrum2 = new ArrayList<>();
            for (final int pickedSignalIndexSpectrum2 : querySpectrum.pickSignals(spectrum.getShift(i, dim1), dim2,
                                                                                  shiftTol)) {
                passed = true;
                // @TODO maybe consider further parameters to check ? e.g. intensity
                if (checkMultiplicity) {
                    passed = querySpectrum.getMultiplicity(pickedSignalIndexSpectrum2)
                                          .equals(spectrum.getMultiplicity(i));
                }
                if (checkEquivalencesCount) {
                    passed = querySpectrum.getEquivalencesCount(pickedSignalIndexSpectrum2)
                            <= spectrum.getEquivalencesCount(i);
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
                            < spectrum.getEquivalencesCount(i); k++) {
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
     * @param spectrum1 first spectrum
     * @param spectrum2 second spectrum (query)
     * @param shiftTols tolerance values [ppm] per each dimension used during spectra shift
     *                  comparisons
     *
     * @return Assignments with signal indices of spectrum1 and matched indices
     * in spectrum2 for each dimension; null if the number of
     * dimensions in both spectra is not the same or is different than the number of given
     * shift tolerances
     *
     * @see #matchSpectra(Spectrum, Spectrum, int, int, double, boolean, boolean)
     */
    public static Assignment matchSpectra(final Spectrum spectrum1, final Spectrum spectrum2, final double[] shiftTols,
                                          final boolean checkMultiplicity, final boolean checkEquivalencesCount) {
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
                                                             checkMultiplicity, checkEquivalencesCount).getAssignments(
                    0));
        }

        return matchAssignment;
    }
}
