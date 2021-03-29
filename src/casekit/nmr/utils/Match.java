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
     * @see #matchSpectra(Spectrum, Spectrum, int, int, double)
     */
    public static Double[] getDeviations(final Spectrum spectrum1, final Spectrum spectrum2, final int dim1,
                                         final int dim2, final double shiftTol) {
        final Double[] deviations = new Double[spectrum1.getSignalCount()];
        final Assignment matchAssignments = Match.matchSpectra(spectrum1, spectrum2, dim1, dim2, shiftTol);
        Signal matchedSignalInSpectrum2;
        for (int i = 0; i
                < spectrum1.getSignalCount(); i++) {
            if (matchAssignments.getAssignment(0, i)
                    == -1) {
                deviations[i] = null;
            } else {
                matchedSignalInSpectrum2 = spectrum2.getSignal(matchAssignments.getAssignment(0, i));
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
     * @see #getDeviations(Spectrum, Spectrum, int, int, double)
     * @see #calculateAverageDeviation(Double[])
     */
    public static Double calculateAverageDeviation(final Spectrum spectrum1, final Spectrum spectrum2, final int dim1,
                                                   final int dim2, final double shiftTol) {
        return Match.calculateAverageDeviation(Match.getDeviations(spectrum1, spectrum2, dim1, dim2, shiftTol));
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
     * @see #getDeviations(Spectrum, Spectrum, int, int, double)
     * @see #calculateAverageDeviation(Double[])
     */
    public static Double calculateRMSD(final Spectrum spectrum1, final Spectrum spectrum2, final int dim1,
                                       final int dim2, final double shiftTol) {
        return Match.calculateRMSD(Match.getDeviations(spectrum1, spectrum2, dim1, dim2, shiftTol));
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
                                          final int dim2, final double shiftTol) {
        if (!Match.checkDimensions(spectrum, querySpectrum, dim1, dim2)) {
            return null;
        }
        final Assignment matchAssignments = new Assignment(spectrum);
        final Set<Integer> assigned = new HashSet<>();
        List<Integer> pickedSignalIndicesSpectrum2;

        for (int i = 0; i
                < spectrum.getSignalCount(); i++) {
            if (spectrum.getShift(i, dim1)
                    == null) {
                continue;
            }

            // @TODO add solvent deviation value for picking closest signal(s)
            pickedSignalIndicesSpectrum2 = new ArrayList<>();
            for (final int pickedSignalIndexSpectrum2 : querySpectrum.pickClosestSignal(spectrum.getShift(i, dim1),
                                                                                        dim2, shiftTol)) {
                // @TODO maybe consider further parameters to check ? e.g. intensity
                if (querySpectrum.getMultiplicity(pickedSignalIndexSpectrum2)
                                 .equals(spectrum.getMultiplicity(i))
                        && querySpectrum.getEquivalencesCount(pickedSignalIndexSpectrum2)
                        <= spectrum.getEquivalencesCount(i)) {
                    pickedSignalIndicesSpectrum2.add(pickedSignalIndexSpectrum2);
                }
            }
            for (final int pickedSignalIndexSpectrum2 : pickedSignalIndicesSpectrum2) {
                if (!assigned.contains(pickedSignalIndexSpectrum2)) {
                    // add signal to list of already assigned signals
                    assigned.add(pickedSignalIndexSpectrum2);
                    // set picked signal index in assignment object
                    matchAssignments.setAssignment(0, i, pickedSignalIndexSpectrum2);

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
     * @see #matchSpectra(Spectrum, Spectrum, int, int, double)
     */
    public static Assignment matchSpectra(final Spectrum spectrum1, final Spectrum spectrum2,
                                          final double[] shiftTols) {
        if ((spectrum1.getNDim()
                != spectrum2.getNDim())
                || (spectrum1.getNDim()
                != shiftTols.length)) {
            return null;
        }
        final Assignment matchAssignment = new Assignment(spectrum1);
        for (int dim = 0; dim
                < spectrum1.getNDim(); dim++) {
            matchAssignment.setAssignments(dim, Match.matchSpectra(spectrum1, spectrum2, dim, dim, shiftTols[dim])
                                                     .getAssignments(0));
        }

        return matchAssignment;
    }


    // might be useful in future to correct matches between spectra

    //    /**
    //     * Corrects a match list regarding a given shift list and an atom container.
    //     * This is useful when two ore more shift values (e.g. DEPT shifts) match
    //     * with the same atom in the atom container. So the purpose here is to
    //     * enable more unambiguous matches. This method first looks for unambiguous
    //     * matches and calculates the median of the difference values between the
    //     * shift list values and the shifts of atom container. Then, all shift list
    //     * values are adjusted (+/-) with this median value.
    //     *
    //     * @param shiftList1 Shift value list to search in
    //     * @param shiftList2 Shift value list to match in shiftList1
    //     * @param matchesInShiftList1 Match list to correct
    //     * @param tol Tolerance value
    //     * @return
    //     */
    //    public static ArrayList<Integer> correctShiftMatches(final ArrayList<Double> shiftList1, final ArrayList<Double> shiftList2, final ArrayList<Integer> matchesInShiftList1, final double tol) {
    //
    //        int matchIndex;
    //        // get differences of unique matches between query shift and ac shifts
    //        ArrayList<Double> diffs = new ArrayList<>();
    //        final HashSet<Integer> uniqueMatchIndicesSet = new HashSet<>(matchesInShiftList1);
    //        for (final int uniqueMatchIndex : uniqueMatchIndicesSet) {
    //            if (Collections.frequency(matchesInShiftList1, uniqueMatchIndex) == 1) {
    //                matchIndex = matchesInShiftList1.indexOf(uniqueMatchIndex);
    //                if (matchesInShiftList1.get(matchIndex) >= 0) {
    //                    diffs.add(shiftList2.get(matchIndex) - shiftList1.get(matchesInShiftList1.get(matchIndex)));
    //                }
    //            }
    //        }
    //        // calculate the median of found unique match differences
    //        if (diffs.size() > 0) {
    //            final double median = casekit.casekit.nmr.Utils.getMedian(diffs);
    //            // add or subtract the median of the differences to all shift list values (input) and match again then
    //            for (int i = 0; i < shiftList2.size(); i++) {
    //                shiftList2.set(i, shiftList2.get(i) - median);
    //            }
    //            // rematch
    //            return casekit.casekit.nmr.Utils.findShiftMatches(shiftList1, shiftList2, tol);
    //        }
    //
    //        return matchesInShiftList1;
    //    }
}