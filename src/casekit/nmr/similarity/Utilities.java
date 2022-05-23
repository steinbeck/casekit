package casekit.nmr.similarity;

import casekit.nmr.elucidation.Constants;
import casekit.nmr.elucidation.model.Detections;
import casekit.nmr.model.Assignment;
import casekit.nmr.model.Signal;
import casekit.nmr.model.Spectrum;
import casekit.nmr.similarity.model.Distance;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.*;

public class Utilities {


    /**
     * @param spectrum1                   first spectrum (possible subspectrum)
     * @param spectrum2                   second spectrum
     * @param dim1                        dim in first spectrum
     * @param dim2                        dim in second spectrum
     * @param shiftTolerance              shift tolerance
     * @param checkMultiplicity           whether to check multiplicity
     * @param checkEquivalencesCount      whether to check equivalences
     * @param allowLowerEquivalencesCount whether to allow lower equivalences
     *
     * @return
     */
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
                              .equalsIgnoreCase(signal2.getMultiplicity()));
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

    /**
     * @param spectrum1                   first spectrum (possible subspectrum)
     * @param spectrum2                   second spectrum
     * @param dim1                        dim in first spectrum
     * @param dim2                        dim in second spectrum
     * @param shiftTolerance              shift tolerance
     * @param checkMultiplicity           whether to check multiplicity
     * @param checkEquivalencesCount      whether to check equivalences
     * @param allowLowerEquivalencesCount whether to allow lower equivalences
     * @param structure                   structure belonging to first spectrum
     * @param assignment                  assignments between structure and first spectrum
     * @param detections                  detections to use as structural filter within given structure
     *
     * @return
     */
    public static List<Distance> buildDistanceList(final Spectrum spectrum1, final Spectrum spectrum2, final int dim1,
                                                   final int dim2, final double shiftTolerance,
                                                   final boolean checkMultiplicity,
                                                   final boolean checkEquivalencesCount,
                                                   final boolean allowLowerEquivalencesCount,
                                                   final IAtomContainer structure, final Assignment assignment,
                                                   final Detections detections) {
        final List<Distance> distanceList = new ArrayList<>();
        Double distanceValue;
        List<Integer> hybridizations;
        Set<String> forbiddenNeighbors, setNeighbors, setNeighborsTemp;
        IAtom atom;
        boolean skip;
        for (int i = 0; i
                < spectrum1.getSignalCount(); i++) {
            for (int j = 0; j
                    < spectrum2.getSignalCount(); j++) {
                // check spectral constraints
                distanceValue = getDistanceValue(spectrum1.getSignal(i), spectrum2.getSignal(j), dim1, dim2,
                                                 checkMultiplicity, checkEquivalencesCount, allowLowerEquivalencesCount,
                                                 shiftTolerance);
                if (distanceValue
                        == null) {
                    continue;
                }
                skip = false;
                // check structural constraints
                forbiddenNeighbors = detections.getForbiddenNeighbors()
                                               .get(j)
                                               .keySet();
                setNeighbors = detections.getSetNeighbors()
                                         .get(j)
                                         .keySet();
                hybridizations = detections.getDetectedHybridizations()
                                           .get(j);
                for (int equiv = 0; equiv
                        < assignment.getAssignment(0, i).length; equiv++) {
                    atom = structure.getAtom(assignment.getAssignment(0, i, equiv));
                    // if certain hybridizations are given and the atom's hybridization is known
                    if (!hybridizations.isEmpty()
                            && Constants.hybridizationConversionMap.containsKey(atom.getHybridization()
                                                                                    .name())) {
                        skip = !hybridizations.contains(Constants.hybridizationConversionMap.get(atom.getHybridization()
                                                                                                     .name()));
                        if (skip) {
                            break;
                        }
                    }
                    setNeighborsTemp = new HashSet<>(setNeighbors);
                    for (final IAtom neighborAtom : structure.getConnectedAtomsList(atom)) {
                        skip = forbiddenNeighbors.contains(neighborAtom.getSymbol());
                        if (skip) {
                            break;
                        }
                        setNeighborsTemp.remove(neighborAtom.getSymbol());
                    }
                    if (!setNeighborsTemp.isEmpty()) {
                        skip = true;
                    }
                    if (skip) {
                        break;
                    }
                }
                if (skip) {
                    continue;
                }

                // if passed all check then add to distance list
                distanceList.add(new Distance(i, j, distanceValue));
            }
        }
        distanceList.sort(Comparator.comparingDouble(Distance::getValue));

        return distanceList;
    }
}
