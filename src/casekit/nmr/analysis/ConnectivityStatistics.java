package casekit.nmr.analysis;

import casekit.nmr.model.DataSet;
import casekit.nmr.model.Spectrum;
import casekit.nmr.utils.Utils;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

public class ConnectivityStatistics {

    /**
     * @param dataSetList
     * @param nucleus
     * @param connectivityStatistics multiplicity -> hybridization -> shift (int) -> connected atom symbol -> connected atom hybridization -> connected atom protons count -> occurrence
     */
    public static void buildConnectivityStatistics(final List<DataSet> dataSetList, final String nucleus,
                                                   final Map<String, Map<String, Map<Integer, Map<String, Map<String, Map<Integer, Integer>>>>>> connectivityStatistics) {
        final String atomType = Utils.getAtomTypeFromNucleus(nucleus);
        for (final DataSet dataSet : dataSetList) {
            if (!dataSet.getSpectrum()
                        .getNuclei()[0].equals(nucleus)) {
                continue;
            }
            buildConnectivityStatistics(dataSet, atomType, connectivityStatistics);
        }
    }

    /**
     * @param dataSet
     * @param atomType
     * @param connectivityStatistics multiplicity -> hybridization -> shift (int) -> connected atom symbol -> connected atom hybridization -> connected atom protons count -> occurrence
     */
    public static void buildConnectivityStatistics(final DataSet dataSet, final String atomType,
                                                   final Map<String, Map<String, Map<Integer, Map<String, Map<String, Map<Integer, Integer>>>>>> connectivityStatistics) {
        final IAtomContainer structure = dataSet.getStructure()
                                                .toAtomContainer();
        final Spectrum spectrum = dataSet.getSpectrum()
                                         .toSpectrum();
        int shift, atomIndex;
        IAtom atom;
        String multiplicity, hybridization, connectedAtomType, connectedAtomHybridization;
        for (int signalIndex = 0; signalIndex
                < spectrum.getSignalCount(); signalIndex++) {
            shift = spectrum.getShift(signalIndex, 0)
                            .intValue();
            for (int equivalenceIndex = 0; equivalenceIndex
                    < dataSet.getAssignment()
                             .getAssignment(0, signalIndex).length; equivalenceIndex++) {
                atomIndex = dataSet.getAssignment()
                                   .getAssignment(0, signalIndex, equivalenceIndex);
                atom = structure.getAtom(atomIndex);
                if (atom.getSymbol()
                        .equals(atomType)) {
                    multiplicity = Utils.getMultiplicityFromProtonsCount(atom.getImplicitHydrogenCount());
                    if (multiplicity
                            == null) {
                        continue;
                    }
                    multiplicity = multiplicity.toLowerCase();
                    hybridization = atom.getHybridization()
                                        .name();
                    connectivityStatistics.putIfAbsent(multiplicity, new ConcurrentHashMap<>());
                    connectivityStatistics.get(multiplicity)
                                          .putIfAbsent(hybridization, new ConcurrentHashMap<>());
                    // check for connected hetero atoms
                    for (final IAtom connectedAtom : structure.getConnectedAtomsList(atom)) {
                        if (connectedAtom.getSymbol()
                                         .equals("H")) {
                            continue;
                        }
                        connectedAtomType = connectedAtom.getSymbol();
                        if (connectedAtom.getHybridization()
                                == null) {
                            continue;
                        }
                        connectedAtomHybridization = connectedAtom.getHybridization()
                                                                  .name();
                        connectivityStatistics.get(multiplicity)
                                              .get(hybridization)
                                              .putIfAbsent(shift, new ConcurrentHashMap<>());
                        connectivityStatistics.get(multiplicity)
                                              .get(hybridization)
                                              .get(shift)
                                              .putIfAbsent(connectedAtomType, new ConcurrentHashMap<>());
                        connectivityStatistics.get(multiplicity)
                                              .get(hybridization)
                                              .get(shift)
                                              .get(connectedAtomType)
                                              .putIfAbsent(connectedAtomHybridization, new ConcurrentHashMap<>());
                        connectivityStatistics.get(multiplicity)
                                              .get(hybridization)
                                              .get(shift)
                                              .get(connectedAtomType)
                                              .get(connectedAtomHybridization)
                                              .putIfAbsent(connectedAtom.getImplicitHydrogenCount(), 0);
                        connectivityStatistics.get(multiplicity)
                                              .get(hybridization)
                                              .get(shift)
                                              .get(connectedAtomType)
                                              .get(connectedAtomHybridization)
                                              .put(connectedAtom.getImplicitHydrogenCount(), connectivityStatistics.get(
                                                                                                                           multiplicity)
                                                                                                                   .get(hybridization)
                                                                                                                   .get(shift)
                                                                                                                   .get(connectedAtomType)
                                                                                                                   .get(connectedAtomHybridization)
                                                                                                                   .get(connectedAtom.getImplicitHydrogenCount())
                                                      + 1);
                    }
                }
            }
        }
    }

    /**
     * @param connectivityStatistics   multiplicity -> hybridization -> shift (int) -> connected atom symbol -> connected atom hybridization -> connected atom protons count -> occurrence
     * @param multiplicity
     * @param hybridization
     * @param shift
     * @param molecularFormulaElements
     *
     * @return
     */
    public static Map<String, Map<String, Map<Integer, Integer>>> extractConnectivities(
            final Map<String, Map<String, Map<Integer, Map<String, Map<String, Map<Integer, Integer>>>>>> connectivityStatistics,
            final String multiplicity, final String hybridization, final int shift,
            final Set<String> molecularFormulaElements) {
        final Map<String, Map<String, Map<Integer, Integer>>> extractedConnectivities = new HashMap<>();
        if (connectivityStatistics.containsKey(multiplicity)
                && connectivityStatistics.get(multiplicity)
                                         .containsKey(hybridization)
                && connectivityStatistics.get(multiplicity)
                                         .get(hybridization)
                                         .containsKey(shift)) {
            for (final Map.Entry<String, Map<String, Map<Integer, Integer>>> entry : connectivityStatistics.get(
                                                                                                                   multiplicity)
                                                                                                           .get(hybridization)
                                                                                                           .get(shift)
                                                                                                           .entrySet()) {
                if (molecularFormulaElements.contains(entry.getKey())) {
                    extractedConnectivities.put(entry.getKey(), entry.getValue());
                }
            }
        }

        return extractedConnectivities;
    }

    public static Map<String, Set<Integer>> filterExtractedConnectivities(
            final Map<String, Map<String, Map<Integer, Integer>>> extractedConnectivities,
            final double thresholdElementCount) {
        final Map<String, Integer> totalCounts = getTotalCounts(extractedConnectivities);
        final int totalCountsSum = getTotalCount(totalCounts);

        final Map<String, Set<Integer>> filteredExtractedConnectivities = new HashMap<>();
        extractedConnectivities.keySet()
                               .forEach(neighborAtomType -> {
                                   extractedConnectivities.get(neighborAtomType)
                                                          .keySet()
                                                          .forEach(neighborHybridization -> {
                                                              for (final Map.Entry<Integer, Integer> entryProtonsCount : extractedConnectivities.get(
                                                                                                                                                        neighborAtomType)
                                                                                                                                                .get(neighborHybridization)
                                                                                                                                                .entrySet()) {
                                                                  if (entryProtonsCount.getValue()
                                                                          / (double) totalCountsSum
                                                                          >= thresholdElementCount) {
                                                                      filteredExtractedConnectivities.putIfAbsent(
                                                                              neighborAtomType, new HashSet<>());
                                                                      filteredExtractedConnectivities.get(
                                                                                                             neighborAtomType)
                                                                                                     .add(entryProtonsCount.getKey());
                                                                  }
                                                              }
                                                          });
                               });
        
        return filteredExtractedConnectivities;
    }

    private static Map<String, Integer> getTotalCounts(
            final Map<String, Map<String, Map<Integer, Integer>>> extractedConnectivities) {
        final Map<String, Integer> totalCounts = new HashMap<>();
        for (final String key1 : extractedConnectivities.keySet()) {
            totalCounts.putIfAbsent(key1, 0);
            for (final String key2 : extractedConnectivities.get(key1)
                                                            .keySet()) {
                for (final Map.Entry<Integer, Integer> countsEntry : extractedConnectivities.get(key1)
                                                                                            .get(key2)
                                                                                            .entrySet()) {
                    totalCounts.put(key1, totalCounts.get(key1)
                            + countsEntry.getValue());
                }
            }
        }

        return totalCounts;
    }

    private static int getTotalCount(final Map<String, Integer> totalCounts) {
        return totalCounts.values()
                          .stream()
                          .reduce(0, (total, current) -> total += current);
    }
}
