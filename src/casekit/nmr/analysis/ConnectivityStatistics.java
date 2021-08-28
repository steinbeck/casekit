package casekit.nmr.analysis;

import casekit.nmr.model.DataSet;
import casekit.nmr.model.Spectrum;
import casekit.nmr.utils.Utils;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

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

    public static List<String> extractNeighborAtomTypes(
            final Map<String, Map<String, Map<Integer, Integer>>> extractedConnectivities, final double threshold) {
        final Map<String, Integer> totalCounts = new HashMap<>();
        for (final String connectedAtomType : extractedConnectivities.keySet()) {
            totalCounts.putIfAbsent(connectedAtomType, 0);
            for (final String connectedAtomHybridization : extractedConnectivities.get(connectedAtomType)
                                                                                  .keySet()) {
                for (final Map.Entry<Integer, Integer> countsEntry : extractedConnectivities.get(connectedAtomType)
                                                                                            .get(connectedAtomHybridization)
                                                                                            .entrySet()) {
                    totalCounts.put(connectedAtomType, totalCounts.get(connectedAtomType)
                            + countsEntry.getValue());
                }
            }
        }
        final int totalCountsSum = totalCounts.values()
                                              .stream()
                                              .reduce(0, (total, current) -> total += current);
        return totalCounts.keySet()
                          .stream()
                          .filter(atomType -> (totalCounts.get(atomType)
                                  / (double) totalCountsSum)
                                  >= threshold)
                          .collect(Collectors.toList());
    }

    public static Map<String, Map<Integer, Integer>> extractNeighborHybridizations(
            final Map<String, Map<String, Map<Integer, Integer>>> extractedConnectivities,
            final String neighborAtomType, final double threshold) {
        if (!extractedConnectivities.containsKey(neighborAtomType)) {
            return new HashMap<>();
        }
        final Map<String, Integer> countsPerHybridization = new HashMap<>();
        extractedConnectivities.get(neighborAtomType)
                               .keySet()
                               .forEach(hybridizationNeighbor -> {
                                   countsPerHybridization.put(hybridizationNeighbor, extractedConnectivities.get(
                                           neighborAtomType)
                                                                                                            .get(hybridizationNeighbor)
                                                                                                            .keySet()
                                                                                                            .stream()
                                                                                                            .reduce(0,
                                                                                                                    (protonsCountSum, protonsCount) -> protonsCountSum += extractedConnectivities.get(
                                                                                                                            neighborAtomType)
                                                                                                                                                                                                 .get(hybridizationNeighbor)
                                                                                                                                                                                                 .get(protonsCount)));
                               });
        final int totalCount = countsPerHybridization.keySet()
                                                     .stream()
                                                     .map(countsPerHybridization::get)

                                                     .reduce(0, (sum, current) -> sum += current);
        final List<String> allowedHybridizationList = countsPerHybridization.keySet()
                                                                            .stream()
                                                                            .filter(hybridizationNeighbor -> countsPerHybridization.get(
                                                                                    hybridizationNeighbor)
                                                                                    / (double) totalCount
                                                                                    >= threshold)
                                                                            .collect(Collectors.toList());

        final Map<String, Map<Integer, Integer>> extractedNeighborHybridizationMap = new HashMap<>();
        for (final String allowedHybridization : allowedHybridizationList) {
            extractedNeighborHybridizationMap.put(allowedHybridization, extractedConnectivities.get(neighborAtomType)
                                                                                               .get(allowedHybridization));
        }

        return extractedNeighborHybridizationMap;
    }
}
