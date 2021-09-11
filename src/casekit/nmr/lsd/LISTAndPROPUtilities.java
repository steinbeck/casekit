package casekit.nmr.lsd;

import casekit.nmr.model.nmrium.Correlation;
import casekit.nmr.utils.Statistics;

import java.util.*;

public class LISTAndPROPUtilities {

    public static void insertELEM(final StringBuilder stringBuilder, final Map<String, String> listMap,
                                  final Set<String> atomTypesByMf) {
        final Set<String> atomTypes = new HashSet<>(atomTypesByMf);
        atomTypes.remove("H");
        for (final String atomType : atomTypes) {
            listMap.put(atomType, "L"
                    + (listMap.size()
                    + 1));
            stringBuilder.append("ELEM")
                         .append(" ")
                         .append(listMap.get(atomType))
                         .append(" ")
                         .append(atomType)
                         .append("\n");
        }
    }

    public static void insertNoHeteroHeteroBonds(final StringBuilder stringBuilder, final Map<String, String> listMap) {
        // create hetero atom list automatically to forbid hetero-hetero bonds
        stringBuilder.append("HETE L1")
                     .append("; hetero atoms\n");
        stringBuilder.append("PROP L1 0 L1 -; no hetero-hetero bonds\n");
        listMap.put("HETE", "L1");
    }

    public static void insertCarbonCombinationLISTs(final StringBuilder stringBuilder,
                                                    final Map<String, String> listMap,
                                                    final List<Correlation> correlationList,
                                                    final Map<Integer, Object[]> indicesMap,
                                                    final Map<Integer, Map<String, Map<String, Set<Integer>>>> detectedConnectivities) {
        final Map<String, Set<Integer>> atomIndicesMap = new HashMap<>();
        Correlation correlation;
        int indexInPyLSD;
        String listKey;
        for (int i = 0; i
                < correlationList.size(); i++) {
            for (int k = 1; k
                    < indicesMap.get(i).length; k++) {
                correlation = correlationList.get(i);
                if (!correlation.getAtomType()
                                .equals("C")
                        || correlation.getHybridization()
                                      .size()
                        != 1
                        || correlation.getProtonsCount()
                                      .size()
                        != 1) {
                    continue;
                }
                listKey = correlation.getAtomType()
                        + "_"
                        + correlation.getHybridization()
                                     .get(0)
                        + "_"
                        + correlation.getProtonsCount()
                                     .get(0);
                indexInPyLSD = (int) indicesMap.get(i)[k];
                atomIndicesMap.putIfAbsent(listKey, new HashSet<>());
                atomIndicesMap.get(listKey)
                              .add(indexInPyLSD);
            }
        }
        String[] split;
        for (final Map.Entry<String, Set<Integer>> combinationEntry : atomIndicesMap.entrySet()) {
            stringBuilder.append("LIST ")
                         .append("L")
                         .append(listMap.size()
                                         + 1);

            for (final Integer pyLSDAtomIndex : combinationEntry.getValue()) {
                stringBuilder.append(" ")
                             .append(pyLSDAtomIndex);
            }
            split = combinationEntry.getKey()
                                    .split("_");
            stringBuilder.append("; ")
                         .append(split[0])
                         .append("H")
                         .append(split[2])
                         .append(", ")
                         .append(split[1])
                         .append("\n");

            listMap.put(combinationEntry.getKey(), "L"
                    + (listMap.size()
                    + 1));
        }
    }

    public static void insertForbiddenConnectionLISTsAndPROPs(final StringBuilder stringBuilder,
                                                              final Map<String, String> listMap,
                                                              final List<Correlation> correlationList,
                                                              final Map<Integer, Object[]> indicesMap,
                                                              final Map<Integer, Map<String, Map<String, Set<Integer>>>> detectedConnectivities,
                                                              final Set<String> atomTypesByMf) {
        // insert ELEM for each heavy atom type in MF
        LISTAndPROPUtilities.insertELEM(stringBuilder, listMap, atomTypesByMf);
        // insert list combinations of carbon and hybridization states
        LISTAndPROPUtilities.insertCarbonCombinationLISTs(stringBuilder, listMap, correlationList, indicesMap,
                                                          detectedConnectivities);
        Correlation correlation;
        String atomType, listKey;
        int indexInPyLSD;
        Map<String, Map<String, Set<Integer>>> connectivities;
        Set<String> nonNeighborAtomTypes, neighborAtomTypes;
        Map<String, Map<Integer, Set<Integer>>> forbiddenNeighborHybridizationsAndProtonCounts;
        for (int i = 0; i
                < correlationList.size(); i++) {
            correlation = correlationList.get(i);
            atomType = correlation.getAtomType();
            connectivities = detectedConnectivities.get(i);
            // consider carbons here only, because of having complete connectivity information
            if (!atomType.equals("C")
                    || connectivities
                    == null
                    || connectivities.isEmpty()) {
                continue;
            }
            // define atom types of non-neighbors
            neighborAtomTypes = connectivities.keySet();
            nonNeighborAtomTypes = new HashSet<>(atomTypesByMf);
            nonNeighborAtomTypes.removeAll(neighborAtomTypes);
            nonNeighborAtomTypes.remove("H");
            forbiddenNeighborHybridizationsAndProtonCounts = Utilities.buildForbiddenNeighborHybridizationsAndProtonCounts(
                    connectivities, neighborAtomTypes);

            // put in the extracted information per correlation
            for (int k = 1; k
                    < indicesMap.get(i).length; k++) {
                indexInPyLSD = (int) indicesMap.get(i)[k];
                // forbid bonds to whole element groups
                for (final String nonNeighborAtomType : nonNeighborAtomTypes) {
                    stringBuilder.append("PROP ")
                                 .append(indexInPyLSD)
                                 .append(" 0 ")
                                 .append(listMap.get(nonNeighborAtomType))
                                 .append(" -")
                                 .append("; no bonds between ")
                                 .append(indexInPyLSD)
                                 .append(" (")
                                 .append(atomType)
                                 .append(", ")
                                 .append(Statistics.roundDouble(correlation.getSignal()
                                                                           .getDelta(), 2))
                                 .append(") and ")
                                 .append(listMap.get(nonNeighborAtomType))
                                 .append(" (")
                                 .append(nonNeighborAtomType)
                                 .append(")")
                                 .append("\n");
                }
                // forbid bonds to possible neighbors with certain hybridization states and proton counts
                for (final String neighborAtomType : forbiddenNeighborHybridizationsAndProtonCounts.keySet()) {
                    for (final int forbiddenNeighborHybridization : forbiddenNeighborHybridizationsAndProtonCounts.get(
                            neighborAtomType)
                                                                                                                  .keySet()) {
                        if (!neighborAtomType.equals("C")) {
                            continue;
                        }
                        for (final Integer forbiddenProtonsCount : forbiddenNeighborHybridizationsAndProtonCounts.get(
                                neighborAtomType)
                                                                                                                 .get(forbiddenNeighborHybridization)) {
                            listKey = neighborAtomType
                                    + "_SP"
                                    + forbiddenNeighborHybridization
                                    + "_"
                                    + forbiddenProtonsCount;
                            if (listMap.containsKey(listKey)) {
                                stringBuilder.append("PROP ")
                                             .append(indexInPyLSD)
                                             .append(" 0 ")
                                             .append(listMap.get(listKey))
                                             .append(" -")
                                             .append("; no bonds between ")
                                             .append(indexInPyLSD)
                                             .append(" (")
                                             .append(atomType)
                                             .append(", ")
                                             .append(Statistics.roundDouble(correlation.getSignal()
                                                                                       .getDelta(), 2))
                                             .append(") and ")
                                             .append(listMap.get(listKey))
                                             .append(" (")
                                             .append(neighborAtomType)
                                             .append(", SP")
                                             .append(forbiddenNeighborHybridization)
                                             .append(", ")
                                             .append(forbiddenProtonsCount)
                                             .append("H")
                                             .append(")")
                                             .append("\n");
                            }
                        }
                    }
                }
            }
        }
    }
}
