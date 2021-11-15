package casekit.nmr.lsd;

import casekit.nmr.model.Signal;
import casekit.nmr.model.nmrium.Correlation;
import casekit.nmr.utils.Statistics;
import casekit.nmr.utils.Utils;

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

    private static String buildListKey(final String atomType, final int hybridization, final int protonsCount) {
        return atomType
                + "_"
                + hybridization
                + "_"
                + protonsCount;
    }

    public static void insertHeavyAtomCombinationLISTs(final StringBuilder stringBuilder,
                                                       final Map<String, String> listMap,
                                                       final List<Correlation> correlationList,
                                                       final Map<Integer, Object[]> indicesMap) {
        final Map<String, Set<Integer>> atomIndicesMap = new HashMap<>();
        Correlation correlation;
        int indexInPyLSD;
        String listKey;
        for (int i = 0; i
                < correlationList.size(); i++) {
            for (int k = 1; k
                    < indicesMap.get(i).length; k++) {
                correlation = correlationList.get(i);
                if (correlation.getAtomType()
                               .equals("H")
                        || correlation.getHybridization()
                                      .size()
                        != 1
                        || correlation.getProtonsCount()
                                      .size()
                        != 1) {
                    continue;
                }
                listKey = buildListKey(correlation.getAtomType(), correlation.getHybridization()
                                                                             .get(0), correlation.getProtonsCount()
                                                                                                 .get(0));
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
                         .append(", SP")
                         .append(split[1])
                         .append("\n");

            listMap.put(combinationEntry.getKey(), "L"
                    + (listMap.size()
                    + 1));

            System.out.println("listMap: "
                                       + listMap);
        }
    }

    public static void insertForbiddenConnectionLISTsAndPROPs(final StringBuilder stringBuilder,
                                                              final Map<String, String> listMap,
                                                              final List<Correlation> correlationList,
                                                              final Map<Integer, Object[]> indicesMap,
                                                              final Map<Integer, Map<String, Map<Integer, Set<Integer>>>> detectedConnectivities,
                                                              final Map<Integer, Map<String, Map<Integer, Set<Integer>>>> forbiddenNeighbors) {
        Correlation correlation;
        Signal signal;
        String atomType, listKey;
        int indexInPyLSD;
        Map<String, Map<Integer, Set<Integer>>> connectivitiesTemp;
        Map<String, Map<Integer, Set<Integer>>> forbiddenNeighborsTemp;
        for (int i = 0; i
                < correlationList.size(); i++) {
            correlation = correlationList.get(i);
            signal = Utils.extractSignalFromCorrelation(correlation);
            atomType = correlation.getAtomType();
            connectivitiesTemp = detectedConnectivities.get(i);
            // consider carbons here only, because of having complete connectivity information
            if (!atomType.equals("C")
                    || connectivitiesTemp
                    == null
                    || connectivitiesTemp.isEmpty()
                    || signal
                    == null) {
                continue;
            }
            forbiddenNeighborsTemp = forbiddenNeighbors.get(i);

            // put in the extracted information per correlation
            for (int k = 1; k
                    < indicesMap.get(i).length; k++) {
                indexInPyLSD = (int) indicesMap.get(i)[k];
                for (final String neighborAtomType : forbiddenNeighborsTemp.keySet()) {
                    // forbid bonds to whole element groups if there is an empty map for an atom type
                    if (forbiddenNeighborsTemp.get(neighborAtomType)
                                              .isEmpty()) {
                        stringBuilder.append("PROP ")
                                     .append(indexInPyLSD)
                                     .append(" 0 ")
                                     .append(listMap.get(neighborAtomType))
                                     .append(" -")
                                     .append("; no bonds between ")
                                     .append(indexInPyLSD)
                                     .append(" (")
                                     .append(atomType)
                                     .append(", ")
                                     .append(Statistics.roundDouble(signal.getShift(0), 2))
                                     .append(") and ")
                                     .append(listMap.get(neighborAtomType))
                                     .append(" (")
                                     .append(neighborAtomType)
                                     .append(")")
                                     .append("\n");
                    } else {
                        // forbid bonds to possible neighbors with certain hybridization states and proton counts
                        for (final int forbiddenNeighborHybridization : forbiddenNeighborsTemp.get(neighborAtomType)
                                                                                              .keySet()) {
                            for (final int forbiddenProtonsCount : forbiddenNeighborsTemp.get(neighborAtomType)
                                                                                         .get(forbiddenNeighborHybridization)) {
                                listKey = buildListKey(neighborAtomType, forbiddenNeighborHybridization,
                                                       forbiddenProtonsCount);
                                System.out.println("forbidden: "
                                                           + listKey);
                                if (listMap.containsKey(listKey)) {
                                    System.out.println("-> HUHU");
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
                                                 .append(Statistics.roundDouble(signal.getShift(0), 2))
                                                 .append(", SP")
                                                 .append(forbiddenNeighborHybridization)
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

    public static void insertSetConnectionLISTsAndPROPs(final StringBuilder stringBuilder,
                                                        final Map<String, String> listMap,
                                                        final List<Correlation> correlationList,
                                                        final Map<Integer, Object[]> indicesMap,
                                                        final Map<Integer, Map<String, Map<Integer, Set<Integer>>>> setNeighbors) {
        Correlation correlation;
        Signal signal;
        String atomType, listKey;
        int indexInPyLSD;
        Map<String, Map<Integer, Set<Integer>>> setNeighborsTemp;
        for (int i = 0; i
                < correlationList.size(); i++) {
            if (setNeighbors.containsKey(i)) {
                correlation = correlationList.get(i);
                signal = Utils.extractSignalFromCorrelation(correlation);
                atomType = correlation.getAtomType();
                setNeighborsTemp = setNeighbors.get(i);

                // put in the extracted information per correlation
                for (int k = 1; k
                        < indicesMap.get(i).length; k++) {
                    indexInPyLSD = (int) indicesMap.get(i)[k];
                    for (final String neighborAtomType : setNeighborsTemp.keySet()) {
                        // forbid bonds to whole element groups if there is an empty map for an atom type
                        if (setNeighborsTemp.get(neighborAtomType)
                                            .isEmpty()) {
                            stringBuilder.append("PROP ")
                                         .append(indexInPyLSD)
                                         .append(" 1 ")
                                         .append(listMap.get(neighborAtomType))
                                         .append(" +")
                                         .append("; at least one bond between ")
                                         .append(indexInPyLSD)
                                         .append(" (")
                                         .append(atomType)
                                         .append(", ")
                                         .append(signal
                                                         != null
                                                 ? Statistics.roundDouble(signal.getShift(0), 2)
                                                 : "?")
                                         .append(") and ")
                                         .append(listMap.get(neighborAtomType))
                                         .append(" (")
                                         .append(neighborAtomType)
                                         .append(")")
                                         .append("\n");
                        } else {
                            for (final int setNeighborHybridization : setNeighborsTemp.get(neighborAtomType)
                                                                                      .keySet()) {
                                for (final int setProtonsCount : setNeighborsTemp.get(neighborAtomType)
                                                                                 .get(setNeighborHybridization)) {
                                    listKey = buildListKey(neighborAtomType, setNeighborHybridization, setProtonsCount);
                                    System.out.println("set: "
                                                               + listKey);
                                    if (listMap.containsKey(listKey)) {
                                        System.out.println("-> HUHU");
                                        stringBuilder.append("PROP ")
                                                     .append(indexInPyLSD)
                                                     .append(" 1 ")
                                                     .append(listMap.get(listKey))
                                                     .append(" +")
                                                     .append("; at least one bond between ")
                                                     .append(indexInPyLSD)
                                                     .append(" (")
                                                     .append(atomType)
                                                     .append(", ")
                                                     .append(signal
                                                                     != null
                                                             ? Statistics.roundDouble(signal.getShift(0), 2)
                                                             : "?")
                                                     .append(") and ")
                                                     .append(listMap.get(listKey))
                                                     .append(" (")
                                                     .append(neighborAtomType)
                                                     .append(", SP")
                                                     .append(setNeighborHybridization)
                                                     .append(", ")
                                                     .append(setProtonsCount)
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
    }
}
