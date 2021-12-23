package casekit.nmr.lsd;

import casekit.nmr.model.Signal;
import casekit.nmr.model.nmrium.Correlation;
import casekit.nmr.utils.Statistics;
import casekit.nmr.utils.Utils;

import java.util.*;

public class LISTAndPROPUtilities {

    public static void insertELEM(final StringBuilder stringBuilder, final Map<String, Object[]> listMap,
                                  final Set<String> atomTypesByMf) {
        final Set<String> atomTypes = new HashSet<>(atomTypesByMf);
        atomTypes.remove("H");
        for (final String atomType : atomTypes) {
            listMap.put(atomType, new Object[]{"L"
                                                       + (listMap.size()
                    + 1)});
            stringBuilder.append("ELEM")
                         .append(" ")
                         .append(listMap.get(atomType)[0])
                         .append(" ")
                         .append(atomType)
                         .append("\n");
        }
    }

    public static void insertNoHeteroHeteroBonds(final StringBuilder stringBuilder,
                                                 final Map<String, Object[]> listMap) {
        // create hetero atom list automatically to forbid hetero-hetero bonds
        listMap.put("HETE", new Object[]{"L"
                                                 + (listMap.size()
                + 1)});
        stringBuilder.append("HETE ")
                     .append((String) listMap.get("HETE")[0])
                     .append("; hetero atoms\n");
        stringBuilder.append("PROP L1 0 L1 -; no hetero-hetero bonds\n");
    }

    private static String buildListKey(final String atomType, final List<Integer> hybridizations,
                                       final List<Integer> protonsCounts) {
        return atomType
                + "_"
                + (!hybridizations.isEmpty()
                   ? hybridizations
                   : "*")
                + "_"
                + (!protonsCounts.isEmpty()
                   ? protonsCounts
                   : "*");
    }

    public static void insertHeavyAtomCombinationLISTs(final StringBuilder stringBuilder,
                                                       final Map<String, Object[]> listMap,
                                                       final List<Correlation> correlationList,
                                                       final Map<Integer, Object[]> indicesMap) {
        final Map<String, Set<Integer>> atomIndicesMap = new LinkedHashMap<>();
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
                listKey = buildListKey(correlation.getAtomType(), new ArrayList<>(), //correlation.getHybridization(),
                                       correlation.getProtonsCount());
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

            for (final int pyLSDAtomIndex : combinationEntry.getValue()) {
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

            listMap.put(combinationEntry.getKey(), new Object[]{"L"
                                                                        + (listMap.size()
                    + 1), combinationEntry.getValue().size()});
        }
    }


    private static boolean checkSkipPROPInsertion(final Map<String, Object[]> listMap,
                                                  final Map<String, Integer> usedPropsCount, final String listKey) {
        if (!listMap.containsKey(listKey)) {
            return true;
        }
        // LSD crashes if we try to use more atoms with specific hybridization or/and proton count
        // thus count and limit that usage as following:
        // list name -> number of usages
        usedPropsCount.putIfAbsent((String) listMap.get(listKey)[0], 1);
        return listMap.get(listKey).length
                > 1
                && usedPropsCount.get((String) listMap.get(listKey)[0])
                > (int) listMap.get(listKey)[1];
    }

    private static void insertPROP(final StringBuilder stringBuilder, final Map<String, Object[]> listMap,
                                   final String atomType, final Signal signal, final String neighborAtomType,
                                   final int indexInPyLSD, final String listKey, final String mode) {
        stringBuilder.append("PROP ")
                     .append(indexInPyLSD)
                     .append(mode.equals("forbid")
                             ? " 0 "
                             : " 1 ")
                     .append(listMap.get(listKey)[0])
                     .append(mode.equals("forbid")
                             ? " -"
                             : " +")
                     .append(mode.equals("forbid")
                             ? "; no bonds between "
                             : "; at least one bond between ")
                     .append(indexInPyLSD)
                     .append(" (")
                     .append(atomType)
                     .append(", ")
                     .append(signal
                                     != null
                             ? Statistics.roundDouble(signal.getShift(0), 2)
                             : "?")
                     .append(") and ")
                     .append(Arrays.toString(listMap.get(listKey)))
                     .append("\n");
    }

    public static void insertConnectionLISTsAndPROPs(final StringBuilder stringBuilder,
                                                     final Map<String, Object[]> listMap,
                                                     final List<Correlation> correlationList,
                                                     final Map<Integer, Object[]> indicesMap,
                                                     final Map<Integer, Map<String, Map<Integer, Set<Integer>>>> neighbors,
                                                     final String mode) {
        Correlation correlation;
        Signal signal;
        String atomType, listKey;
        int indexInPyLSD;
        Map<String, Map<Integer, Set<Integer>>> neighborsTemp;
        final Map<String, Integer> usedPropsCount = new HashMap<>();
        for (int i = 0; i
                < correlationList.size(); i++) {
            if (neighbors.containsKey(i)) {
                correlation = correlationList.get(i);
                signal = Utils.extractSignalFromCorrelation(correlation);
                atomType = correlation.getAtomType();
                neighborsTemp = neighbors.get(i);

                // put in the extracted information per correlation and equivalent
                for (int k = 1; k
                        < indicesMap.get(i).length; k++) {
                    indexInPyLSD = (int) indicesMap.get(i)[k];
                    for (final String neighborAtomType : neighborsTemp.keySet()) {
                        // forbid/set bonds to whole element groups if there is an empty map for an atom type
                        if (neighborsTemp.get(neighborAtomType)
                                         .isEmpty()) {
                            insertPROP(stringBuilder, listMap, atomType, signal, neighborAtomType, indexInPyLSD,
                                       neighborAtomType, mode);
                        } else {
                            for (final int neighborHybridization : neighborsTemp.get(neighborAtomType)
                                                                                .keySet()) {
                                for (final int protonsCount : neighborsTemp.get(neighborAtomType)
                                                                           .get(neighborHybridization)) {
                                    listKey = buildListKey(neighborAtomType, neighborHybridization
                                                                                     == -1
                                                                             ? new ArrayList<>()
                                                                             : List.of(neighborHybridization),
                                                           List.of(protonsCount));
                                    if (checkSkipPROPInsertion(listMap, usedPropsCount, listKey)) {
                                        continue;
                                    }
                                    if (listMap.containsKey(listKey)) {
                                        insertPROP(stringBuilder, listMap, atomType, signal, neighborAtomType,
                                                   indexInPyLSD, listKey, mode);
                                        usedPropsCount.put((String) listMap.get(listKey)[0],
                                                           usedPropsCount.get((String) listMap.get(listKey)[0])
                                                                   + 1);
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
