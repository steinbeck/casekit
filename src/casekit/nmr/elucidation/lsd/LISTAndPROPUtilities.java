package casekit.nmr.elucidation.lsd;

import casekit.nmr.elucidation.model.MolecularConnectivity;
import casekit.nmr.model.Signal;
import casekit.nmr.utils.Statistics;

import java.util.*;
import java.util.stream.Collectors;

public class LISTAndPROPUtilities {

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

    public static void insertGeneralLISTs(final StringBuilder stringBuilder, final Map<String, Object[]> listMap,
                                          final Map<Integer, List<MolecularConnectivity>> molecularConnectivityMap,
                                          final Set<String> atomTypesByMf) {
        final Set<String> atomTypes = new HashSet<>(atomTypesByMf);
        atomTypes.remove("H");
        Set<Integer> elementIndices;
        for (final String atomType : atomTypes) {
            listMap.put(atomType, new Object[]{"L"
                                                       + (listMap.size()
                    + 1)});
            elementIndices = new HashSet<>();
            for (final int correlationIndex : molecularConnectivityMap.keySet()) {
                elementIndices.addAll(molecularConnectivityMap.get(correlationIndex)
                                                              .stream()
                                                              .filter(molecularConnectivity -> molecularConnectivity.getAtomType()
                                                                                                                    .equals(atomType))
                                                              .map(MolecularConnectivity::getIndex)
                                                              .collect(Collectors.toSet()));
            }

            stringBuilder.append("LIST ")
                         .append(listMap.get(atomType)[0]);
            for (final int elementIndex : elementIndices) {
                stringBuilder.append(" ")
                             .append(elementIndex);
            }
            stringBuilder.append("; list of all ")
                         .append(atomType)
                         .append(" atoms")
                         .append("\n");
        }
    }

    private static String buildListKey(final String atomType, final Set<Integer> hybridizations,
                                       final Set<Integer> protonsCounts) {
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
                                                       final Map<Integer, List<MolecularConnectivity>> molecularConnectivityMap) {
        final Map<String, Set<Integer>> atomIndicesMap = new LinkedHashMap<>();
        String listKey;
        for (final int correlationIndex : molecularConnectivityMap.keySet()) {
            for (final MolecularConnectivity molecularConnectivity : molecularConnectivityMap.get(correlationIndex)) {
                if (molecularConnectivity.getAtomType()
                                         .equals("H")
                        //                        || molecularConnectivity.getHybridizations()
                        //                                                .size()
                        //                        != 1
                        || molecularConnectivity.getProtonCounts()
                                                .size()
                        != 1) {
                    continue;
                }
                listKey = buildListKey(molecularConnectivity.getAtomType(), new HashSet<>(),
                                       // correlation.getHybridization(),
                                       molecularConnectivity.getProtonCounts());
                atomIndicesMap.putIfAbsent(listKey, new HashSet<>());
                atomIndicesMap.get(listKey)
                              .add(molecularConnectivity.getIndex());
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
                                                  final Map<String, Integer> usedPropsCount, final String listKey,
                                                  final String mode) {
        if (!listMap.containsKey(listKey)) {
            return true;
        }
        // LSD crashes if we try to use more atoms with specific hybridization or/and proton count (in mode "allow")
        // thus count and limit that usage as following:
        // list name -> number of usages
        usedPropsCount.putIfAbsent((String) listMap.get(listKey)[0], 1);
        return mode.equals("allow")
                && listMap.get(listKey).length
                > 1
                && usedPropsCount.get((String) listMap.get(listKey)[0])
                > (int) listMap.get(listKey)[1];
    }

    private static void insertPROP(final StringBuilder stringBuilder, final Map<String, Object[]> listMap,
                                   final String atomType, final Signal signal, final int indexInPyLSD,
                                   final String listKey, final String mode) {
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
                                                     final Map<Integer, List<MolecularConnectivity>> molecularConnectivityMap,
                                                     final String mode) {
        String listKey;
        Map<String, Map<Integer, Set<Integer>>> neighborsTemp;
        final Map<String, Integer> usedPropsCount = new HashMap<>();
        for (final int correlationIndex : molecularConnectivityMap.keySet()) {
            for (final MolecularConnectivity molecularConnectivity : molecularConnectivityMap.get(correlationIndex)) {
                if (mode.equals("forbid")) {
                    neighborsTemp = molecularConnectivity.getForbiddenNeighbors();
                } else if (mode.equals("allow")) {
                    neighborsTemp = molecularConnectivity.getSetNeighbors();
                } else {
                    neighborsTemp = null;
                }
                if (neighborsTemp
                        != null) {
                    for (final String neighborAtomType : neighborsTemp.keySet()) {
                        // forbid/set bonds to whole element groups if there is an empty map for an atom type
                        if (neighborsTemp.get(neighborAtomType)
                                         .isEmpty()) {
                            insertPROP(stringBuilder, listMap, molecularConnectivity.getAtomType(),
                                       molecularConnectivity.getSignal(), molecularConnectivity.getIndex(),
                                       neighborAtomType, mode);
                        } else {
                            for (final int neighborHybridization : neighborsTemp.get(neighborAtomType)
                                                                                .keySet()) {
                                for (final int protonsCount : neighborsTemp.get(neighborAtomType)
                                                                           .get(neighborHybridization)) {
                                    listKey = buildListKey(neighborAtomType, neighborHybridization
                                                                                     == -1
                                                                             ? new HashSet<>()
                                                                             : Set.of(neighborHybridization),
                                                           Set.of(protonsCount));
                                    if (checkSkipPROPInsertion(listMap, usedPropsCount, listKey, mode)) {
                                        continue;
                                    }
                                    if (listMap.containsKey(listKey)) {
                                        insertPROP(stringBuilder, listMap, molecularConnectivity.getAtomType(),
                                                   molecularConnectivity.getSignal(), molecularConnectivity.getIndex(),
                                                   listKey, mode);
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
