package casekit.nmr.lsd;

import java.util.*;
import java.util.stream.Collectors;

public class Utilities {

    //    public static void reduceDefaultHybridizationsAndProtonCountsOfHeteroAtoms(final List<Correlation> correlationList,
    //                                                                               final Map<Integer, Map<String, Map<String, Set<Integer>>>> detectedConnectivities) {
    //        final Map<String, Set<String>> allowedNeighborAtomHybridizations = buildAllowedNeighborAtomHybridizations(
    //                correlationList, detectedConnectivities);
    //        final Map<String, Set<Integer>> allowedNeighborAtomProtonCounts = buildAllowedNeighborAtomProtonCounts(
    //                correlationList, detectedConnectivities);
    //        // hetero atoms can bond to carbons only, due to that we can use further connectivity information
    //        // do not allow bond between carbon and hetero atoms in certain hybridization states and proton counts
    //        for (final Correlation correlation : correlationList) {
    //            // ignore C and H atoms
    //            if (correlation.getAtomType()
    //                           .equals("C")
    //                    || correlation.getAtomType()
    //                                  .equals("H")) {
    //                continue;
    //            }
    //            // but only if we have seen the hetero atom type in connectivity statistics
    //            // and hybridization states or protons count was not set beforehand
    //            if (correlation.getHybridization()
    //                           .isEmpty()) {
    //                correlation.getHybridization()
    //                           .addAll(allowedNeighborAtomHybridizations.get(correlation.getAtomType()));
    //            } else if (correlation.getEdited()
    //                    != null
    //                    && correlation.getEdited()
    //                                  .containsKey("hybridization")
    //                    && !correlation.getEdited()
    //                                   .get("hybridization")
    //                    && allowedNeighborAtomHybridizations.containsKey(correlation.getAtomType())) {
    //                correlation.getHybridization()
    //                           .retainAll(allowedNeighborAtomHybridizations.get(correlation.getAtomType()));
    //            }
    //            if (correlation.getProtonsCount()
    //                           .isEmpty()) {
    //                correlation.getProtonsCount()
    //                           .addAll(allowedNeighborAtomProtonCounts.get(correlation.getAtomType()));
    //            } else if (correlation.getEdited()
    //                    != null
    //                    && correlation.getEdited()
    //                                  .containsKey("protonsCount")
    //                    && !correlation.getEdited()
    //                                   .get("protonsCount")
    //                    && allowedNeighborAtomProtonCounts.containsKey(correlation.getAtomType())) {
    //                correlation.getProtonsCount()
    //                           .retainAll(allowedNeighborAtomProtonCounts.get(correlation.getAtomType()));
    //            }
    //        }
    //    }
    //
    //    public static Map<String, Set<String>> buildAllowedNeighborAtomHybridizations(
    //            final List<Correlation> correlationList,
    //            final Map<Integer, Map<String, Map<String, Set<Integer>>>> detectedConnectivities) {
    //        final Map<String, Set<String>> allowedHeteroAtomHybridizations = new HashMap<>();
    //        for (final Map.Entry<Integer, Map<String, Map<String, Set<Integer>>>> correlationEntry : detectedConnectivities.entrySet()) {
    //            if (!correlationList.get(correlationEntry.getKey())
    //                                .getAtomType()
    //                                .equals("C")
    //                    && !correlationList.get(correlationEntry.getKey())
    //                                       .getAtomType()
    //                                       .equals("H")) {
    //                continue;
    //            }
    //            for (final Map.Entry<String, Map<String, Set<Integer>>> neighborAtomTypeEntry : correlationEntry.getValue()
    //                                                                                                            .entrySet()) {
    //                allowedHeteroAtomHybridizations.putIfAbsent(neighborAtomTypeEntry.getKey(), new HashSet<>());
    //                allowedHeteroAtomHybridizations.get(neighborAtomTypeEntry.getKey())
    //                                               .addAll(neighborAtomTypeEntry.getValue()
    //                                                                            .keySet());
    //            }
    //        }
    //
    //        return allowedHeteroAtomHybridizations;
    //    }
    //
    //    public static Map<String, Set<Integer>> buildAllowedNeighborAtomProtonCounts(
    //            final List<Correlation> correlationList,
    //            final Map<Integer, Map<String, Map<String, Set<Integer>>>> detectedConnectivities) {
    //        final Map<String, Set<Integer>> allowedHeteroAtomProtonCounts = new HashMap<>();
    //        for (final Map.Entry<Integer, Map<String, Map<String, Set<Integer>>>> correlationEntry : detectedConnectivities.entrySet()) {
    //            if (!correlationList.get(correlationEntry.getKey())
    //                                .getAtomType()
    //                                .equals("C")
    //                    && !correlationList.get(correlationEntry.getKey())
    //                                       .getAtomType()
    //                                       .equals("H")) {
    //                continue;
    //            }
    //            for (final Map.Entry<String, Map<String, Set<Integer>>> neighborAtomTypeEntry : correlationEntry.getValue()
    //                                                                                                            .entrySet()) {
    //                allowedHeteroAtomProtonCounts.putIfAbsent(neighborAtomTypeEntry.getKey(), new HashSet<>());
    //                for (final Map.Entry<String, Set<Integer>> neighborHybridizationEntry : neighborAtomTypeEntry.getValue()
    //                                                                                                             .entrySet()) {
    //                    allowedHeteroAtomProtonCounts.get(neighborAtomTypeEntry.getKey())
    //                                                 .addAll(neighborHybridizationEntry.getValue());
    //                }
    //            }
    //        }
    //
    //        return allowedHeteroAtomProtonCounts;
    //    }

    public static Map<String, Set<Integer>> buildForbiddenNeighbors(final Map<String, Set<Integer>> connectivities,
                                                                    final Set<String> possibleNeighborAtomTypes) {

        // define forbidden neighbors (carbons only)
        // or put just an empty map which stands for the whole element
        final Map<String, Set<Integer>> forbiddenNeighbors = new HashMap<>();
        for (final String possibleNeighborAtomType : possibleNeighborAtomTypes) {
            if (possibleNeighborAtomType.equals("H")) {
                continue;
            }
            forbiddenNeighbors.put(possibleNeighborAtomType, new HashSet<>());
            if (connectivities.containsKey(possibleNeighborAtomType)) {

                forbiddenNeighbors.get(possibleNeighborAtomType)
                                  .addAll(Arrays.stream(
                                                        Constants.defaultProtonsCountPerValencyMap.get(possibleNeighborAtomType))
                                                .boxed()
                                                .collect(Collectors.toSet()));

                for (final int neighborProtonCount : connectivities.get(possibleNeighborAtomType)) {
                    // remove found protons count from list of forbidden ones
                    forbiddenNeighbors.get(possibleNeighborAtomType)
                                      .remove(neighborProtonCount);
                    if (forbiddenNeighbors.get(possibleNeighborAtomType)
                                          .isEmpty()) {
                        forbiddenNeighbors.remove(possibleNeighborAtomType);
                        break;
                    }
                }
            }
        }

        return forbiddenNeighbors;
    }
}
