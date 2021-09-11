package casekit.nmr.lsd;

import casekit.nmr.model.nmrium.Correlation;

import java.util.*;
import java.util.stream.Collectors;

public class Utilities {

    public static void reduceDefaultHybridizationsAndProtonCountsOfHeteroAtoms(final List<Correlation> correlationList,
                                                                               final Map<Integer, Map<String, Map<String, Set<Integer>>>> detectedConnectivities) {
        final Map<String, Set<String>> allowedHeteroAtomHybridizations = buildAllowedHeteroAtomHybridizations(
                correlationList, detectedConnectivities);
        final Map<String, Set<Integer>> allowedHeteroAtomProtonCounts = buildAllowedHeteroAtomProtonCounts(
                correlationList, detectedConnectivities);
        // hetero atoms can bond to carbons only, due to that we can use more connectivity information
        // do not allow bond between carbon and hetero atoms in certain hybridization states and proton counts
        for (final Correlation correlation : correlationList) {
            // ignore C and H atoms
            if (correlation.getAtomType()
                           .equals("C")
                    || correlation.getAtomType()
                                  .equals("H")) {
                continue;
            }
            // but only if we have seen the hetero atom type in connectivity statistics
            // and hybridization states or protons count was not set beforehand
            if (correlation.getEdited()
                    != null
                    && correlation.getEdited()
                                  .containsKey("hybridization")
                    && correlation.getEdited()
                                  .get("hybridization")
                    && allowedHeteroAtomHybridizations.containsKey(correlation.getAtomType())) {
                correlation.getHybridization()
                           .retainAll(allowedHeteroAtomHybridizations.get(correlation.getAtomType()));
            }
            if (correlation.getEdited()
                    != null
                    && correlation.getEdited()
                                  .containsKey("protonsCount")
                    && correlation.getEdited()
                                  .get("protonsCount")
                    && allowedHeteroAtomProtonCounts.containsKey(correlation.getAtomType())) {
                correlation.getProtonsCount()
                           .retainAll(allowedHeteroAtomProtonCounts.get(correlation.getAtomType()));
            }
        }
    }

    public static Map<String, Set<String>> buildAllowedHeteroAtomHybridizations(final List<Correlation> correlationList,
                                                                                final Map<Integer, Map<String, Map<String, Set<Integer>>>> detectedConnectivities) {
        final Map<String, Set<String>> allowedHeteroAtomHybridizations = new HashMap<>();
        for (final Map.Entry<Integer, Map<String, Map<String, Set<Integer>>>> correlationEntry : detectedConnectivities.entrySet()) {
            if (!correlationList.get(correlationEntry.getKey())
                                .getAtomType()
                                .equals("C")
                    && !correlationList.get(correlationEntry.getKey())
                                       .getAtomType()
                                       .equals("H")) {
                continue;
            }
            for (final Map.Entry<String, Map<String, Set<Integer>>> neighborAtomTypeEntry : correlationEntry.getValue()
                                                                                                            .entrySet()) {
                allowedHeteroAtomHybridizations.putIfAbsent(neighborAtomTypeEntry.getKey(), new HashSet<>());
                allowedHeteroAtomHybridizations.get(neighborAtomTypeEntry.getKey())
                                               .addAll(neighborAtomTypeEntry.getValue()
                                                                            .keySet());
            }
        }

        return allowedHeteroAtomHybridizations;
    }

    public static Map<String, Set<Integer>> buildAllowedHeteroAtomProtonCounts(final List<Correlation> correlationList,
                                                                               final Map<Integer, Map<String, Map<String, Set<Integer>>>> detectedConnectivities) {
        final Map<String, Set<Integer>> allowedHeteroAtomProtonCounts = new HashMap<>();
        for (final Map.Entry<Integer, Map<String, Map<String, Set<Integer>>>> correlationEntry : detectedConnectivities.entrySet()) {
            if (!correlationList.get(correlationEntry.getKey())
                                .getAtomType()
                                .equals("C")
                    && !correlationList.get(correlationEntry.getKey())
                                       .getAtomType()
                                       .equals("H")) {
                continue;
            }
            for (final Map.Entry<String, Map<String, Set<Integer>>> neighborAtomTypeEntry : correlationEntry.getValue()
                                                                                                            .entrySet()) {
                allowedHeteroAtomProtonCounts.putIfAbsent(neighborAtomTypeEntry.getKey(), new HashSet<>());
                for (final Map.Entry<String, Set<Integer>> neighborHybridizationEntry : neighborAtomTypeEntry.getValue()
                                                                                                             .entrySet()) {
                    allowedHeteroAtomProtonCounts.get(neighborAtomTypeEntry.getKey())
                                                 .addAll(neighborHybridizationEntry.getValue());
                }
            }
        }

        return allowedHeteroAtomProtonCounts;
    }

    public static Map<String, Map<Integer, Set<Integer>>> buildForbiddenNeighborHybridizationsAndProtonCounts(
            final Map<String, Map<String, Set<Integer>>> connectivities, final Set<String> neighborAtomTypes) {

        // define forbidden hybridizations and proton counts (carbons only) of possible neighbors
        final Map<String, Map<Integer, Set<Integer>>> forbiddenNeighborHybridizationsAndProtonCounts = new HashMap<>();
        for (final String neighborAtomType : neighborAtomTypes) {
            forbiddenNeighborHybridizationsAndProtonCounts.put(neighborAtomType, new HashMap<>());
            for (final int defaultHybridization : Arrays.stream(Constants.defaultHybridizationMap.get(neighborAtomType))
                                                        .boxed()
                                                        .collect(Collectors.toList())) {
                forbiddenNeighborHybridizationsAndProtonCounts.get(neighborAtomType)
                                                              .put(defaultHybridization, Arrays.stream(
                                                                      Constants.defaultProtonsCountPerValencyMap.get(
                                                                              neighborAtomType))
                                                                                               .boxed()
                                                                                               .collect(
                                                                                                       Collectors.toSet()));
            }
            for (final String neighborHybridization : connectivities.get(neighborAtomType)
                                                                    .keySet()) {
                // remove found protons count per hybridzations from list of forbidden ones
                for (final int forbiddenNeighborHybridization : new HashSet<>(
                        forbiddenNeighborHybridizationsAndProtonCounts.get(neighborAtomType)
                                                                      .keySet())) {
                    forbiddenNeighborHybridizationsAndProtonCounts.get(neighborAtomType)
                                                                  .get(forbiddenNeighborHybridization)
                                                                  .removeAll(connectivities.get(neighborAtomType)
                                                                                           .get(neighborHybridization));
                    if (forbiddenNeighborHybridizationsAndProtonCounts.get(neighborAtomType)
                                                                      .get(forbiddenNeighborHybridization)
                                                                      .isEmpty()) {
                        forbiddenNeighborHybridizationsAndProtonCounts.get(neighborAtomType)
                                                                      .remove(forbiddenNeighborHybridization);
                    }
                }
                if (forbiddenNeighborHybridizationsAndProtonCounts.get(neighborAtomType)
                                                                  .isEmpty()) {
                    forbiddenNeighborHybridizationsAndProtonCounts.remove(neighborAtomType);
                    break;
                }
            }
        }

        return forbiddenNeighborHybridizationsAndProtonCounts;
    }
}
