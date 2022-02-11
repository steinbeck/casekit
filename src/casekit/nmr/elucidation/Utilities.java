package casekit.nmr.elucidation;

import casekit.nmr.elucidation.model.Detections;
import casekit.nmr.elucidation.model.Grouping;
import casekit.nmr.elucidation.model.MolecularConnectivity;
import casekit.nmr.model.PathLength;
import casekit.nmr.model.Signal;
import casekit.nmr.model.nmrium.Correlation;
import casekit.nmr.model.nmrium.Link;
import casekit.nmr.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

public class Utilities {

    public static void reduceDefaultHybridizationsAndProtonCountsOfHeteroAtoms(final List<Correlation> correlationList,
                                                                               final Map<Integer, Map<String, Map<Integer, Set<Integer>>>> detectedConnectivities,
                                                                               final Map<Integer, List<Integer>> detectedHybridizations) {
        if (detectedConnectivities
                == null
                || detectedConnectivities.isEmpty()) {
            return;
        }
        final Map<String, Set<Integer>> allowedNeighborAtomHybridizations = buildAllowedNeighborAtomHybridizations(
                correlationList, detectedConnectivities);
        final Map<String, Set<Integer>> allowedNeighborAtomProtonCounts = buildAllowedNeighborAtomProtonCounts(
                correlationList, detectedConnectivities);
        // hetero atoms can bond to carbons only, due to that we can use further connectivity information
        // do not allow bond between carbon and hetero atoms in certain hybridization states and proton counts
        Correlation correlation;
        for (int i = 0; i
                < correlationList.size(); i++) {
            correlation = correlationList.get(i);
            // ignore C and H atoms
            if (correlation.getAtomType()
                           .equals("C")
                    || correlation.getAtomType()
                                  .equals("H")) {
                continue;
            }
            final Set<Integer> hybridizationsToAdd = allowedNeighborAtomHybridizations.containsKey(
                    correlation.getAtomType())
                                                     ? allowedNeighborAtomHybridizations.get(correlation.getAtomType())
                                                     : Arrays.stream(Constants.defaultHybridizationMap.get(
                                                                     correlation.getAtomType()))
                                                             .boxed()
                                                             .collect(Collectors.toSet());
            final Set<Integer> protonCountsToAdd = allowedNeighborAtomProtonCounts.containsKey(
                    correlation.getAtomType())
                                                   ? allowedNeighborAtomProtonCounts.get(correlation.getAtomType())
                                                   : Arrays.stream(Constants.defaultProtonsCountPerValencyMap.get(
                                                                   correlation.getAtomType()))
                                                           .boxed()
                                                           .collect(Collectors.toSet());
            // but only if we have seen the hetero atom type in connectivity statistics
            // and hybridization states or protons count was not set beforehand
            if (correlation.getHybridization()
                           .isEmpty()) {
                correlation.getHybridization()
                           .addAll(hybridizationsToAdd);
            } else if (correlation.getEdited()
                    != null
                    && correlation.getEdited()
                                  .containsKey("hybridization")
                    && !correlation.getEdited()
                                   .get("hybridization")
                    && allowedNeighborAtomHybridizations.containsKey(correlation.getAtomType())) {
                correlation.getHybridization()
                           .retainAll(hybridizationsToAdd);
            }
            if (correlation.getProtonsCount()
                           .isEmpty()) {
                correlation.getProtonsCount()
                           .addAll(protonCountsToAdd);
            } else if (correlation.getEdited()
                    != null
                    && correlation.getEdited()
                                  .containsKey("protonsCount")
                    && !correlation.getEdited()
                                   .get("protonsCount")) {
                correlation.getProtonsCount()
                           .retainAll(protonCountsToAdd);
            }
            detectedHybridizations.putIfAbsent(i, new ArrayList<>());
            detectedHybridizations.get(i)
                                  .addAll(correlation.getHybridization());
        }
    }

    public static Map<String, Set<Integer>> buildAllowedNeighborAtomHybridizations(
            final List<Correlation> correlationList,
            final Map<Integer, Map<String, Map<Integer, Set<Integer>>>> detectedConnectivities) {
        final Map<String, Set<Integer>> allowedHeteroAtomHybridizations = new HashMap<>();
        for (final Map.Entry<Integer, Map<String, Map<Integer, Set<Integer>>>> correlationEntry : detectedConnectivities.entrySet()) {
            if (!correlationList.get(correlationEntry.getKey())
                                .getAtomType()
                                .equals("C")
                    && !correlationList.get(correlationEntry.getKey())
                                       .getAtomType()
                                       .equals("H")) {
                continue;
            }
            for (final Map.Entry<String, Map<Integer, Set<Integer>>> neighborAtomTypeEntry : correlationEntry.getValue()
                                                                                                             .entrySet()) {
                allowedHeteroAtomHybridizations.putIfAbsent(neighborAtomTypeEntry.getKey(), new HashSet<>());
                allowedHeteroAtomHybridizations.get(neighborAtomTypeEntry.getKey())
                                               .addAll(neighborAtomTypeEntry.getValue()
                                                                            .keySet());
            }
        }

        return allowedHeteroAtomHybridizations;
    }

    public static Map<String, Set<Integer>> buildAllowedNeighborAtomProtonCounts(
            final List<Correlation> correlationList,
            final Map<Integer, Map<String, Map<Integer, Set<Integer>>>> detectedConnectivities) {
        final Map<String, Set<Integer>> allowedHeteroAtomProtonCounts = new HashMap<>();
        for (final Map.Entry<Integer, Map<String, Map<Integer, Set<Integer>>>> correlationEntry : detectedConnectivities.entrySet()) {
            if (!correlationList.get(correlationEntry.getKey())
                                .getAtomType()
                                .equals("C")
                    && !correlationList.get(correlationEntry.getKey())
                                       .getAtomType()
                                       .equals("H")) {
                continue;
            }
            for (final Map.Entry<String, Map<Integer, Set<Integer>>> neighborAtomTypeEntry : correlationEntry.getValue()
                                                                                                             .entrySet()) {
                allowedHeteroAtomProtonCounts.putIfAbsent(neighborAtomTypeEntry.getKey(), new HashSet<>());
                for (final Map.Entry<Integer, Set<Integer>> neighborHybridizationEntry : neighborAtomTypeEntry.getValue()
                                                                                                              .entrySet()) {
                    allowedHeteroAtomProtonCounts.get(neighborAtomTypeEntry.getKey())
                                                 .addAll(neighborHybridizationEntry.getValue());
                }
            }
        }

        return allowedHeteroAtomProtonCounts;
    }

    public static Map<String, Map<Integer, Set<Integer>>> buildForbiddenNeighbors(
            final Map<String, Map<Integer, Set<Integer>>> connectivities, final Set<String> possibleNeighborAtomTypes) {

        // define forbidden neighbors (for carbons only)
        // or put just an empty map which means the whole element is forbidden
        final Map<String, Map<Integer, Set<Integer>>> forbiddenNeighbors = new HashMap<>();
        for (final String possibleNeighborAtomType : possibleNeighborAtomTypes) {
            if (possibleNeighborAtomType.equals("H")) {
                continue;
            }
            forbiddenNeighbors.put(possibleNeighborAtomType, new HashMap<>());
            if (connectivities.containsKey(possibleNeighborAtomType)) {
                for (final int defaultHybridization : Arrays.stream(
                                                                    Constants.defaultHybridizationMap.get(possibleNeighborAtomType))
                                                            .boxed()
                                                            .collect(Collectors.toList())) {
                    forbiddenNeighbors.get(possibleNeighborAtomType)
                                      .put(defaultHybridization, Arrays.stream(
                                                                               Constants.defaultProtonsCountPerValencyMap.get(possibleNeighborAtomType))
                                                                       .boxed()
                                                                       .collect(Collectors.toSet()));
                }
                for (final int possibleNeighborHybridization : connectivities.get(possibleNeighborAtomType)
                                                                             .keySet()) {
                    // remove found protons count per hybridzations from list of forbidden ones
                    for (final int forbiddenNeighborHybridization : new HashSet<>(
                            forbiddenNeighbors.get(possibleNeighborAtomType)
                                              .keySet())) {
                        forbiddenNeighbors.get(possibleNeighborAtomType)
                                          .get(forbiddenNeighborHybridization)
                                          .removeAll(connectivities.get(possibleNeighborAtomType)
                                                                   .get(possibleNeighborHybridization));
                        if (forbiddenNeighbors.get(possibleNeighborAtomType)
                                              .get(forbiddenNeighborHybridization)
                                              .isEmpty()) {
                            forbiddenNeighbors.get(possibleNeighborAtomType)
                                              .remove(forbiddenNeighborHybridization);
                        }
                    }
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


    public static Map<Integer, Set<Integer>> buildFixedNeighborsByINADEQUATE(final List<Correlation> correlationList) {
        final Map<Integer, Set<Integer>> fixedNeighbors = new HashMap<>();
        final Set<String> uniqueSet = new HashSet<>();
        Correlation correlation;
        for (int i = 0; i
                < correlationList.size(); i++) {
            correlation = correlationList.get(i);
            // @TODO for now use INADEQUATE information of atoms without equivalences only
            if (correlation.getEquivalence()
                    != 1) {
                continue;
            }
            for (final Link link : correlation.getLink()) {
                if (link.getExperimentType()
                        .equals("inadequate")) {
                    for (final int matchIndex : link.getMatch()) {
                        // insert BOND pair once only and not if equivalences exist
                        if (!uniqueSet.contains(i
                                                        + " "
                                                        + matchIndex)
                                && correlationList.get(matchIndex)
                                                  .getEquivalence()
                                == 1) {
                            fixedNeighbors.putIfAbsent(i, new HashSet<>());
                            fixedNeighbors.get(i)
                                          .add(matchIndex);
                            uniqueSet.add(i
                                                  + " "
                                                  + matchIndex);
                            uniqueSet.add(matchIndex
                                                  + " "
                                                  + i);
                        }
                    }
                }
            }
        }

        return fixedNeighbors;
    }

    private static boolean hasMatch(final Correlation correlation1, final Correlation correlation2,
                                    final double tolerance) {
        final Signal signal1 = Utils.extractSignalFromCorrelation(correlation1);
        final Signal signal2 = Utils.extractSignalFromCorrelation(correlation2);
        if (signal1
                == null
                || signal2
                == null) {
            return false;
        }
        int dim1 = -1;
        int dim2 = -1;
        String atomType;
        for (int i = 0; i
                < signal1.getNuclei().length; i++) {
            atomType = Utils.getAtomTypeFromNucleus(signal1.getNuclei()[i]);
            if (atomType.equals(correlation1.getAtomType())) {
                dim1 = i;
                break;
            }
        }
        for (int i = 0; i
                < signal2.getNuclei().length; i++) {
            atomType = Utils.getAtomTypeFromNucleus(signal2.getNuclei()[i]);
            if (atomType.equals(correlation2.getAtomType())) {
                dim2 = i;
                break;
            }
        }
        if (dim1
                == -1
                || dim2
                == -1) {
            return false;
        }

        final double shift1 = signal1.getShift(dim1);
        final double shift2 = signal2.getShift(dim2);

        return Math.abs(shift1
                                - shift2)
                <= tolerance;

    }

    private static Map<String, Map<Integer, List<Integer>>> findGroups(final List<Correlation> correlationList,
                                                                       final Map<String, Double> tolerances) {
        // cluster group index -> list of correlation index pair
        final Map<String, Map<Integer, List<Integer>>> groups = new HashMap<>();
        int groupIndex = 0;
        final Set<Integer> inserted = new HashSet<>();
        int foundGroupIndex;
        for (int i = 0; i
                < correlationList.size(); i++) {
            final Correlation correlation = correlationList.get(i);
            if (inserted.contains(i)
                    || correlation.isPseudo()) {
                continue;
            }
            groups.putIfAbsent(correlation.getAtomType(), new HashMap<>());
            // if we have a match somewhere then add the correlation index into to group
            // if not then create a new group
            foundGroupIndex = -1;
            for (final Map.Entry<Integer, List<Integer>> groupEntry : groups.get(correlation.getAtomType())
                                                                            .entrySet()) {
                if (groupEntry.getValue()
                              .stream()
                              .anyMatch(correlationIndex -> hasMatch(correlation, correlationList.get(correlationIndex),
                                                                     tolerances.get(correlation.getAtomType())))) {
                    foundGroupIndex = groupEntry.getKey();
                    break;
                }
            }
            if (foundGroupIndex
                    != -1) {
                groups.get(correlation.getAtomType())
                      .get(foundGroupIndex)
                      .add(i);
                inserted.add(i);
            } else {
                groups.get(correlation.getAtomType())
                      .put(groupIndex, new ArrayList<>());
                groups.get(correlation.getAtomType())
                      .get(groupIndex)
                      .add(i);
                inserted.add(i);
                groupIndex++;
            }
        }

        return groups;
    }

    private static Map<String, Map<Integer, Integer>> transformGroups(
            final Map<String, Map<Integer, List<Integer>>> groups) {
        final Map<String, Map<Integer, Integer>> transformedGroups = new HashMap<>();
        for (final Map.Entry<String, Map<Integer, List<Integer>>> atomTypeEntry : groups.entrySet()) {
            transformedGroups.put(atomTypeEntry.getKey(), new HashMap<>());
            for (final Map.Entry<Integer, List<Integer>> groupEntry : atomTypeEntry.getValue()
                                                                                   .entrySet()) {
                for (final int correlationIndex : groupEntry.getValue()) {
                    transformedGroups.get(atomTypeEntry.getKey())
                                     .put(correlationIndex, groupEntry.getKey());
                }
            }
        }

        return transformedGroups;
    }

    public static Grouping buildGroups(final List<Correlation> correlationList, final Map<String, Double> tolerances) {
        final Map<String, Map<Integer, List<Integer>>> groups = findGroups(correlationList, tolerances);

        return new Grouping(tolerances, groups, transformGroups(groups));
    }

    private static List<Integer> getProtonCounts(final List<Correlation> correlationList, final int index) {
        final Correlation correlation = correlationList.get(index);
        if (correlation.getProtonsCount()
                != null
                && !correlation.getProtonsCount()
                               .isEmpty()) {
            // if protonCounts is already given
            return correlation.getProtonsCount();
        }
        final List<Integer> protonCounts = new ArrayList<>();
        final int[] defaultProtonCounts = Constants.defaultProtonsCountPerValencyMap.get(
                Constants.defaultAtomLabelMap.get(correlation.getAtomType()));
        for (final int defaultProtonCount : defaultProtonCounts) {
            protonCounts.add(defaultProtonCount);
        }

        return protonCounts;
    }

    private static List<Integer> getHybridizations(final List<Correlation> correlationList, final int index,
                                                   final Map<Integer, List<Integer>> detectedHybridizations) {
        final Correlation correlation = correlationList.get(index);
        List<Integer> hybridizations = new ArrayList<>();
        if (correlation.getHybridization()
                != null
                && !correlation.getHybridization()
                               .isEmpty()) {
            // if hybridization is already given
            return correlation.getHybridization();
        } else {
            // if hybridization is not given then use the detected ones
            if (detectedHybridizations.containsKey(index)) {
                hybridizations = new ArrayList<>(detectedHybridizations.get(index));
            }
            if (hybridizations.isEmpty()
                    && correlation.getAtomType()
                                  .equals("C")
                    && correlation.getProtonsCount()
                                  .size()
                    == 1
                    && correlation.getProtonsCount()
                                  .get(0)
                    >= 2) {
                // a carbon with at least two protons can only be SP2 or SP3
                hybridizations.add(2);
                hybridizations.add(3);
            }
        }
        if (hybridizations.isEmpty()) {
            for (int i = 0; i
                    < Constants.defaultHybridizationMap.get(correlation.getAtomType()).length; i++) {
                hybridizations.add(Constants.defaultHybridizationMap.get(correlation.getAtomType())[i]);
            }
        }

        return hybridizations;
    }

    //    private static List<Integer> getAttachedProtonIndices(final Correlation correlation) {
    //        return correlation.getLink()
    //                          .stream()
    //                          .filter(link -> link.getExperimentType()
    //                                              .equals("hsqc")
    //                                  || link.getExperimentType()
    //                                         .equals("hmqc"))
    //                          .map(Link::getMatch)
    //                          .reduce(new ArrayList<>(), ((l, curr) -> {
    //                              l.addAll(curr);
    //                              return l;
    //                          }));
    //    }

    //    public static Map<Integer, Integer[]> buildIndicesMap(final List<Correlation> correlationList) {
    //        // index in correlation data -> [indices in PyLSD file...]
    //        final Map<Integer, Integer[]> indicesMap = new HashMap<>();
    //        // init element indices within correlations with same order as in correlation data input
    //        int heavyAtomIndexInPyLSDFile = 1;
    //        int protonIndexInPyLSDFile = 1;
    //        int protonsToInsert;
    //        Correlation correlation;
    //        List<Integer> attachedProtonIndices;
    //        Integer[] arrayToInsert;
    //        for (int i = 0; i
    //                < correlationList.size(); i++) {
    //            correlation = correlationList.get(i);
    //            // set entry for each correlation with consideration of equivalences
    //            if (!correlation.getAtomType()
    //                            .equals("H")) {
    //                // insert for protons
    //                attachedProtonIndices = getAttachedProtonIndices(correlation);
    //                protonsToInsert = correlation.getEquivalence();
    //                arrayToInsert = new Integer[protonsToInsert];
    //                for (int j = 0; j
    //                        < arrayToInsert.length; j++) {
    //                    arrayToInsert[j] = protonIndexInPyLSDFile;
    //                    protonIndexInPyLSDFile++;
    //                }
    //                for (final int attachedProtonIndex : attachedProtonIndices) {
    //                    indicesMap.put(attachedProtonIndex, arrayToInsert);
    //                }
    //                // insert for heavy atom itself
    //                indicesMap.put(i, new Integer[correlation.getEquivalence()]);
    //                for (int j = 0; j
    //                        < correlation.getEquivalence(); j++) {
    //                    indicesMap.get(i)[j] = heavyAtomIndexInPyLSDFile;
    //                    heavyAtomIndexInPyLSDFile++;
    //                }
    //            }
    //        }
    //
    //        System.out.println("\n -> indicesMap:\n");
    //        for (final int index : indicesMap.keySet()) {
    //            System.out.println("i: "
    //                                       + index
    //                                       + " -> "
    //                                       + Arrays.toString(indicesMap.get(index)));
    //        }
    //        System.out.println("\n");
    //
    //        return indicesMap;
    //    }

    public static Map<Integer, Integer[]> buildIndicesMap(final List<MolecularConnectivity> molecularConnectivityList) {
        // index in correlation data -> [indices in PyLSD file...]
        final Map<Integer, Integer[]> indicesMap = new HashMap<>();
        // init element indices within correlations with same order as in correlation data input
        int heavyAtomIndexInPyLSDFile = 1;
        int protonIndexInPyLSDFile = 1;
        int protonsToInsert;
        Integer[] arrayToInsert;
        for (final MolecularConnectivity molecularConnectivity : molecularConnectivityList) {
            // set entry for each correlation with consideration of equivalences
            if (!molecularConnectivity.getAtomType()
                                      .equals("H")) {
                if (molecularConnectivity.getHsqc()
                        != null) {
                    // insert for protons
                    protonsToInsert = molecularConnectivity.getEquivalence();
                    arrayToInsert = new Integer[protonsToInsert];
                    for (int j = 0; j
                            < arrayToInsert.length; j++) {
                        arrayToInsert[j] = protonIndexInPyLSDFile;
                        protonIndexInPyLSDFile++;
                    }
                    for (final int attachedProtonIndex : molecularConnectivity.getHsqc()) {
                        indicesMap.put(attachedProtonIndex, arrayToInsert);
                    }
                }
                // insert for heavy atom itself
                indicesMap.put(molecularConnectivity.getIndex(), new Integer[molecularConnectivity.getEquivalence()]);
                for (int j = 0; j
                        < molecularConnectivity.getEquivalence(); j++) {
                    indicesMap.get(molecularConnectivity.getIndex())[j] = heavyAtomIndexInPyLSDFile;
                    heavyAtomIndexInPyLSDFile++;
                }
            }
        }

        System.out.println("\n -> indicesMap:\n");
        for (final int index : indicesMap.keySet()) {
            System.out.println("i: "
                                       + index
                                       + " -> "
                                       + Arrays.toString(indicesMap.get(index)));
        }
        System.out.println("\n");

        return indicesMap;
    }

    public static MolecularConnectivity findMolecularConnectivityByIndex(
            final Map<Integer, List<MolecularConnectivity>> molecularConnectivityMap, final String atomType,
            final boolean exclude, final int index) {
        for (final int correlationIndex : molecularConnectivityMap.keySet()) {
            for (final MolecularConnectivity molecularConnectivity : molecularConnectivityMap.get(correlationIndex)) {
                if (((exclude
                        && !molecularConnectivity.getAtomType()
                                                 .equals(atomType))
                        || (!exclude
                        && molecularConnectivity.getAtomType()
                                                .equals(atomType)))
                        && molecularConnectivity.getIndex()
                        == index) {
                    return molecularConnectivity;
                }
            }
        }

        return null;
    }

    private static List<Integer> buildPossibilities(final Map<Integer, Integer[]> indicesMap,
                                                    final List<MolecularConnectivity> molecularConnectivityList,
                                                    final int correlationIndex, final Grouping grouping) {
        final MolecularConnectivity molecularConnectivity = molecularConnectivityList.get(correlationIndex);
        // add possible indices from grouping
        final int groupIndex;
        final Set<Integer> possibilities = new HashSet<>();
        if (grouping.getTransformedGroups()
                    .containsKey(molecularConnectivity.getAtomType())
                && grouping.getTransformedGroups()
                           .get(molecularConnectivity.getAtomType())
                           .containsKey(correlationIndex)) {
            groupIndex = grouping.getTransformedGroups()
                                 .get(molecularConnectivity.getAtomType())
                                 .get(correlationIndex);
            for (final int groupMemberIndex : grouping.getGroups()
                                                      .get(molecularConnectivity.getAtomType())
                                                      .get(groupIndex)) {

                //                if (indicesMap.containsKey(groupMemberIndex)) {
                // add equivalence indices of group members
                possibilities.addAll(Arrays.asList(indicesMap.get(groupMemberIndex)));
                //                }
            }
        } else {
            // add for equivalences only
            possibilities.addAll(Arrays.asList(indicesMap.get(correlationIndex)));
        }

        return new ArrayList<>(possibilities);
    }

    //    private static void addMolecularConnectivity(
    //            final Map<Integer, List<MolecularConnectivity>> molecularConnectivityMap,
    //            final MolecularConnectivity molecularConnectivity, final int index) {
    //        molecularConnectivityMap.putIfAbsent(molecularConnectivity.getIndex(), new ArrayList<>());
    //        // do not insert duplicates
    //        if (molecularConnectivityMap.get(molecularConnectivity.getIndex())
    //                                    .stream()
    //                                    .noneMatch(molecularConnectivityTemp -> molecularConnectivityTemp.getIndex()
    //                                            == index)) {
    //            final MolecularConnectivity newMolecularConnectivity = Utils.cloneObject(molecularConnectivity,
    //                                                                                     MolecularConnectivity.class)
    //            newMolecularConnectivity.setIndex(index);
    //            molecularConnectivity.setGroupMembers(
    //                    buildGroupMembers(indicesMap, correlationList, correlationIndex, grouping));
    //
    //            molecularConnectivityMap.get(molecularConnectivity.getIndex())
    //                                    .add(newMolecularConnectivity);
    //        }
    //    }

    public static Map<Integer, List<MolecularConnectivity>> buildMolecularConnectivityMap(
            final List<MolecularConnectivity> molecularConnectivityList, final Detections detections,
            final Grouping grouping, final Map<String, Integer[]> defaultBondDistances) {

        final Map<Integer, Integer[]> indicesMap = buildIndicesMap(molecularConnectivityList);
        // correlation index -> [MolecularConnectivity]
        final Map<Integer, List<MolecularConnectivity>> molecularConnectivityMap = new HashMap<>();
        MolecularConnectivity newMolecularConnectivity;

        int index, correlationIndex, protonCorrelationIndex, protonIndex;
        PathLength pathLength;
        for (final MolecularConnectivity molecularConnectivity : molecularConnectivityList) {
            correlationIndex = molecularConnectivity.getIndex();
            // skip in case of non-linked proton which has no index in indices map
            if (!indicesMap.containsKey(correlationIndex)) {
                continue;
            }
            for (int k = 0; k
                    < indicesMap.get(correlationIndex).length; k++) {
                index = indicesMap.get(correlationIndex)[k];
                newMolecularConnectivity = new MolecularConnectivity();
                newMolecularConnectivity.setIndex(index);
                newMolecularConnectivity.setAtomType(molecularConnectivity.getAtomType());
                newMolecularConnectivity.setEquivalence(1);
                newMolecularConnectivity.setPseudo(molecularConnectivity.isPseudo());
                newMolecularConnectivity.setProtonCounts(molecularConnectivity.getProtonCounts());
                newMolecularConnectivity.setHybridizations(molecularConnectivity.getHybridizations());
                newMolecularConnectivity.setSignal(molecularConnectivity.getSignal());

                if (!molecularConnectivity.getAtomType()
                                          .equals("H")
                        && molecularConnectivity.getHsqc()
                        != null) {
                    // using the first proton correlation index from HSQC list is enough because show will direct to same heavy atom index
                    protonIndex = indicesMap.get(molecularConnectivity.getHsqc()
                                                                      .get(0))[k];
                    newMolecularConnectivity.setHsqc(new ArrayList<>());
                    newMolecularConnectivity.getHsqc()
                                            .add(protonIndex);
                }
                if (!molecularConnectivity.getAtomType()
                                          .equals("H")
                        && molecularConnectivity.getHmbc()
                        != null) {
                    pathLength = molecularConnectivity.getSignal()
                                                      .getPathLength();
                    for (final Map.Entry<Integer, Integer[]> entry : molecularConnectivity.getHmbc()
                                                                                          .entrySet()) {
                        protonCorrelationIndex = entry.getKey();
                        //                        // ignore linked H atoms without any attachment to a heavy atom
                        //                        if (molecularConnectivityList.get(entry.getKey())
                        //                                           .getAtomType()
                        //                                           .equals("H")
                        //                                && correlationList.get(matchIndex)
                        //                                                  .getAttachment()
                        //                                                  .keySet()
                        //                                                  .isEmpty()) {
                        //                            continue;
                        //                        }
                        if (indicesMap.containsKey(protonCorrelationIndex)) {
                            for (int l = 0; l
                                    < indicesMap.get(protonCorrelationIndex).length; l++) {
                                protonIndex = indicesMap.get(protonCorrelationIndex)[l];
                                if (newMolecularConnectivity.getHmbc()
                                        == null) {
                                    newMolecularConnectivity.setHmbc(new HashMap<>());
                                }
                                newMolecularConnectivity.getHmbc()
                                                        .put(protonIndex, pathLength
                                                                                  == null
                                                                          ? defaultBondDistances.get("hmbc")
                                                                          : new Integer[]{pathLength.getMin(),
                                                                                          pathLength.getMax()});
                            }
                        }
                    }
                }
                if (molecularConnectivity.getAtomType()
                                         .equals("H")
                        && molecularConnectivity.getCosy()
                        != null) {
                    pathLength = molecularConnectivity.getSignal()
                                                      .getPathLength();
                    for (final Map.Entry<Integer, Integer[]> entry : molecularConnectivity.getCosy()
                                                                                          .entrySet()) {
                        protonCorrelationIndex = entry.getKey();
                        if (indicesMap.containsKey(protonCorrelationIndex)
                                && k
                                < indicesMap.get(protonCorrelationIndex).length) {
                            protonIndex = indicesMap.get(protonCorrelationIndex)[k];
                            if (newMolecularConnectivity.getCosy()
                                    == null) {
                                newMolecularConnectivity.setCosy(new HashMap<>());
                            }
                            newMolecularConnectivity.getCosy()
                                                    .put(protonIndex, pathLength
                                                                              == null
                                                                      ? defaultBondDistances.get("cosy")
                                                                      : new Integer[]{pathLength.getMin(),
                                                                                      pathLength.getMax()});
                        }
                    }
                }

                molecularConnectivityMap.putIfAbsent(correlationIndex, new ArrayList<>());
                molecularConnectivityMap.get(correlationIndex)
                                        .add(newMolecularConnectivity);
            }
            // set detections
            for (final MolecularConnectivity molecularConnectivityTemp : molecularConnectivityMap.get(
                    correlationIndex)) {
                if (detections.getForbiddenNeighbors()
                              .containsKey(correlationIndex)) {
                    molecularConnectivityTemp.setForbiddenNeighbors(detections.getForbiddenNeighbors()
                                                                              .get(correlationIndex));
                }
                if (detections.getSetNeighbors()
                              .containsKey(correlationIndex)) {
                    molecularConnectivityTemp.setSetNeighbors(detections.getSetNeighbors()
                                                                        .get(correlationIndex));
                }
                molecularConnectivityTemp.setGroupMembers(
                        buildPossibilities(indicesMap, molecularConnectivityList, correlationIndex, grouping));
            }
            // fill in fixed neighbors
            if (molecularConnectivity.getEquivalence()
                    == 1
                    && detections.getFixedNeighbors()
                                 .containsKey(correlationIndex)) {
                MolecularConnectivity molecularConnectivityTemp;
                for (final int correlationIndex2 : detections.getFixedNeighbors()
                                                             .get(correlationIndex)) {
                    molecularConnectivityTemp = molecularConnectivityList.get(correlationIndex2);
                    // use fixed neighbor information of atoms without equivalence equals 1 only
                    if (molecularConnectivityTemp.getEquivalence()
                            > 1) {
                        continue;
                    }
                    index = indicesMap.get(correlationIndex)[0];
                    newMolecularConnectivity = findMolecularConnectivityByIndex(molecularConnectivityMap,
                                                                                molecularConnectivity.getAtomType(),
                                                                                false, index);
                    if (newMolecularConnectivity.getFixedNeighbors()
                            == null) {
                        newMolecularConnectivity.setFixedNeighbors(new ArrayList<>());
                    }
                    newMolecularConnectivity.getFixedNeighbors()
                                            .add(indicesMap.get(correlationIndex2)[0]);
                }
            }
        }

        // filter out HMBC or COSY correlation to itself
        for (final int correlationIndexTemp : molecularConnectivityMap.keySet()) {
            for (final MolecularConnectivity molecularConnectivityTemp : molecularConnectivityMap.get(
                    correlationIndexTemp)) {
                if (molecularConnectivityTemp.getHsqc()
                        != null) {
                    if (molecularConnectivityTemp.getHmbc()
                            != null) {
                        molecularConnectivityTemp.getHmbc()
                                                 .remove(molecularConnectivityTemp.getHsqc());
                    }
                    if (molecularConnectivityTemp.getCosy()
                            != null) {
                        molecularConnectivityTemp.getCosy()
                                                 .remove(molecularConnectivityTemp.getHsqc());
                    }
                }
            }
        }

        return molecularConnectivityMap;
    }

    public static List<MolecularConnectivity> buildMolecularConnectivityList(final List<Correlation> correlationList,
                                                                             final Detections detections,
                                                                             final Grouping grouping,
                                                                             final Map<String, Integer[]> defaultBondDistances) {
        final List<MolecularConnectivity> molecularConnectivityList = new ArrayList<>();
        Correlation correlation;
        int groupIndex;
        Map<String, Object> signal2DMap;
        Map<String, Object> pathLengthMap;
        MolecularConnectivity molecularConnectivity;
        for (int correlationIndex = 0; correlationIndex
                < correlationList.size(); correlationIndex++) {
            correlation = correlationList.get(correlationIndex);
            molecularConnectivity = new MolecularConnectivity();
            molecularConnectivity.setIndex(correlationIndex);
            molecularConnectivity.setAtomType(correlation.getAtomType());
            molecularConnectivity.setSignal(Utils.extractSignalFromCorrelation(correlation));
            molecularConnectivity.setEquivalence(correlation.getEquivalence());
            molecularConnectivity.setPseudo(correlation.isPseudo());

            if (!correlation.getAtomType()
                            .equals("H")) {
                molecularConnectivity.setProtonCounts(getProtonCounts(correlationList, correlationIndex));
                molecularConnectivity.setHybridizations(
                        getHybridizations(correlationList, correlationIndex, detections.getDetectedHybridizations()));
            }
            if (grouping.getGroups()
                        .containsKey(correlation.getAtomType())
                    && grouping.getTransformedGroups()
                               .get(correlation.getAtomType())
                               .containsKey(correlationIndex)) {
                groupIndex = grouping.getTransformedGroups()
                                     .get(correlation.getAtomType())
                                     .get(correlationIndex);
                molecularConnectivity.setGroupMembers(grouping.getGroups()
                                                              .get(correlation.getAtomType())
                                                              .get(groupIndex));
            }
            for (final Link link : correlation.getLink()) {
                if (link.getExperimentType()
                        .equals("hsqc")
                        || link.getExperimentType()
                               .equals("hmqc")
                        && !correlation.getAtomType()
                                       .equals("H")) {
                    if (molecularConnectivity.getHsqc()
                            == null) {
                        molecularConnectivity.setHsqc(new ArrayList<>());
                    }
                    for (final int matchIndex : link.getMatch()) {
                        molecularConnectivity.getHsqc()
                                             .add(matchIndex);
                    }
                } else if (link.getExperimentType()
                               .equals("hmbc")
                        || link.getExperimentType()
                               .equals("cosy")) {
                    if (link.getExperimentType()
                            .equals("hmbc")
                            && correlation.getAtomType()
                                          .equals("H")) {
                        continue;
                    }
                    if (link.getExperimentType()
                            .equals("cosy")
                            && !correlation.getAtomType()
                                           .equals("H")) {
                        continue;
                    }
                    // ignore H atoms without any attachment to a heavy atom
                    if (correlationList.get(correlationIndex)
                                       .getAtomType()
                                       .equals("H")
                            && correlationList.get(correlationIndex)
                                              .getAttachment()
                                              .keySet()
                                              .isEmpty()) {
                        continue;
                    }
                    signal2DMap = (Map<String, Object>) link.getSignal();
                    if (signal2DMap
                            != null
                            && signal2DMap.containsKey("pathLength")) {
                        pathLengthMap = (Map<String, Object>) signal2DMap.get("pathLength");
                    } else {
                        pathLengthMap = null;
                    }
                    for (final int matchIndex : link.getMatch()) {
                        // ignore linked H atoms without any attachment to a heavy atom
                        if (correlationList.get(matchIndex)
                                           .getAtomType()
                                           .equals("H")
                                && correlationList.get(matchIndex)
                                                  .getAttachment()
                                                  .keySet()
                                                  .isEmpty()) {
                            continue;
                        }
                        if (link.getExperimentType()
                                .equals("hmbc")) {
                            if (molecularConnectivity.getHmbc()
                                    == null) {
                                molecularConnectivity.setHmbc(new HashMap<>());
                            }
                            molecularConnectivity.getHmbc()
                                                 .put(matchIndex, pathLengthMap
                                                                          == null
                                                                  ? defaultBondDistances.get("hmbc")
                                                                  : new Integer[]{(int) pathLengthMap.get("min"),
                                                                                  (int) pathLengthMap.get("max")});
                        } else {
                            if (molecularConnectivity.getCosy()
                                    == null) {
                                molecularConnectivity.setCosy(new HashMap<>());
                            }
                            molecularConnectivity.getCosy()
                                                 .put(matchIndex, pathLengthMap
                                                                          == null
                                                                  ? defaultBondDistances.get("cosy")
                                                                  : new Integer[]{(int) pathLengthMap.get("min"),
                                                                                  (int) pathLengthMap.get("max")});
                        }
                    }
                }
            }

            // set detections
            if (detections.getForbiddenNeighbors()
                          .containsKey(correlationIndex)) {
                molecularConnectivity.setForbiddenNeighbors(detections.getForbiddenNeighbors()
                                                                      .get(correlationIndex));
            }
            if (detections.getSetNeighbors()
                          .containsKey(correlationIndex)) {
                molecularConnectivity.setSetNeighbors(detections.getSetNeighbors()
                                                                .get(correlationIndex));
            }
            // fill in fixed neighbors
            if (correlation.getEquivalence()
                    == 1
                    && detections.getFixedNeighbors()
                                 .containsKey(correlationIndex)) {
                molecularConnectivity.setFixedNeighbors(new ArrayList<>());
                molecularConnectivity.getFixedNeighbors()
                                     .add(correlationIndex);
            }

            molecularConnectivityList.add(molecularConnectivity);
        }

        return molecularConnectivityList;
    }

    public static MolecularConnectivity getHeavyAtomMolecularConnectivity(
            final Map<Integer, List<MolecularConnectivity>> molecularConnectivityMap, final int protonIndex) {
        for (final Map.Entry<Integer, List<MolecularConnectivity>> entry : molecularConnectivityMap.entrySet()) {
            for (final MolecularConnectivity molecularConnectivityTemp : entry.getValue()
                                                                              .stream()
                                                                              .filter(mc -> !mc.getAtomType()
                                                                                               .equals("H"))
                                                                              .collect(Collectors.toSet())) {
                if (molecularConnectivityTemp.getHsqc()
                        != null
                        && molecularConnectivityTemp.getHsqc()
                                                    .contains(protonIndex)) {
                    return molecularConnectivityTemp;
                }
            }
        }

        return null;
    }


    //    private static void buildCombinations(final Stack<Object[]> correlationListStack,
    //                                          final List<List<Correlation>> correlationListList, final Grouping grouping) {
    //        if (correlationListStack.isEmpty()) {
    //            return;
    //        }
    //        final Object[] objects = correlationListStack.pop();
    //        final List<Correlation> correlationList = (List<Correlation>) objects[0];
    //        final int correlationIndex = (int) objects[1];
    //        final Set<String> alreadySwapped = (Set<String>) objects[2];
    //        if (correlationIndex
    //                >= correlationList.size()) {
    //            correlationListList.add(correlationList);
    //            buildCombinations(correlationListStack, correlationListList, grouping);
    //            return;
    //        }
    //        final Correlation correlation = correlationList.get(correlationIndex);
    //        if (!grouping.getGroups()
    //                     .containsKey(correlation.getAtomType())
    //                || correlation.getAtomType()
    //                              .equals("H")) {
    //            correlationListStack.push(new Object[]{correlationList, correlationIndex
    //                    + 1, alreadySwapped});
    //            buildCombinations(correlationListStack, correlationListList, grouping);
    //            return;
    //        }
    //        final List<Link> linkList = correlation.getLink()
    //                                               .stream()
    //                                               .filter(link -> link.getExperimentType()
    //                                                                   .equals("hsqc")
    //                                                       || link.getExperimentType()
    //                                                              .equals("hmqc"))
    //                                               .collect(Collectors.toList());
    //        final List<Integer> linkIndicesList = new ArrayList<>();
    //        for (final Link link : linkList) {
    //            linkIndicesList.add(correlation.getLink()
    //                                           .indexOf(link));
    //        }
    //        for (int l = 0; l
    //                < linkList.size(); l++) {
    //            final Link link = linkList.get(l);
    //            final int protonIndex = link.getMatch()
    //                                        .get(0);
    //            final int protonGroupIndex = grouping.getTransformedGroups()
    //                                                 .get("H")
    //                                                 .get(protonIndex);
    //            final List<Integer> protonGroupMemberList = new ArrayList<>(grouping.getGroups()
    //                                                                                .get("H")
    //                                                                                .get(protonGroupIndex));
    //
    //            for (final int protonGroupMemberIndex : protonGroupMemberList) {
    //                if (protonGroupMemberIndex
    //                        == protonIndex) {
    //                    continue;
    //                }
    //                System.out.println("\n\n swap at: "
    //                                           + correlationIndex
    //                                           + " -> "
    //                                           + l
    //                                           + " ---> "
    //                                           + protonIndex
    //                                           + " <-> "
    //                                           + protonGroupMemberIndex);
    //                final List<Integer> indicesToSwap = new ArrayList<>();
    //                indicesToSwap.add(protonIndex);
    //                indicesToSwap.add(protonGroupMemberIndex);
    //                indicesToSwap.sort(Integer::compare);
    //                System.out.println("-> indicesToSwap: "
    //                                           + indicesToSwap);
    //                final String swapKey = indicesToSwap.stream()
    //                                                    .map(String::valueOf)
    //                                                    .collect(Collectors.joining("_"));
    //                System.out.println("-> swapKey: "
    //                                           + swapKey);
    //                if (alreadySwapped.contains(swapKey)) {
    //                    continue;
    //                }
    //
    //
    //                // replace current link by new one in heavy atom correlation
    //                final Link newLink = Utils.cloneObject(link, Link.class);
    //                newLink.setId(Utils.generateID());
    //                newLink.setAtomType(new String[]{link.getAtomType()[0],
    //                                                 correlation.getAtomType()}); // H and current heavy atom type
    //                newLink.getMatch()
    //                       .set(0, protonGroupMemberIndex);
    //                final Correlation clonedCorrelation = Utils.cloneObject(correlation, Correlation.class);
    //                clonedCorrelation.getLink()
    //                                 .set(linkIndicesList.get(l), newLink);
    //                System.out.println(" ------> replaced proton "
    //                                           + protonIndex
    //                                           + " by "
    //                                           + protonGroupMemberIndex);
    //                // remove link from link list in current proton correlation
    //                final Correlation clonedCorrelationProton = Utils.cloneObject(correlationList.get(protonIndex),
    //                                                                              Correlation.class);
    //                System.out.println(" ------> remove heavy atom "
    //                                           + correlationIndex
    //                                           + " from "
    //                                           + protonIndex);
    //                System.out.println(" --> link array size before: "
    //                                           + clonedCorrelationProton.getLink()
    //                                                                    .size());
    //                clonedCorrelationProton.setLink(clonedCorrelationProton.getLink()
    //                                                                       .stream()
    //                                                                       .filter(linkTemp -> !(linkTemp.getExperimentID()
    //                                                                                                     .equals(link.getExperimentID())
    //                                                                               && Objects.equals(
    //                                                                               ((Map<String, Object>) linkTemp.getSignal()).get(
    //                                                                                       "id"),
    //                                                                               ((Map<String, Object>) link.getSignal()).get(
    //                                                                                       "id"))))
    //                                                                       .collect(Collectors.toList()));
    //                System.out.println(" --> link array size after:  "
    //                                           + clonedCorrelationProton.getLink()
    //                                                                    .size());
    //                // add new link in proton group member correlation
    //                final Correlation clonedCorrelationGroupMemberProton = Utils.cloneObject(
    //                        correlationList.get(protonGroupMemberIndex), Correlation.class);
    //                final Link clonedNewLink = Utils.cloneObject(newLink, Link.class);
    //                clonedNewLink.setId(Utils.generateID());
    //                clonedNewLink.setAxis("x");
    //                clonedNewLink.getMatch()
    //                             .set(0, correlationIndex);
    //                clonedCorrelationGroupMemberProton.getLink()
    //                                                  .add(clonedNewLink);
    //                System.out.println(" ------> added heavy atom "
    //                                           + correlationIndex
    //                                           + " to "
    //                                           + protonGroupMemberIndex);
    //
    //                // clone the current correlation list and set the modified correlations
    //                final List<Correlation> clonedCorrelationList = Utils.cloneList(correlationList, Correlation.class);
    //                clonedCorrelationList.set(correlationIndex, clonedCorrelation);
    //                clonedCorrelationList.set(protonGroupMemberIndex, clonedCorrelationGroupMemberProton);
    //
    //
    //                // if proton group member is attached to another heavy atom(s) then replace it there by current proton
    //                final List<Link> linkListGroupMemberProton = clonedCorrelationGroupMemberProton.getLink()
    //                                                                                               .stream()
    //                                                                                               .filter(linkTemp -> (linkTemp.getExperimentType()
    //                                                                                                                            .equals("hsqc")
    //                                                                                                       || linkTemp.getExperimentType()
    //                                                                                                                  .equals("hmqc"))
    //                                                                                                       && !linkTemp.getId()
    //                                                                                                                   .equals(clonedNewLink.getId()))
    //                                                                                               .collect(
    //                                                                                                       Collectors.toList());
    //                for (final Link linkGroupMemberProton : linkListGroupMemberProton) {
    //                    final int heavyAtomIndexGroupMemberProton = linkGroupMemberProton.getMatch()
    //                                                                                     .get(0);
    //                    final Correlation clonedLinkedHeavyAtomCorrelationOfGroupMemberProton = Utils.cloneObject(
    //                            correlationList.get(heavyAtomIndexGroupMemberProton), Correlation.class);
    //                    System.out.println(" ------> remove group member proton "
    //                                               + protonGroupMemberIndex
    //                                               + " from "
    //                                               + linkGroupMemberProton.getMatch()
    //                                                                      .get(0));
    //                    System.out.println(" --> link array size before: "
    //                                               + clonedLinkedHeavyAtomCorrelationOfGroupMemberProton.getLink()
    //                                                                                                    .size());
    //                    // remove link to previously added proton group member
    //                    clonedLinkedHeavyAtomCorrelationOfGroupMemberProton.setLink(
    //                            clonedLinkedHeavyAtomCorrelationOfGroupMemberProton.getLink()
    //                                                                               .stream()
    //                                                                               .filter(linkTemp -> !(linkTemp.getExperimentID()
    //                                                                                                             .equals(linkGroupMemberProton.getExperimentID())
    //                                                                                       && Objects.equals(
    //                                                                                       ((Map<String, Object>) linkTemp.getSignal()).get(
    //                                                                                               "id"),
    //                                                                                       ((Map<String, Object>) linkGroupMemberProton.getSignal()).get(
    //                                                                                               "id"))))
    //                                                                               .collect(Collectors.toList()));
    //                    System.out.println(" --> link array size after:  "
    //                                               + clonedLinkedHeavyAtomCorrelationOfGroupMemberProton.getLink()
    //                                                                                                    .size());
    //                    // add new link to current proton to heavy atom correlation and to current proton correlation
    //                    Link clonedLink = Utils.cloneObject(link, Link.class);
    //                    clonedLink.getMatch()
    //                              .set(0, protonIndex);
    //                    clonedLink.setId(Utils.generateID());
    //                    clonedLink.setAxis("y");
    //                    clonedLink.setAtomType(new String[]{link.getAtomType()[0],
    //                                                        clonedLinkedHeavyAtomCorrelationOfGroupMemberProton.getAtomType()});
    //                    clonedLinkedHeavyAtomCorrelationOfGroupMemberProton.getLink()
    //                                                                       .add(clonedLink);
    //                    System.out.println(" ------> added proton "
    //                                               + protonIndex
    //                                               + " to "
    //                                               + heavyAtomIndexGroupMemberProton);
    //                    clonedLink = Utils.cloneObject(link, Link.class);
    //                    clonedLink.getMatch()
    //                              .set(0, heavyAtomIndexGroupMemberProton);
    //                    clonedLink.setId(Utils.generateID());
    //                    clonedCorrelationProton.getLink()
    //                                           .add(clonedLink);
    //                    System.out.println(" ------> added heavy atom "
    //                                               + heavyAtomIndexGroupMemberProton
    //                                               + " to "
    //                                               + protonIndex);
    //
    //
    //                    // add it to new correlation list
    //                    clonedCorrelationList.set(linkGroupMemberProton.getMatch()
    //                                                                   .get(0),
    //                                              clonedLinkedHeavyAtomCorrelationOfGroupMemberProton);
    //                }
    //                clonedCorrelationList.set(protonIndex, clonedCorrelationProton);
    //
    //                alreadySwapped.add(swapKey);
    //
    //                correlationListStack.push(new Object[]{clonedCorrelationList, correlationIndex
    //                        + 1, alreadySwapped});
    //            }
    //        }
    //        correlationListStack.push(new Object[]{correlationList, correlationIndex
    //                + 1, alreadySwapped});
    //        buildCombinations(correlationListStack, correlationListList, grouping);
    //    }

    public static List<Map<Integer, List<MolecularConnectivity>>> buildMolecularConnectivityMapCombinationList(
            final List<Correlation> correlationList, final Detections detections, final Grouping grouping,
            final Map<String, Integer[]> defaultBondDistances) {
        final List<Map<Integer, List<MolecularConnectivity>>> molecularConnectivityMapCombinationList = new ArrayList<>();

        final List<MolecularConnectivity> originalMolecularConnectivityList = buildMolecularConnectivityList(
                correlationList, detections, grouping, defaultBondDistances);
        System.out.println("\n -> originalMolecularConnectivityList: \n");
        System.out.println(originalMolecularConnectivityList);

        final Map<Integer, List<MolecularConnectivity>> originalMolecularConnectivityMap = buildMolecularConnectivityMap(
                originalMolecularConnectivityList, detections, grouping, defaultBondDistances);

        molecularConnectivityMapCombinationList.add(originalMolecularConnectivityMap);

        System.out.println("\n -> originalMolecularConnectivityMap: \n");
        System.out.println(originalMolecularConnectivityMap);


        //        final Stack<Object[]> stack = new Stack<>();
        //        // experiment -> [swapKey]
        //        final Map<String, Set<String>> swapped = new HashMap<>();
        //        stack.push(new Object[]{correlationList, 0, swapped});
        //        final List<List<Correlation>> correlationListList = new ArrayList<>();
        //        buildCombinations(stack, correlationListList, grouping);
        //
        //        for (final List<Correlation> correlationListTemp : correlationListList) {
        //            molecularConnectivityMapCombinationList.add(
        //                    buildMolecularConnectivityMapWithIndices(correlationListTemp, detections, grouping,
        //                                                             defaultBondDistances));
        //        }
        //        System.out.println("\n\n -> correlationListList size: "
        //                                   + correlationListList.size()
        //                                   + "\n\n");

        return molecularConnectivityMapCombinationList;
    }
}
