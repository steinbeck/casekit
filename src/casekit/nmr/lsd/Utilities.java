package casekit.nmr.lsd;

import casekit.io.FileSystem;
import casekit.nmr.lsd.model.Detections;
import casekit.nmr.lsd.model.Grouping;
import casekit.nmr.lsd.model.MolecularConnectivity;
import casekit.nmr.model.Signal;
import casekit.nmr.model.nmrium.Correlation;
import casekit.nmr.model.nmrium.Link;
import casekit.nmr.utils.Statistics;
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

    public static boolean writeNeighborsFile(final String pathToNeighborsFile, final List<Correlation> correlationList,
                                             final Map<Integer, Object[]> indicesMap,
                                             final Map<Integer, Map<String, Map<Integer, Set<Integer>>>> neighbors) {
        final StringBuilder stringBuilder = new StringBuilder();
        Correlation correlation;
        Signal signal;
        String atomType;
        int indexInPyLSD;
        int sstrIndex = 1, sstrIndexCorrelation;
        Map<String, Map<Integer, Set<Integer>>> neighborsTemp;
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
                        for (final Map.Entry<Integer, Set<Integer>> entryPerHybridization : neighborsTemp.get(
                                                                                                                 neighborAtomType)
                                                                                                         .entrySet()) {
                            sstrIndexCorrelation = sstrIndex;
                            stringBuilder.append(
                                    casekit.nmr.lsd.inputfile.Utilities.buildSSTR(sstrIndexCorrelation, atomType,
                                                                                  correlation.getHybridization(),
                                                                                  correlation.getProtonsCount()));
                            stringBuilder.append("; ")
                                         .append(atomType)
                                         .append(" at ")
                                         .append(signal
                                                         != null
                                                 ? Statistics.roundDouble(signal.getShift(0), 2)
                                                 : "?")
                                         .append(" (")
                                         .append(indexInPyLSD)
                                         .append(")")
                                         .append("\n");
                            stringBuilder.append("ASGN S")
                                         .append(sstrIndexCorrelation)
                                         .append(" ")
                                         .append(indexInPyLSD)
                                         .append("\n");
                            sstrIndex++;

                            final List<Integer> tempList = new ArrayList<>();
                            if (entryPerHybridization.getKey()
                                    != -1) {
                                tempList.add(entryPerHybridization.getKey());
                            }
                            stringBuilder.append(
                                                 casekit.nmr.lsd.inputfile.Utilities.buildSSTR(sstrIndex, neighborAtomType, tempList,
                                                                                               new ArrayList<>(
                                                                                                       entryPerHybridization.getValue())))
                                         .append("\n");
                            stringBuilder.append("LINK S")
                                         .append(sstrIndexCorrelation)
                                         .append(" S")
                                         .append(sstrIndex)
                                         .append("\n")
                                         .append("\n");
                            sstrIndex++;
                        }
                    }
                }
            }
        }

        System.out.println(stringBuilder);


        return !stringBuilder.toString()
                             .isEmpty()
                && FileSystem.writeFile(pathToNeighborsFile, stringBuilder.toString());
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

    public static boolean hasMatch(final Correlation correlation1, final Correlation correlation2,
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

    public static Map<String, Map<Integer, Set<Integer>>> findGroups(final List<Correlation> correlationList,
                                                                     final Map<String, Double> tolerances) {
        // cluster group index -> list of correlation index pair
        final Map<String, Map<Integer, Set<Integer>>> groups = new HashMap<>();
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
            for (final Map.Entry<Integer, Set<Integer>> groupEntry : groups.get(correlation.getAtomType())
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
                      .put(groupIndex, new HashSet<>());
                groups.get(correlation.getAtomType())
                      .get(groupIndex)
                      .add(i);
                inserted.add(i);
                groupIndex++;
            }
        }

        return groups;
    }

    public static Map<String, Map<Integer, Integer>> transformGroups(
            final Map<String, Map<Integer, Set<Integer>>> groups) {
        final Map<String, Map<Integer, Integer>> transformedGroups = new HashMap<>();
        for (final Map.Entry<String, Map<Integer, Set<Integer>>> atomTypeEntry : groups.entrySet()) {
            transformedGroups.put(atomTypeEntry.getKey(), new HashMap<>());
            for (final Map.Entry<Integer, Set<Integer>> groupEntry : atomTypeEntry.getValue()
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
        final Map<String, Map<Integer, Set<Integer>>> groups = findGroups(correlationList, tolerances);

        return new Grouping(tolerances, groups, transformGroups(groups));
    }

    private static Set<Integer> getProtonCounts(final List<Correlation> correlationList, final int index) {
        final Correlation correlation = correlationList.get(index);
        if (correlation.getProtonsCount()
                != null
                && !correlation.getProtonsCount()
                               .isEmpty()) {
            // if protonCounts is already given
            return new HashSet<>(correlation.getProtonsCount());
        }
        final Set<Integer> protonCounts = new HashSet<>();
        for (int i = 0; i
                < Constants.defaultProtonsCountPerValencyMap.get(
                Constants.defaultAtomLabelMap.get(correlation.getAtomType())).length; i++) {
            protonCounts.add(Constants.defaultProtonsCountPerValencyMap.get(
                    Constants.defaultAtomLabelMap.get(correlation.getAtomType()))[i]);
        }

        return protonCounts;
    }

    private static Set<Integer> getHybridizations(final List<Correlation> correlationList, final int index,
                                                  final Map<Integer, List<Integer>> detectedHybridizations) {
        final Correlation correlation = correlationList.get(index);
        Set<Integer> hybridizations = new HashSet<>();
        if (correlation.getHybridization()
                != null
                && !correlation.getHybridization()
                               .isEmpty()) {
            // if hybridization is already given
            return new HashSet<>(correlation.getHybridization());
        } else {
            // if hybridization is not given then use the detected ones
            if (detectedHybridizations.containsKey(index)) {
                hybridizations = new HashSet<>(detectedHybridizations.get(index));
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

    private static Map<Integer, Integer[]> buildIndicesMap(final List<Correlation> correlationList) {
        // index in correlation data -> [indices in PyLSD file...]
        final Map<Integer, Integer[]> indicesMap = new HashMap<>();
        // init element indices within correlations with same order as in correlation data input
        int heavyAtomIndexInPyLSDFile = 1;
        int protonIndexInPyLSDFile = 1;
        int protonsToInsert, protonsCount;
        Correlation correlation;
        for (int i = 0; i
                < correlationList.size(); i++) {
            correlation = correlationList.get(i);
            // set entry for each correlation with consideration of equivalences
            if (correlation.getAtomType()
                           .equals("H")) {
                protonsToInsert = 0;
                for (final Link link : correlation.getLink()) {
                    if (link.getExperimentType()
                            .equals("hsqc")
                            || link.getExperimentType()
                                   .equals("hmqc")) {
                        for (final int matchIndex : link.getMatch()) {
                            protonsCount = correlationList.get(matchIndex)
                                                          .getProtonsCount()
                                                          .get(0);
                            protonsToInsert += (correlation.getEquivalence()
                                    / (double) protonsCount)
                                    * correlationList.get(matchIndex)
                                                     .getAttachment()
                                                     .get("H")
                                                     .size();
                        }
                    }
                }
                indicesMap.put(i, new Integer[protonsToInsert]);
                for (int j = 0; j
                        < protonsToInsert; j++) {
                    indicesMap.get(i)[j] = protonIndexInPyLSDFile;
                    protonIndexInPyLSDFile++;
                }
            } else {
                indicesMap.put(i, new Integer[correlation.getEquivalence()]);
                for (int j = 0; j
                        < correlation.getEquivalence(); j++) {
                    indicesMap.get(i)[j] = heavyAtomIndexInPyLSDFile;
                    heavyAtomIndexInPyLSDFile++;
                }
            }
        }

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

    private static Set<Integer> buildGroupMembers(final Map<Integer, Integer[]> indicesMap,
                                                  final List<Correlation> correlationList, final int correlationIndex,
                                                  final Grouping grouping) {
        final Correlation correlation = correlationList.get(correlationIndex);
        // add possible indices from grouping
        final int groupIndex;
        final Set<Integer> possibilities = new HashSet<>();
        if (grouping.getTransformedGroups()
                    .containsKey(correlation.getAtomType())) {
            groupIndex = grouping.getTransformedGroups()
                                 .get(correlation.getAtomType())
                                 .get(correlationIndex);
            for (final int groupCorrelationIndex : grouping.getGroups()
                                                           .get(correlation.getAtomType())
                                                           .get(groupIndex)) {
                // add equivalence indices of group members
                for (int k = 0; k
                        < indicesMap.get(groupCorrelationIndex).length; k++) {
                    possibilities.add(indicesMap.get(groupCorrelationIndex)[k]);
                }
            }
        } else {
            // add for equivalences only
            for (int k = 0; k
                    < indicesMap.get(correlationIndex).length; k++) {
                possibilities.add(indicesMap.get(correlationIndex)[k]);
            }
        }

        return possibilities;
    }

    private static void addMolecularConnectivity(
            final Map<Integer, List<MolecularConnectivity>> molecularConnectivityMap,
            final Map<Integer, Integer[]> indicesMap, final List<Correlation> correlationList,
            final int correlationIndex, final int index, final Detections detections, final Grouping grouping) {
        final Correlation correlation = correlationList.get(correlationIndex);
        molecularConnectivityMap.putIfAbsent(correlationIndex, new ArrayList<>());
        if (molecularConnectivityMap.get(correlationIndex)
                                    .stream()
                                    .noneMatch(molecularConnectivityTemp -> molecularConnectivityTemp.getIndex()
                                            == index)) {
            final MolecularConnectivity molecularConnectivity = new MolecularConnectivity();
            molecularConnectivity.setIndex(index);
            molecularConnectivity.setAtomType(correlation.getAtomType());
            molecularConnectivity.setSignal(Utils.extractSignalFromCorrelation(correlation));

            if (!correlationList.get(correlationIndex)
                                .getAtomType()
                                .equals("H")) {
                molecularConnectivity.setProtonCounts(getProtonCounts(correlationList, correlationIndex));
                molecularConnectivity.setHybridizations(
                        getHybridizations(correlationList, correlationIndex, detections.getDetectedHybridizations()));
            }
            molecularConnectivity.setGroupMembers(
                    buildGroupMembers(indicesMap, correlationList, correlationIndex, grouping));
            molecularConnectivityMap.get(correlationIndex)
                                    .add(molecularConnectivity);
        }
    }

    public static Map<Integer, List<MolecularConnectivity>> buildMolecularConnectivityMap(
            final List<Correlation> correlationList, final Detections detections, final Grouping grouping,
            final Map<String, Integer[]> defaultBondDistances) {

        final Map<Integer, Integer[]> indicesMap = buildIndicesMap(correlationList);
        // correlation index -> [MolecularConnectivity]
        final Map<Integer, List<MolecularConnectivity>> molecularConnectivityMap = new HashMap<>();
        Correlation correlation, correlation2;
        int index, protonIndex;
        Map<String, Object> signal2DMap;
        Map<String, Object> pathLengthMap;
        MolecularConnectivity molecularConnectivity;
        for (int correlationIndex = 0; correlationIndex
                < correlationList.size(); correlationIndex++) {
            correlation = correlationList.get(correlationIndex);
            for (int k = 0; k
                    < indicesMap.get(correlationIndex).length; k++) {
                index = indicesMap.get(correlationIndex)[k];
                addMolecularConnectivity(molecularConnectivityMap, indicesMap, correlationList, correlationIndex, index,
                                         detections, grouping);
                molecularConnectivity = findMolecularConnectivityByIndex(molecularConnectivityMap,
                                                                         correlation.getAtomType(), false, index);
                for (final Link link : correlation.getLink()) {
                    if (link.getExperimentType()
                            .equals("hsqc")
                            || link.getExperimentType()
                                   .equals("hmqc")
                            && !correlation.getAtomType()
                                           .equals("H")) {
                        for (final int matchIndex : link.getMatch()) {
                            protonIndex = indicesMap.get(matchIndex)[k];
                            if (molecularConnectivity.getHsqc()
                                    == null) {
                                molecularConnectivity.setHsqc(new HashSet<>());
                            }
                            molecularConnectivity.getHsqc()
                                                 .add(protonIndex);
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
                                for (int l = 0; l
                                        < indicesMap.get(matchIndex).length; l++) {
                                    protonIndex = indicesMap.get(matchIndex)[l];
                                    if (molecularConnectivity.getHmbc()
                                            == null) {
                                        molecularConnectivity.setHmbc(new HashMap<>());
                                    }
                                    molecularConnectivity.getHmbc()
                                                         .put(protonIndex, pathLengthMap
                                                                                   == null
                                                                           ? defaultBondDistances.get("hmbc")
                                                                           : new Integer[]{(int) pathLengthMap.get(
                                                                                   "min"), (int) pathLengthMap.get(
                                                                                   "max")});
                                }
                            } else {
                                if (k
                                        < indicesMap.get(matchIndex).length) {
                                    protonIndex = indicesMap.get(matchIndex)[k];
                                    if (molecularConnectivity.getCosy()
                                            == null) {
                                        molecularConnectivity.setCosy(new HashMap<>());
                                    }
                                    molecularConnectivity.getCosy()
                                                         .put(protonIndex, pathLengthMap
                                                                                   == null
                                                                           ? defaultBondDistances.get("cosy")
                                                                           : new Integer[]{(int) pathLengthMap.get(
                                                                                   "min"), (int) pathLengthMap.get(
                                                                                   "max")});
                                }
                            }
                        }
                    }
                }
            }
            // set detections
            if (!correlation.getAtomType()
                            .equals("H")
                    && !correlation.isPseudo()) {
                for (final MolecularConnectivity molecularConnectivityTemp : molecularConnectivityMap.get(
                        correlationIndex)) {
                    molecularConnectivityTemp.setForbiddenNeighbors(detections.getForbiddenNeighbors()
                                                                              .get(correlationIndex));
                    molecularConnectivityTemp.setSetNeighbors(detections.getSetNeighbors()
                                                                        .get(correlationIndex));
                }
            }
            // fill in fixed neighbors
            if (correlation.getEquivalence()
                    == 1
                    && detections.getFixedNeighbors()
                                 .containsKey(correlationIndex)) {
                for (final int correlationIndex2 : detections.getFixedNeighbors()
                                                             .get(correlationIndex)) {
                    correlation2 = correlationList.get(correlationIndex2);
                    // use fixed neighbor information of atoms without equivalence equals 1 only
                    if (correlation2.getEquivalence()
                            > 1) {
                        continue;
                    }
                    index = indicesMap.get(correlationIndex)[0];
                    molecularConnectivity = findMolecularConnectivityByIndex(molecularConnectivityMap,
                                                                             correlation.getAtomType(), false, index);
                    if (molecularConnectivity.getFixedNeighbors()
                            == null) {
                        molecularConnectivity.setFixedNeighbors(new HashSet<>());
                    }
                    molecularConnectivity.getFixedNeighbors()
                                         .add(indicesMap.get(correlationIndex2)[0]);
                }
            }
        }
        // filter out HMBC or COSY correlation to itself
        for (final int correlationIndex : molecularConnectivityMap.keySet()) {
            for (final MolecularConnectivity molecularConnectivityTemp : molecularConnectivityMap.get(
                    correlationIndex)) {
                if (molecularConnectivityTemp.getHsqc()
                        != null) {
                    if (molecularConnectivityTemp.getHmbc()
                            != null) {
                        for (final int hmbcIndexInPyLSD : molecularConnectivityTemp.getHsqc()) {
                            molecularConnectivityTemp.getHmbc()
                                                     .remove(hmbcIndexInPyLSD);
                        }
                    }
                    if (molecularConnectivityTemp.getCosy()
                            != null) {
                        for (final int cosyIndexInPyLSD : molecularConnectivityTemp.getHsqc()) {
                            molecularConnectivityTemp.getCosy()
                                                     .remove(cosyIndexInPyLSD);
                        }
                    }
                }
            }
        }

        return molecularConnectivityMap;
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
}
