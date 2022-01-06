package casekit.nmr.lsd.inputfile;

import casekit.nmr.lsd.Constants;
import casekit.nmr.lsd.model.Detections;
import casekit.nmr.lsd.model.ElucidationOptions;
import casekit.nmr.model.Signal;
import casekit.nmr.model.nmrium.Correlation;
import casekit.nmr.model.nmrium.Correlations;
import casekit.nmr.model.nmrium.Link;
import casekit.nmr.utils.Statistics;
import casekit.nmr.utils.Utils;

import java.text.SimpleDateFormat;
import java.util.*;

public class PyLSDInputFileBuilder {

    private static String buildHeader() {
        final StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append("; PyLSD input file created by casekit (https://github.com/michaelwenk/casekit)\n");
        final SimpleDateFormat formatter = new SimpleDateFormat("yyyy-MM-dd 'at' HH:mm:ss z");
        final Date date = new Date(System.currentTimeMillis());
        stringBuilder.append("; ")
                     .append(formatter.format(date));

        return stringBuilder.toString();
    }

    private static String buildFORM(final String mf, final Map<String, Integer> elementCounts) {
        final StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append("; Molecular Formula: ")
                     .append(mf)
                     .append("\n");
        stringBuilder.append("FORM ");
        elementCounts.forEach((elem, count) -> stringBuilder.append(elem)
                                                            .append(" ")
                                                            .append(count)
                                                            .append(" "));

        return stringBuilder.toString();
    }

    private static String buildPIEC() {
        return "PIEC 1";
    }

    private static String buildELIM(final int elimP1, final int elimP2) {
        return "ELIM "
                + elimP1
                + " "
                + elimP2;
    }

    private static Map<Integer, Object[]> buildIndicesMap(final List<Correlation> correlationList) {
        // index in correlation data -> [atom type, indices in PyLSD file...]
        final Map<Integer, Object[]> indicesMap = new HashMap<>();
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
                indicesMap.put(i, new Object[1
                        + protonsToInsert]);
                indicesMap.get(i)[0] = correlation.getAtomType();
                for (int j = 1; j
                        <= protonsToInsert; j++) {
                    indicesMap.get(i)[j] = protonIndexInPyLSDFile;
                    protonIndexInPyLSDFile++;
                }
            } else {
                indicesMap.put(i, new Object[1
                        + correlation.getEquivalence()]);
                indicesMap.get(i)[0] = correlation.getAtomType();
                for (int j = 1; j
                        <= correlation.getEquivalence(); j++) {
                    indicesMap.get(i)[j] = heavyAtomIndexInPyLSDFile;
                    heavyAtomIndexInPyLSDFile++;
                }
            }
        }
        //        System.out.println("\n\n");
        //        for (final Map.Entry<Integer, Object[]> entry : indicesMap.entrySet()) {
        //            System.out.println(entry.getKey()
        //                                       + ": "
        //                                       + entry.getValue()[0]);
        //            for (int i = 1; i
        //                    < entry.getValue().length; i++) {
        //                System.out.println(entry.getValue()[i]);
        //            }
        //        }
        //        System.out.println("\n\n");

        return indicesMap;
    }

    private static String buildMULT(final List<Correlation> correlationList, final int index,
                                    final Map<Integer, Object[]> indicesMap,
                                    final Map<Integer, List<Integer>> detectedHybridizations) {
        final Correlation correlation = correlationList.get(index);
        if (correlation.getAtomType()
                       .equals("H")) {
            return null;
        }
        final StringBuilder stringBuilder = new StringBuilder();
        List<Integer> hybridizations = new ArrayList<>();
        final StringBuilder hybridizationStringBuilder;
        final StringBuilder attachedProtonsCountStringBuilder;

        if (correlation.getHybridization()
                != null
                && !correlation.getHybridization()
                               .isEmpty()) {
            // if hybridization is already given
            hybridizations.addAll(correlation.getHybridization());
        } else {
            // if hybridization is not given then use the detected ones
            if (detectedHybridizations.containsKey(index)) {
                hybridizations = detectedHybridizations.get(index);
            }
            if (hybridizations.isEmpty()
                    && correlation.getAtomType()
                                  .equals("C")
                    && !correlation.getProtonsCount()
                                   .isEmpty()
                    && correlation.getProtonsCount()
                                  .get(0)
                    >= 2) {
                // a carbon with at least two protons can only be SP2 or SP3
                hybridizations.add(2);
                hybridizations.add(3);
            }
        }
        if (hybridizations.isEmpty()) {
            hybridizationStringBuilder = new StringBuilder();
            if (Constants.defaultHybridizationMap.get(correlation.getAtomType()).length
                    > 1) {
                hybridizationStringBuilder.append("(");
            }
            for (int i = 0; i
                    < Constants.defaultHybridizationMap.get(correlation.getAtomType()).length; i++) {
                hybridizationStringBuilder.append(Constants.defaultHybridizationMap.get(correlation.getAtomType())[i]);
                if (i
                        < Constants.defaultHybridizationMap.get(correlation.getAtomType()).length
                        - 1) {
                    hybridizationStringBuilder.append(" ");
                }
            }
            if (Constants.defaultHybridizationMap.get(correlation.getAtomType()).length
                    > 1) {
                hybridizationStringBuilder.append(")");
            }
        } else {
            hybridizationStringBuilder = new StringBuilder();
            if (hybridizations.size()
                    > 1) {
                hybridizationStringBuilder.append("(");
            }
            for (int k = 0; k
                    < hybridizations.size(); k++) {
                hybridizationStringBuilder.append(hybridizations.get(k));
                if (k
                        < hybridizations.size()
                        - 1) {
                    hybridizationStringBuilder.append(" ");
                }
            }
            if (hybridizations.size()
                    > 1) {
                hybridizationStringBuilder.append(")");
            }
        }
        // set attached protons count
        attachedProtonsCountStringBuilder = new StringBuilder();
        // if protons count is given
        if (correlation.getProtonsCount()
                != null
                && !correlation.getProtonsCount()
                               .isEmpty()) {
            if (correlation.getProtonsCount()
                           .size()
                    == 1) {
                attachedProtonsCountStringBuilder.append(correlation.getProtonsCount()
                                                                    .get(0));
            } else {
                attachedProtonsCountStringBuilder.append("(");
                for (final int protonsCount : correlation.getProtonsCount()) {
                    attachedProtonsCountStringBuilder.append(protonsCount)
                                                     .append(" ");
                }
                attachedProtonsCountStringBuilder.deleteCharAt(attachedProtonsCountStringBuilder.length()
                                                                       - 1);
                attachedProtonsCountStringBuilder.append(")");
            }
        } else { // if protons count is not given then set it to default value
            if (Constants.defaultProtonsCountPerValencyMap.get(
                    Constants.defaultAtomLabelMap.get(correlation.getAtomType())).length
                    > 1) {
                attachedProtonsCountStringBuilder.append("(");
            }
            for (int i = 0; i
                    < Constants.defaultProtonsCountPerValencyMap.get(
                    Constants.defaultAtomLabelMap.get(correlation.getAtomType())).length; i++) {
                attachedProtonsCountStringBuilder.append(Constants.defaultProtonsCountPerValencyMap.get(
                        Constants.defaultAtomLabelMap.get(correlation.getAtomType()))[i]);
                if (i
                        < Constants.defaultProtonsCountPerValencyMap.get(
                        Constants.defaultAtomLabelMap.get(correlation.getAtomType())).length
                        - 1) {
                    attachedProtonsCountStringBuilder.append(" ");
                }
            }
            if (Constants.defaultProtonsCountPerValencyMap.get(
                    Constants.defaultAtomLabelMap.get(correlation.getAtomType())).length
                    > 1) {
                attachedProtonsCountStringBuilder.append(")");
            }
        }
        for (int j = 1; j
                < indicesMap.get(index).length; j++) {
            stringBuilder.append("MULT ")
                         .append(indicesMap.get(index)[j])
                         .append(" ")
                         .append(correlation.getAtomType())
                         .append(" ")
                         .append(hybridizationStringBuilder)
                         .append(" ")
                         .append(attachedProtonsCountStringBuilder);
            if (!correlation.isPseudo()) {
                stringBuilder.append("; ")
                             .append(buildShiftString(correlationList, correlation));
            }
            if (j
                    >= 2) {
                stringBuilder.append("; equivalent to ")
                             .append(indicesMap.get(index)[1]);
            }
            stringBuilder.append("\n");
        }

        return stringBuilder.toString();
    }

    private static String buildShiftString(final List<Correlation> correlationList, final Correlation correlation) {

        final Signal signal = Utils.extractSignalFromCorrelation(correlation);
        if (signal
                == null) {
            return "?";
        }

        String heavyAtomShiftString = "";
        if (correlation.getAtomType()
                       .equals("H")) {
            final ArrayList<String> bondHeavyAtomTypes = new ArrayList<>(correlation.getAttachment()
                                                                                    .keySet());
            if (!bondHeavyAtomTypes.isEmpty()) {
                final Optional<Integer> firstOptional = correlation.getAttachment()
                                                                   .get(bondHeavyAtomTypes.get(0))
                                                                   .stream()
                                                                   .findFirst();
                if (firstOptional.isPresent()) {
                    heavyAtomShiftString = " ("
                            + buildShiftString(correlationList, correlationList.get(firstOptional.get()))
                            + ")";
                }
            }

        }

        return correlation.isPseudo()
               ? "?"
               : Statistics.roundDouble(signal.getShift(0), 3)
                       + heavyAtomShiftString;
    }

    private static String buildShiftsComment(final List<Correlation> correlationList, final Correlation correlation1,
                                             final Correlation correlation2) {
        return "; "
                + correlation1.getAtomType()
                + ": "
                + buildShiftString(correlationList, correlation1)
                + " -> "
                + correlation2.getAtomType()
                + ": "
                + buildShiftString(correlationList, correlation2);
    }

    private static String buildHSQC(final List<Correlation> correlationList, final int index,
                                    final Map<Integer, Object[]> indicesMap) {
        final Correlation correlation = correlationList.get(index);
        if (correlation.getAtomType()
                       .equals("H")) {
            return null;
        }
        final StringBuilder stringBuilder = new StringBuilder();
        for (final Link link : correlation.getLink()) {
            if (link.getExperimentType()
                    .equals("hsqc")
                    || link.getExperimentType()
                           .equals("hmqc")) {
                for (final int matchIndex : link.getMatch()) {
                    // for each equivalence of heavy atom and attached protons
                    for (int k = 1; k
                            < indicesMap.get(index).length; k++) {
                        stringBuilder.append("HSQC ")
                                     .append(indicesMap.get(index)[k])
                                     .append(" ")
                                     .append(indicesMap.get(matchIndex)[k])
                                     .append(buildShiftsComment(correlationList, correlation,
                                                                correlationList.get(matchIndex)))
                                     .append("\n");
                    }
                }
            }
        }

        return stringBuilder.toString();
    }

    private static String buildPossibilitiesString(final Map<Integer, Object[]> indicesMap, final int index) {
        final StringBuilder possibilitiesStringBuilder = new StringBuilder();
        if (indicesMap.get(index).length
                > 2) {
            possibilitiesStringBuilder.append("(");
        }
        for (int k = 1; k
                < indicesMap.get(index).length; k++) {
            possibilitiesStringBuilder.append((int) indicesMap.get(index)[k]);
            if (k
                    < indicesMap.get(index).length
                    - 1) {
                possibilitiesStringBuilder.append(" ");
            }
        }
        if (indicesMap.get(index).length
                > 2) {
            possibilitiesStringBuilder.append(")");
        }

        return possibilitiesStringBuilder.toString();
    }

    private static String buildHMBC(final List<Correlation> correlationList, final int index,
                                    final Map<Integer, Object[]> indicesMap) {
        final Correlation correlation = correlationList.get(index);
        if (correlation.getAtomType()
                       .equals("H")) {
            return null;
        }

        //        final Set<Integer> group = new HashSet<>();
        //        if (grouping.getTransformedGroups()
        //                    .containsKey(correlation.getAtomType())
        //                && grouping.getTransformedGroups()
        //                           .get(correlation.getAtomType())
        //                           .containsKey(index)) {
        //            final int groupIndex = grouping.getTransformedGroups()
        //                                           .get(correlation.getAtomType())
        //                                           .get(index);
        //            group.addAll(grouping.getGroups()
        //                                 .get(correlation.getAtomType())
        //                                 .get(groupIndex));
        //            System.out.println("\nindex: "
        //                                       + index
        //                                       + " -> groupIndex: "
        //                                       + groupIndex
        //                                       + " -> group: "
        //                                       + group);
        //        }

        final String defaultBondDistanceString = 2
                + " "
                + 3;
        String bondDistanceString;
        Map<String, Object> signal2DMap, pathLengthMap;
        final Set<String> uniqueSet = new LinkedHashSet<>(); // in case of same content exists multiple times
        for (final Link link : correlation.getLink()) {
            if (link.getExperimentType()
                    .equals("hmbc")) {
                for (final int matchIndex : link.getMatch()) {
                    //                    for (int k = 1; k
                    //                            < indicesMap.get(index).length; k++) {
                    for (int l = 1; l
                            < indicesMap.get(matchIndex).length; l++) {
                        //                            // only add an HMBC correlation if there is no direct link via HSQC and the equivalence index is not equal
                        //                            if (!(correlationList.get(matchIndex)
                        //                                                 .getAttachment()
                        //                                                 .containsKey(correlation.getAtomType())
                        //                                    && correlationList.get(matchIndex)
                        //                                                      .getAttachment()
                        //                                                      .get(correlation.getAtomType())
                        //                                                      .contains(index))) {
                        bondDistanceString = null;
                        signal2DMap = (Map<String, Object>) link.getSignal();
                        if (signal2DMap
                                != null
                                && signal2DMap.containsKey("pathLength")) {
                            pathLengthMap = (Map<String, Object>) signal2DMap.get("pathLength");
                            bondDistanceString = pathLengthMap.get("min")
                                    + " "
                                    + pathLengthMap.get("max");
                        }
                        uniqueSet.add(buildPossibilitiesString(indicesMap, index)
                                              //indicesMap.get(index)[k]
                                              + " "
                                              + indicesMap.get(matchIndex)[l]
                                              + " "
                                              + (bondDistanceString
                                                         != null
                                                 ? bondDistanceString
                                                 : defaultBondDistanceString)
                                              + buildShiftsComment(correlationList, correlation,
                                                                   correlationList.get(matchIndex)));
                        //                    }
                        //                        }
                    }
                }
            }
        }

        return uniqueSet.stream()
                        .map(str -> "HMBC "
                                + str
                                + "\n")
                        .reduce("", (strAll, str) -> strAll
                                + str);
    }

    private static String buildCOSY(final List<Correlation> correlationList, final int index,
                                    final Map<Integer, Object[]> indicesMap) {
        final Correlation correlation = correlationList.get(index);
        if (!correlation.getAtomType()
                        .equals("H")) {
            return null;
        }
        final String defaultBondDistanceString = 3
                + " "
                + 4;
        String bondDistanceString;
        Map<String, Object> signal2DMap, pathLengthMap;
        final Set<String> uniqueSet = new LinkedHashSet<>(); // in case of same content exists multiple times
        for (final Link link : correlation.getLink()) {
            if (link.getExperimentType()
                    .equals("cosy")) {
                for (final int matchIndex : link.getMatch()) {
                    //                    // only add a COSY entry if it is not from same signal
                    //                    if (!correlationList.get(matchIndex)
                    //                                        .getId()
                    //                                        .equals(correlation.getId())) {
                    for (int l = 1; l
                            < indicesMap.get(matchIndex).length; l++) {
                        //                        // only allow COSY values between possible equivalent protons and only one another non-equivalent proton
                        //                        if (indicesMap.get(matchIndex).length
                        //                                == 2) {
                        bondDistanceString = null;
                        signal2DMap = (Map<String, Object>) link.getSignal();
                        if (signal2DMap
                                != null
                                && signal2DMap.containsKey("pathLength")) {
                            pathLengthMap = (Map<String, Object>) signal2DMap.get("pathLength");
                            bondDistanceString = pathLengthMap.get("min")
                                    + " "
                                    + pathLengthMap.get("max");
                        }
                        uniqueSet.add(buildPossibilitiesString(indicesMap, index)
                                              //indicesMap.get(index)[k]
                                              + " "
                                              + indicesMap.get(matchIndex)[l]
                                              + " "
                                              + (bondDistanceString
                                                         != null
                                                 ? bondDistanceString
                                                 : defaultBondDistanceString)
                                              + buildShiftsComment(correlationList, correlation,
                                                                   correlationList.get(matchIndex)));
                        //                        }
                    }
                    //                    }
                }
            }
        }

        return uniqueSet.stream()
                        .map(str -> "COSY "
                                + str
                                + "\n")
                        .reduce("", (strAll, str) -> strAll
                                + str);
    }

    private static String buildSHIX(final Correlation correlation, final int index,
                                    final Map<Integer, Object[]> indicesMap) {
        if (correlation.getAtomType()
                       .equals("H")
                || correlation.isPseudo()) {
            return null;
        }
        final Signal signal = Utils.extractSignalFromCorrelation(correlation);
        if (signal
                == null) {
            return null;
        }
        final StringBuilder stringBuilder = new StringBuilder();
        for (int k = 1; k
                < indicesMap.get(index).length; k++) {
            stringBuilder.append("SHIX ")
                         .append(indicesMap.get(index)[k])
                         .append(" ")
                         .append(Statistics.roundDouble(signal.getShift(0), 3))
                         .append("\n");
        }

        return stringBuilder.toString();
    }

    private static String buildSHIH(final Correlation correlation, final int index,
                                    final Map<Integer, Object[]> indicesMap) {
        if (!correlation.getAtomType()
                        .equals("H")
                || correlation.isPseudo()) {
            return null;
        }
        final Signal signal = Utils.extractSignalFromCorrelation(correlation);
        if (signal
                == null) {
            return null;
        }
        final StringBuilder stringBuilder = new StringBuilder();
        // only consider protons which are attached via HSQC/HMQC (pseudo and real links)
        for (final Link link : correlation.getLink()) {
            if ((link.getExperimentType()
                     .equals("hsqc")
                    || link.getExperimentType()
                           .equals("hmqc"))
                    && !link.getMatch()
                            .isEmpty()) { // && !link.isPseudo()
                for (int k = 1; k
                        < indicesMap.get(index).length; k++) {
                    stringBuilder.append("SHIH ")
                                 .append(indicesMap.get(index)[k])
                                 .append(" ")
                                 .append(Statistics.roundDouble(signal.getShift(0), 5))
                                 .append("\n");
                }
            }
        }

        return stringBuilder.toString();
    }

    private static String buildLISTsAndPROPs(final List<Correlation> correlationList,
                                             final Map<Integer, Object[]> indicesMap,
                                             final Map<String, Integer> elementCounts,
                                             final Map<Integer, Map<String, Map<Integer, Set<Integer>>>> forbiddenNeighbors,
                                             final Map<Integer, Map<String, Map<Integer, Set<Integer>>>> setNeighbors,
                                             final boolean allowHeteroHeteroBonds) {
        final StringBuilder stringBuilder = new StringBuilder();
        // list key -> [list name, size]
        final Map<String, Object[]> listMap = new HashMap<>();

        // LIST and PROP for hetero hetero bonds to disallow
        if (!allowHeteroHeteroBonds) {
            LISTAndPROPUtilities.insertNoHeteroHeteroBonds(stringBuilder, listMap);
        }
        // insert ELEM for each heavy atom type in MF
        LISTAndPROPUtilities.insertELEM(stringBuilder, listMap, elementCounts.keySet());
        // insert list combinations of carbon and hybridization states
        LISTAndPROPUtilities.insertHeavyAtomCombinationLISTs(stringBuilder, listMap, correlationList, indicesMap);
        // insert forbidden connection lists and properties
        LISTAndPROPUtilities.insertConnectionLISTsAndPROPs(stringBuilder, listMap, correlationList, indicesMap,
                                                           forbiddenNeighbors, "forbid");
        // insert set connection lists and properties
        LISTAndPROPUtilities.insertConnectionLISTsAndPROPs(stringBuilder, listMap, correlationList, indicesMap,
                                                           setNeighbors, "allow");

        return stringBuilder.toString();
    }

    private static String buildDEFFs(final String[] filterPaths, final String[] pathsToNeighborsFiles) {
        final StringBuilder stringBuilder = new StringBuilder();
        // DEFF -> add filters
        final Map<String, String> filters = new LinkedHashMap<>();
        int counter = 1;
        for (final String filterPath : filterPaths) {
            filters.put("F"
                                + counter, filterPath);
            counter++;
        }
        for (final String pathToNeighborsFiles : pathsToNeighborsFiles) {
            filters.put("F"
                                + counter, pathToNeighborsFiles);
            counter++;
        }

        if (!filters.isEmpty()) {
            stringBuilder.append("; externally defined filters\n");
            filters.forEach((label, filePath) -> stringBuilder.append("DEFF ")
                                                              .append(label)
                                                              .append(" \"")
                                                              .append(filePath)
                                                              .append("\"\n"));
            stringBuilder.append("\n");
        }

        return stringBuilder.toString();
    }

    private static String buildFEXP(final Map<String, Boolean> fexpMap) {
        final StringBuilder stringBuilder = new StringBuilder();

        if (!fexpMap.isEmpty()) {
            stringBuilder.append("FEXP \"");
            int counter = 0;
            for (final String label : fexpMap.keySet()) {
                if (!fexpMap.get(label)) {
                    stringBuilder.append("NOT ");
                }
                stringBuilder.append(label);
                if (counter
                        < fexpMap.keySet()
                                 .size()
                        - 1) {
                    stringBuilder.append(" and ");
                }
                counter++;
            }
            stringBuilder.append("\"\n");
        }

        return stringBuilder.toString();
    }

    private static String buildDEFFsAndFEXP(final List<Correlation> correlationList,
                                            final Map<Integer, Object[]> indicesMap,
                                            final ElucidationOptions elucidationOptions,
                                            final Map<Integer, Map<String, Map<Integer, Set<Integer>>>> forbiddenNeighbors,
                                            final Map<Integer, Map<String, Map<Integer, Set<Integer>>>> setNeighbors) {
        final StringBuilder stringBuilder = new StringBuilder();
        final Map<String, Boolean> fexpMap = new HashMap<>();
        for (int i = 0; i
                < elucidationOptions.getFilterPaths().length; i++) {
            fexpMap.put("F"
                                + (i
                    + 1), false);
        }
        //        // build and write neighbors files
        //        final List<String> pathsToNeighborsFilesToUse = new ArrayList<>();
        //        if (Utilities.writeNeighborsFile(elucidationOptions.getPathsToNeighborsFiles()[0], correlationList, indicesMap,
        //                                         forbiddenNeighbors)) {
        //            fexpMap.put("F"
        //                                + (fexpMap.size()
        //                    + 1), false);
        //            pathsToNeighborsFilesToUse.add(elucidationOptions.getPathsToNeighborsFiles()[0]);
        //        }
        //        if (Utilities.writeNeighborsFile(elucidationOptions.getPathsToNeighborsFiles()[1], correlationList, indicesMap,
        //                                         setNeighbors)) {
        //            fexpMap.put("F"
        //                                + (fexpMap.size()
        //                    + 1), true);
        //            pathsToNeighborsFilesToUse.add(elucidationOptions.getPathsToNeighborsFiles()[1]);
        //        }
        // build DEFFs
        stringBuilder.append(
                             //                             buildDEFFs(elucidationOptions.getFilterPaths(), pathsToNeighborsFilesToUse.toArray(String[]::new)))
                             buildDEFFs(elucidationOptions.getFilterPaths(), new String[]{}))
                     .append("\n");
        // build FEXP
        stringBuilder.append(buildFEXP(fexpMap))
                     .append("\n");

        return stringBuilder.toString();
    }

    private static String buildBONDByFixedNeighbors(final List<Correlation> correlationList,
                                                    final Map<Integer, Object[]> indicesMap,
                                                    final Map<Integer, Set<Integer>> fixedNeighbors) {
        final StringBuilder stringBuilder = new StringBuilder();

        final Set<String> uniqueSet = new HashSet<>();
        int correlationIndex1;
        Correlation correlation1, correlation2;
        for (final Map.Entry<Integer, Set<Integer>> entry : fixedNeighbors.entrySet()) {
            correlationIndex1 = entry.getKey();
            correlation1 = correlationList.get(correlationIndex1);
            if (correlation1.getEquivalence()
                    > 1) {
                continue;
            }
            for (final int correlationIndex2 : entry.getValue()) {
                correlation2 = correlationList.get(correlationIndex2);
                // @TODO for now use fixed neighbor information of atoms without equivalences only
                if (correlation2.getEquivalence()
                        > 1) {
                    continue;
                }
                // insert BOND pair once only and not if equivalences exist
                if (!uniqueSet.contains(indicesMap.get(correlationIndex1)[1]
                                                + " "
                                                + indicesMap.get(correlationIndex2)[1])) {
                    stringBuilder.append("BOND ")
                                 .append(indicesMap.get(correlationIndex1)[1])
                                 .append(" ")
                                 .append(indicesMap.get(correlationIndex2)[1])
                                 .append(buildShiftsComment(correlationList, correlation1, correlation2))
                                 .append("\n");
                    uniqueSet.add(indicesMap.get(correlationIndex1)[1]
                                          + " "
                                          + indicesMap.get(correlationIndex2)[1]);
                    uniqueSet.add(indicesMap.get(correlationIndex2)[1]
                                          + " "
                                          + indicesMap.get(correlationIndex1)[1]);
                }
            }
        }

        return stringBuilder.toString();
    }

    private static String buildBOND(final List<Correlation> correlationList, final Map<Integer, Object[]> indicesMap,
                                    final Map<Integer, Set<Integer>> fixedNeighbors) {
        return buildBONDByFixedNeighbors(correlationList, indicesMap, fixedNeighbors)
                + "\n";
    }

    public static String buildPyLSDInputFileContent(final Correlations correlations, final String mf,
                                                    final Detections detections,
                                                    final ElucidationOptions elucidationOptions) {
        if (mf
                != null) {
            final List<Correlation> correlationList = correlations.getValues();
            final Map<String, Integer> elementCounts = new LinkedHashMap<>(Utils.getMolecularFormulaElementCounts(mf));
            final StringBuilder stringBuilder = new StringBuilder();
            // create header
            stringBuilder.append(buildHeader())
                         .append("\n\n");
            // FORM
            stringBuilder.append(buildFORM(mf, elementCounts))
                         .append("\n\n");
            // PIEC
            stringBuilder.append(buildPIEC())
                         .append("\n\n");
            // ELIM
            if (elucidationOptions.isUseElim()) {
                stringBuilder.append(buildELIM(elucidationOptions.getElimP1(), elucidationOptions.getElimP2()))
                             .append("\n\n");
            }

            final Map<String, List<String>> collection = new LinkedHashMap<>();
            collection.put("MULT", new ArrayList<>());
            collection.put("HSQC", new ArrayList<>());
            collection.put("HMBC", new ArrayList<>());
            collection.put("COSY", new ArrayList<>());
            collection.put("SHIX", new ArrayList<>());
            collection.put("SHIH", new ArrayList<>());
            // index in correlation data -> [atom type, index in PyLSD file]
            final Map<Integer, Object[]> indicesMap = buildIndicesMap(correlationList);

            Correlation correlation;
            for (int i = 0; i
                    < correlationList.size(); i++) {
                correlation = correlationList.get(i);
                collection.get("MULT")
                          .add(buildMULT(correlationList, i, indicesMap, detections.getDetectedHybridizations()));
                collection.get("HSQC")
                          .add(buildHSQC(correlationList, i, indicesMap));
                collection.get("HMBC")
                          .add(buildHMBC(correlationList, i, indicesMap));
                collection.get("COSY")
                          .add(buildCOSY(correlationList, i, indicesMap));
                collection.get("SHIX")
                          .add(buildSHIX(correlation, i, indicesMap));
                collection.get("SHIH")
                          .add(buildSHIH(correlation, i, indicesMap));
            }

            collection.keySet()
                      .forEach(key -> {
                          collection.get(key)
                                    .stream()
                                    .filter(Objects::nonNull)
                                    .forEach(stringBuilder::append);
                          stringBuilder.append("\n");
                      });

            // BOND (interpretation, INADEQUATE, previous assignments) -> input fragments
            stringBuilder.append(buildBOND(correlationList, indicesMap, detections.getFixedNeighbors()))
                         .append("\n");

            // LIST PROP for certain limitations or properties of atoms in lists, e.g. hetero hetero bonds allowance
            stringBuilder.append(
                                 buildLISTsAndPROPs(correlationList, indicesMap, elementCounts, detections.getForbiddenNeighbors(),
                                                    detections.getSetNeighbors(), elucidationOptions.isAllowHeteroHeteroBonds()))
                         .append("\n");
            // DEFF and FEXP as filters (good/bad lists)
            stringBuilder.append(buildDEFFsAndFEXP(correlationList, indicesMap, elucidationOptions,
                                                   detections.getForbiddenNeighbors(), detections.getSetNeighbors()))
                         .append("\n");

            return stringBuilder.toString();
        }

        return "";
    }
}
