package casekit.nmr.lsd;

import casekit.nmr.lsd.model.ElucidationOptions;
import casekit.nmr.model.Signal;
import casekit.nmr.model.nmrium.Correlation;
import casekit.nmr.model.nmrium.Data;
import casekit.nmr.model.nmrium.Link;
import casekit.nmr.utils.Statistics;
import casekit.nmr.utils.Utils;

import java.text.SimpleDateFormat;
import java.util.*;

public class PyLSDInputFileBuilder {

    private static String buildHeader() {
        final StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append("; PyLSD input file created by webCASE\n");
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

    private static Map<Integer, Object[]> buildIndicesMap(final List<Correlation> correlationList,
                                                          final Map<String, Integer> elementCounts) {
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
                            if (protonsCount
                                    == 3) {
                                protonsCount = 1;
                            }
                            protonsToInsert += correlationList.get(matchIndex)
                                                              .getEquivalence()
                                    * (protonsCount
                                    / correlationList.get(matchIndex)
                                                     .getAttachment()
                                                     .get("H")
                                                     .size());
                        }
                    }
                }
                indicesMap.put(i, new Object[1
                        + protonsToInsert]);
                indicesMap.get(i)[0] = correlation.getAtomType();
                for (int j = 0; j
                        < protonsToInsert; j++) {
                    indicesMap.get(i)[1
                            + j] = protonIndexInPyLSDFile;
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

        return indicesMap;
    }

    private static String buildMULT(final Correlation correlation, final int index,
                                    final Map<Integer, Object[]> indicesMap,
                                    final Map<Integer, List<Integer>> detectedHybridizations) {
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
                             .append(buildShiftString(correlation));
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

    private static String buildShiftString(final Correlation correlation) {

        final Signal signal = Utils.extractSignalFromCorrelation(correlation);
        if (signal
                == null) {
            return "?";
        }

        return correlation.isPseudo()
               ? "?"
               : String.valueOf(Statistics.roundDouble(signal.getShift(0), 3));
    }

    private static String buildShiftsComment(final Correlation correlation1, final Correlation correlation2) {
        return "; "
                + correlation1.getAtomType()
                + ": "
                + buildShiftString(correlation1)
                + " -> "
                + correlation2.getAtomType()
                + ": "
                + buildShiftString(correlation2);
    }

    private static String buildHSQC(final List<Correlation> correlationList, final int index,
                                    final Map<Integer, Object[]> indicesMap) {
        final Correlation correlation = correlationList.get(index);
        if (correlation.getAtomType()
                       .equals("H")) {
            return null;
        }
        final StringBuilder stringBuilder = new StringBuilder();
        final Map<Integer, Integer> protonEquivalenceIndexMap = new HashMap<>();
        int protonsCount;
        for (final Link link : correlation.getLink()) {
            if (link.getExperimentType()
                    .equals("hsqc")
                    || link.getExperimentType()
                           .equals("hmqc")) {
                for (final int matchIndex : link.getMatch()) {
                    protonEquivalenceIndexMap.putIfAbsent(matchIndex, 1); // k = 1 in indicesMap
                    // for each equivalence of heavy atom and attached protons
                    for (int k = 1; k
                            < indicesMap.get(index).length; k++) {
                        protonsCount = correlation.getProtonsCount()
                                                  .get(0);
                        // consider CH3 same as CH (avoid multiple entries in PyLSD input file)
                        if (protonsCount
                                == 3) {
                            protonsCount = 1;
                        }
                        for (int p = 0; p
                                < Math.min(protonsCount, correlationList.get(matchIndex)
                                                                        .getEquivalence()); p++) {
                            stringBuilder.append("HSQC ")
                                         .append(indicesMap.get(index)[k])
                                         .append(" ")
                                         .append(indicesMap.get(matchIndex)[protonEquivalenceIndexMap.get(matchIndex)])
                                         .append(buildShiftsComment(correlation, correlationList.get(matchIndex)))
                                         .append("\n");
                            protonEquivalenceIndexMap.put(matchIndex, protonEquivalenceIndexMap.get(matchIndex)
                                    + 1);
                        }
                    }
                }
            }
        }

        return stringBuilder.toString();
    }

    private static String buildHMBC(final List<Correlation> correlationList, final int index,
                                    final Map<Integer, Object[]> indicesMap, final int hmbcP3, final int hmbcP4) {
        final Correlation correlation = correlationList.get(index);
        if (correlation.getAtomType()
                       .equals("H")) {
            return null;
        }
        final String defaultBondDistance = hmbcP3
                + " "
                + hmbcP4;
        final Set<String> uniqueSet = new LinkedHashSet<>(); // in case of same content exists multiple times
        for (final Link link : correlation.getLink()) {
            if (link.getExperimentType()
                    .equals("hmbc")) {
                for (final int matchIndex : link.getMatch()) {
                    for (int k = 1; k
                            < indicesMap.get(index).length; k++) {
                        for (int l = 1; l
                                < indicesMap.get(matchIndex).length; l++) {
                            // only add an HMBC correlation if there is no direct link via HSQC and the equivalence index is not equal
                            if (!(correlationList.get(matchIndex)
                                                 .getAttachment()
                                                 .containsKey(correlation.getAtomType())
                                    && correlationList.get(matchIndex)
                                                      .getAttachment()
                                                      .get(correlation.getAtomType())
                                                      .contains(index)
                                    && l
                                    == k)) {
                                uniqueSet.add(indicesMap.get(index)[k]
                                                      + " "
                                                      + indicesMap.get(matchIndex)[l]
                                                      + " "
                                                      + defaultBondDistance
                                                      + buildShiftsComment(correlation,
                                                                           correlationList.get(matchIndex)));
                            }
                        }
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
                                    final Map<Integer, Object[]> indicesMap, final int cosyP3, final int cosyP4) {
        final Correlation correlation = correlationList.get(index);
        if (!correlation.getAtomType()
                        .equals("H")) {
            return null;
        }
        final String defaultBondDistance = cosyP3
                + " "
                + cosyP4;
        final Set<String> uniqueSet = new LinkedHashSet<>(); // in case of same content exists multiple times
        for (final Link link : correlation.getLink()) {
            if (link.getExperimentType()
                    .equals("cosy")) {
                for (final int matchIndex : link.getMatch()) {
                    // only add an COSY correlation if the two signals there is not equivalent
                    if (!correlationList.get(matchIndex)
                                        .getId()
                                        .equals(correlation.getId())) {
                        for (int k = 1; k
                                < indicesMap.get(index).length; k++) {
                            // only allow COSY values between possible equivalent protons and only one another non-equivalent proton
                            if (indicesMap.get(matchIndex).length
                                    == 2) {
                                uniqueSet.add(indicesMap.get(index)[k]
                                                      + " "
                                                      + indicesMap.get(matchIndex)[1]
                                                      + " "
                                                      + defaultBondDistance
                                                      + buildShiftsComment(correlation,
                                                                           correlationList.get(matchIndex)));
                            }
                        }
                    }
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
                                             final Map<Integer, Map<String, Map<Integer, Set<Integer>>>> detectedConnectivities,
                                             final Map<Integer, Map<String, Map<Integer, Set<Integer>>>> forbiddenNeighbors,
                                             final Map<Integer, Map<String, Map<Integer, Set<Integer>>>> setNeighbors,
                                             final boolean allowHeteroHeteroBonds) {
        final StringBuilder stringBuilder = new StringBuilder();
        final Map<String, String> listMap = new HashMap<>();

        // LIST and PROP for hetero hetero bonds to disallow
        if (!allowHeteroHeteroBonds) {
            LISTAndPROPUtilities.insertNoHeteroHeteroBonds(stringBuilder, listMap);
        }
        // insert ELEM for each heavy atom type in MF
        LISTAndPROPUtilities.insertELEM(stringBuilder, listMap, elementCounts.keySet());
        //        // insert list combinations of carbon and hybridization states
        //        LISTAndPROPUtilities.insertHeavyAtomCombinationLISTs(stringBuilder, listMap, correlationList, indicesMap);
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
        stringBuilder.append("; externally defined filters\n");
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
        // build and write neighbors files
        final List<String> pathsToNeighborsFilesToUse = new ArrayList<>();
        if (Utilities.writeNeighborsFile(elucidationOptions.getPathsToNeighborsFiles()[0], correlationList, indicesMap,
                                         forbiddenNeighbors)) {
            fexpMap.put("F"
                                + (fexpMap.size()
                    + 1), false);
            pathsToNeighborsFilesToUse.add(elucidationOptions.getPathsToNeighborsFiles()[0]);
        }
        if (Utilities.writeNeighborsFile(elucidationOptions.getPathsToNeighborsFiles()[1], correlationList, indicesMap,
                                         setNeighbors)) {
            fexpMap.put("F"
                                + (fexpMap.size()
                    + 1), true);
            pathsToNeighborsFilesToUse.add(elucidationOptions.getPathsToNeighborsFiles()[1]);
        }
        // build DEFFs
        stringBuilder.append(
                             buildDEFFs(elucidationOptions.getFilterPaths(), pathsToNeighborsFilesToUse.toArray(String[]::new)))
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
                                 .append(buildShiftsComment(correlation1, correlation2))
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

    public static String buildPyLSDInputFileContent(final Data data, final String mf,
                                                    final Map<Integer, List<Integer>> detectedHybridizations,
                                                    final Map<Integer, Map<String, Map<Integer, Set<Integer>>>> detectedConnectivities,
                                                    final Map<Integer, Map<String, Map<Integer, Set<Integer>>>> forbiddenNeighbors,
                                                    final Map<Integer, Map<String, Map<Integer, Set<Integer>>>> setNeighbors,
                                                    final Map<Integer, Set<Integer>> fixedNeighbors,
                                                    final ElucidationOptions elucidationOptions) {
        final Map<String, Map<String, Object>> state = data.getCorrelations()
                                                           .getState();
        final boolean hasErrors = state.keySet()
                                       .stream()
                                       .anyMatch(atomType -> state.get(atomType)
                                                                  .containsKey("error")
                                               && !((Map<String, Object>) state.get(atomType)
                                                                               .get("error")).isEmpty());
        if (mf
                != null
                && !hasErrors) {
            final List<Correlation> correlationList = data.getCorrelations()
                                                          .getValues();
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
            final Map<Integer, Object[]> indicesMap = buildIndicesMap(correlationList, elementCounts);

            Correlation correlation;
            for (int i = 0; i
                    < correlationList.size(); i++) {
                correlation = correlationList.get(i);
                collection.get("MULT")
                          .add(buildMULT(correlation, i, indicesMap, detectedHybridizations));
                collection.get("HSQC")
                          .add(buildHSQC(correlationList, i, indicesMap));
                collection.get("HMBC")
                          .add(buildHMBC(correlationList, i, indicesMap, elucidationOptions.getHmbcP3(),
                                         elucidationOptions.getHmbcP4()));
                collection.get("COSY")
                          .add(buildCOSY(correlationList, i, indicesMap, elucidationOptions.getCosyP3(),
                                         elucidationOptions.getCosyP4()));
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
            stringBuilder.append(buildBOND(correlationList, indicesMap, fixedNeighbors))
                         .append("\n");

            // LIST PROP for certain limitations or properties of atoms in lists, e.g. hetero hetero bonds allowance
            stringBuilder.append(buildLISTsAndPROPs(correlationList, indicesMap, elementCounts, detectedConnectivities,
                                                    forbiddenNeighbors, setNeighbors,
                                                    elucidationOptions.isAllowHeteroHeteroBonds()))
                         .append("\n");
            // DEFF and FEXP as filters (good/bad lists)
            stringBuilder.append(buildDEFFsAndFEXP(correlationList, indicesMap, elucidationOptions, forbiddenNeighbors,
                                                   setNeighbors))
                         .append("\n");

            return stringBuilder.toString();
        }

        return "";
    }
}
