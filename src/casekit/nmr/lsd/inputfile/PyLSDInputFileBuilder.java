package casekit.nmr.lsd.inputfile;

import casekit.nmr.lsd.Constants;
import casekit.nmr.lsd.Utilities;
import casekit.nmr.lsd.model.Detections;
import casekit.nmr.lsd.model.ElucidationOptions;
import casekit.nmr.lsd.model.Grouping;
import casekit.nmr.lsd.model.MolecularConnectivity;
import casekit.nmr.model.nmrium.Correlation;
import casekit.nmr.model.nmrium.Correlations;
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

    private static String buildPossibilitiesString(final Set<Integer> possibilities) {
        final StringBuilder possibilitiesStringBuilder = new StringBuilder();

        if (possibilities.size()
                > 1) {
            possibilitiesStringBuilder.append("(");
        }
        int counter = 0;
        for (final int possibility : possibilities) {
            possibilitiesStringBuilder.append(possibility);
            if (counter
                    < possibilities.size()
                    - 1) {
                possibilitiesStringBuilder.append(" ");
            }
            counter++;
        }
        if (possibilities.size()
                > 1) {
            possibilitiesStringBuilder.append(")");
        }

        return possibilitiesStringBuilder.toString();
    }

    private static Map<String, StringBuilder> buildStringBuilderMap(
            final Map<Integer, List<MolecularConnectivity>> molecularConnectivityMap) {
        StringBuilder stringBuilder;
        final Map<String, StringBuilder> stringBuilderMap = new HashMap<>();
        stringBuilderMap.put("MULT", new StringBuilder());
        stringBuilderMap.put("HSQC", new StringBuilder());
        stringBuilderMap.put("HMBC", new StringBuilder());
        stringBuilderMap.put("COSY", new StringBuilder());
        stringBuilderMap.put("BOND", new StringBuilder());
        stringBuilderMap.put("SHIX", new StringBuilder());
        stringBuilderMap.put("SHIH", new StringBuilder());
        StringBuilder hybridizationStringBuilder, attachedProtonsCountStringBuilder;
        int counter, firstOfEquivalenceIndexPyLSD;
        Set<Integer> groupMembers;
        MolecularConnectivity molecularConnectivityGroupMember, molecularConnectivityHeavyAtom;
        for (final int correlationIndex : molecularConnectivityMap.keySet()) {
            firstOfEquivalenceIndexPyLSD = -1;
            for (final MolecularConnectivity molecularConnectivity : molecularConnectivityMap.get(correlationIndex)) {
                if (firstOfEquivalenceIndexPyLSD
                        == -1) {
                    firstOfEquivalenceIndexPyLSD = molecularConnectivity.getIndex();
                }
                if (!molecularConnectivity.getAtomType()
                                          .equals("H")) {

                    hybridizationStringBuilder = new StringBuilder();
                    if (molecularConnectivity.getHybridizations()
                                             .size()
                            > 1) {
                        hybridizationStringBuilder.append("(");
                    }
                    counter = 0;
                    for (final int hybrid : molecularConnectivity.getHybridizations()) {
                        hybridizationStringBuilder.append(hybrid);
                        if (counter
                                < molecularConnectivity.getHybridizations()
                                                       .size()
                                - 1) {
                            hybridizationStringBuilder.append(" ");
                        }
                        counter++;
                    }
                    if (molecularConnectivity.getHybridizations()
                                             .size()
                            > 1) {
                        hybridizationStringBuilder.append(")");
                    }
                    attachedProtonsCountStringBuilder = new StringBuilder();
                    if (molecularConnectivity.getProtonCounts()
                                             .size()
                            > 1) {
                        attachedProtonsCountStringBuilder.append("(");
                    }
                    counter = 0;
                    for (final int protonCount : molecularConnectivity.getProtonCounts()) {
                        attachedProtonsCountStringBuilder.append(protonCount);
                        if (counter
                                < molecularConnectivity.getProtonCounts()
                                                       .size()
                                - 1) {
                            attachedProtonsCountStringBuilder.append(" ");
                        }
                        counter++;
                    }
                    if (molecularConnectivity.getProtonCounts()
                                             .size()
                            > 1) {
                        attachedProtonsCountStringBuilder.append(")");
                    }
                    stringBuilder = stringBuilderMap.get("MULT");
                    stringBuilder.append("MULT ")
                                 .append(molecularConnectivity.getIndex())
                                 .append(" ")
                                 .append(Constants.defaultAtomLabelMap.get(molecularConnectivity.getAtomType()))
                                 .append(" ")
                                 .append(hybridizationStringBuilder)
                                 .append(" ")
                                 .append(attachedProtonsCountStringBuilder);
                    stringBuilder.append("; ")
                                 .append(buildShiftString(molecularConnectivityMap, molecularConnectivity));
                    if (molecularConnectivityMap.get(correlationIndex)
                                                .size()
                            > 1
                            && molecularConnectivity.getIndex()
                            != firstOfEquivalenceIndexPyLSD) {
                        stringBuilder.append("; equivalent to ")
                                     .append(firstOfEquivalenceIndexPyLSD);
                    }
                    stringBuilder.append("\n");
                    if (molecularConnectivity.getHsqc()
                            != null) {
                        stringBuilder = stringBuilderMap.get("HSQC");
                        for (final int protonIndexPyLSD : molecularConnectivity.getHsqc()) {
                            stringBuilder.append("HSQC ")
                                         .append(molecularConnectivity.getIndex())
                                         .append(" ")
                                         .append(protonIndexPyLSD)
                                         .append(buildShiftsComment(molecularConnectivityMap, molecularConnectivity,
                                                                    Utilities.findMolecularConnectivityByIndex(
                                                                            molecularConnectivityMap, "H", false,
                                                                            protonIndexPyLSD)))
                                         .append("\n");
                        }
                    }
                    if (molecularConnectivity.getHmbc()
                            != null) {
                        stringBuilder = stringBuilderMap.get("HMBC");
                        for (final int protonIndexInPyLSD : molecularConnectivity.getHmbc()
                                                                                 .keySet()) {
                            // filter out group members which are directly bonded to that proton
                            groupMembers = new HashSet<>(molecularConnectivity.getGroupMembers());
                            for (final int groupMemberIndex : new HashSet<>(groupMembers)) {
                                molecularConnectivityGroupMember = Utilities.findMolecularConnectivityByIndex(
                                        molecularConnectivityMap, molecularConnectivity.getAtomType(), false,
                                        groupMemberIndex);
                                if (molecularConnectivityGroupMember.getHsqc()
                                        != null
                                        && molecularConnectivityGroupMember.getHsqc()
                                                                           .contains(protonIndexInPyLSD)) {
                                    groupMembers.remove(groupMemberIndex);
                                }
                            }
                            if (!groupMembers.isEmpty()) {
                                stringBuilder.append("HMBC ")
                                             .append(buildPossibilitiesString(groupMembers))
                                             .append(" ")
                                             .append(protonIndexInPyLSD)
                                             .append(" ")
                                             .append(molecularConnectivity.getHmbc()
                                                                          .get(protonIndexInPyLSD)[0])
                                             .append(" ")
                                             .append(molecularConnectivity.getHmbc()
                                                                          .get(protonIndexInPyLSD)[1])
                                             .append(buildShiftsComment(molecularConnectivityMap, molecularConnectivity,
                                                                        Utilities.findMolecularConnectivityByIndex(
                                                                                molecularConnectivityMap, "H", false,
                                                                                protonIndexInPyLSD)))
                                             .append("\n");
                            }
                        }
                    }
                    if (molecularConnectivity.getFixedNeighbors()
                            != null) {
                        stringBuilder = stringBuilderMap.get("BOND");
                        for (final int bondedIndexInPyLSD : molecularConnectivity.getFixedNeighbors()) {
                            stringBuilder.append("BOND ")
                                         .append(molecularConnectivity.getIndex())
                                         .append(" ")
                                         .append(bondedIndexInPyLSD)
                                         .append(buildShiftsComment(molecularConnectivityMap, molecularConnectivity,
                                                                    Utilities.findMolecularConnectivityByIndex(
                                                                            molecularConnectivityMap, "H", true,
                                                                            bondedIndexInPyLSD)))
                                         .append("\n");
                        }
                    }
                } else if (molecularConnectivity.getAtomType()
                                                .equals("H")) {
                    if (molecularConnectivity.getCosy()
                            != null) {
                        stringBuilder = stringBuilderMap.get("COSY");
                        for (final int protonIndexInPyLSD : molecularConnectivity.getCosy()
                                                                                 .keySet()) {
                            // filter out group members which
                            groupMembers = new HashSet<>(molecularConnectivity.getGroupMembers());
                            // 1) use only one attached proton of a CH2 group (optional)
                            final Set<Integer> alreadyFoundHeavyAtomIndex = new HashSet<>();
                            for (final int groupMemberIndex : new HashSet<>(groupMembers)) {
                                molecularConnectivityHeavyAtom = Utilities.getHeavyAtomMolecularConnectivity(
                                        molecularConnectivityMap, groupMemberIndex);
                                if (alreadyFoundHeavyAtomIndex.contains(molecularConnectivityHeavyAtom.getIndex())) {
                                    groupMembers.remove(groupMemberIndex);
                                } else {
                                    alreadyFoundHeavyAtomIndex.add(molecularConnectivityHeavyAtom.getIndex());
                                }
                            }
                            // 2) would direct to itself when using COSY correlation
                            molecularConnectivityHeavyAtom = Utilities.getHeavyAtomMolecularConnectivity(
                                    molecularConnectivityMap, protonIndexInPyLSD);
                            if (molecularConnectivityHeavyAtom
                                    != null) {
                                for (final int groupMemberIndex : new HashSet<>(groupMembers)) {
                                    if (molecularConnectivityHeavyAtom.getHsqc()
                                                                      .contains(groupMemberIndex)) {
                                        groupMembers.remove(groupMemberIndex);
                                    }
                                }
                                if (!groupMembers.isEmpty()) {
                                    stringBuilder.append("COSY ")
                                                 .append(buildPossibilitiesString(groupMembers))
                                                 .append(" ")
                                                 .append(protonIndexInPyLSD)
                                                 .append(" ")
                                                 .append(molecularConnectivity.getCosy()
                                                                              .get(protonIndexInPyLSD)[0])
                                                 .append(" ")
                                                 .append(molecularConnectivity.getCosy()
                                                                              .get(protonIndexInPyLSD)[1])
                                                 .append(buildShiftsComment(molecularConnectivityMap,
                                                                            molecularConnectivity,
                                                                            Utilities.findMolecularConnectivityByIndex(
                                                                                    molecularConnectivityMap, "H",
                                                                                    false, protonIndexInPyLSD)))
                                                 .append("\n");
                                }
                            }
                        }
                    }
                }
                if (molecularConnectivity.getSignal()
                        != null) {
                    stringBuilder = stringBuilderMap.get(molecularConnectivity.getAtomType()
                                                                              .equals("H")
                                                         ? "SHIH"
                                                         : "SHIX");
                    stringBuilder.append(molecularConnectivity.getAtomType()
                                                              .equals("H")
                                         ? "SHIH"
                                         : "SHIX")
                                 .append(" ")
                                 .append(molecularConnectivity.getIndex())
                                 .append(" ")
                                 .append(Statistics.roundDouble(molecularConnectivity.getSignal()
                                                                                     .getShift(0), 5))
                                 .append("\n");
                }
            }
        }

        return stringBuilderMap;
    }

    private static String buildShiftString(final Map<Integer, List<MolecularConnectivity>> molecularConnectivityMap,
                                           final MolecularConnectivity molecularConnectivity) {
        if (molecularConnectivity
                == null
                || molecularConnectivity.getSignal()
                == null) {
            return "?";
        }

        final String heavyAtomShiftString = "";
        //        if (molecularConnectivity.getAtomType()
        //                                 .equals("H")) {
        //            MolecularConnectivity heavyAtomMolecularConnectivity = null;
        //            boolean found = false;
        //            for (final int correlationIndex : molecularConnectivityMap.keySet()) {
        //                for (final MolecularConnectivity molecularConnectivityTemp : molecularConnectivityMap.get(
        //                        correlationIndex)) {
        //                    if (molecularConnectivityTemp.getHsqc()
        //                            != null
        //                            && molecularConnectivityTemp.getHsqc()
        //                                                        .contains(molecularConnectivity.getIndex())) {
        //                        heavyAtomMolecularConnectivity = molecularConnectivityTemp;
        //                        found = true;
        //                        break;
        //                    }
        //                }
        //                if (found) {
        //                    break;
        //                }
        //            }
        //            if (heavyAtomMolecularConnectivity
        //                    != null) {
        //                heavyAtomShiftString = " ("
        //                        + buildShiftString(molecularConnectivityMap, heavyAtomMolecularConnectivity)
        //                        + ")";
        //            }
        //        }

        return Statistics.roundDouble(molecularConnectivity.getSignal()
                                                           .getShift(0), 3)
                + heavyAtomShiftString;
    }

    private static String buildShiftsComment(final Map<Integer, List<MolecularConnectivity>> molecularConnectivityMap,
                                             final MolecularConnectivity molecularConnectivity1,
                                             final MolecularConnectivity molecularConnectivity2) {
        return "; "
                + molecularConnectivity1.getAtomType()
                + ": "
                + buildShiftString(molecularConnectivityMap, molecularConnectivity1)
                + " -> "
                + molecularConnectivity2.getAtomType()
                + ": "
                + buildShiftString(molecularConnectivityMap, molecularConnectivity2);
    }

    private static String buildLISTsAndPROPs(final Map<Integer, List<MolecularConnectivity>> molecularConnectivityMap,
                                             final Map<String, Integer> elementCounts,
                                             final boolean allowHeteroHeteroBonds) {
        final StringBuilder stringBuilder = new StringBuilder();
        // list key -> [list name, size]
        final Map<String, Object[]> listMap = new HashMap<>();

        // LIST and PROP for hetero hetero bonds to disallow
        if (!allowHeteroHeteroBonds) {
            LISTAndPROPUtilities.insertNoHeteroHeteroBonds(stringBuilder, listMap);
        }
        // insert LIST for each heavy atom type in MF
        LISTAndPROPUtilities.insertGeneralLISTs(stringBuilder, listMap, molecularConnectivityMap,
                                                elementCounts.keySet());
        // insert list combinations of carbon and hybridization states
        LISTAndPROPUtilities.insertHeavyAtomCombinationLISTs(stringBuilder, listMap, molecularConnectivityMap);
        // insert forbidden connection lists and properties
        LISTAndPROPUtilities.insertConnectionLISTsAndPROPs(stringBuilder, listMap, molecularConnectivityMap, "forbid");
        // insert set connection lists and properties
        LISTAndPROPUtilities.insertConnectionLISTsAndPROPs(stringBuilder, listMap, molecularConnectivityMap, "allow");

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

    private static String buildDEFFsAndFEXP(final ElucidationOptions elucidationOptions) {
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

    public static String buildPyLSDInputFileContent(final Correlations correlations, final String mf,
                                                    final Detections detections, final Grouping grouping,
                                                    final ElucidationOptions elucidationOptions) {
        if (mf
                == null) {
            return "";
        }
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

        final Map<String, Integer[]> defaultBondDistances = new HashMap<>();
        defaultBondDistances.put("hmbc", new Integer[]{2, 3});
        defaultBondDistances.put("cosy", new Integer[]{3, 4});
        final Map<Integer, List<MolecularConnectivity>> molecularConnectivityMap = Utilities.buildMolecularConnectivityMap(
                correlationList, detections, grouping, defaultBondDistances);
        final Map<String, StringBuilder> stringBuilderMap = buildStringBuilderMap(molecularConnectivityMap);
        stringBuilder.append(stringBuilderMap.get("MULT")
                                             .toString())
                     .append("\n");
        stringBuilder.append(stringBuilderMap.get("HSQC")
                                             .toString())
                     .append("\n");

        stringBuilder.append(stringBuilderMap.get("BOND")
                                             .toString())
                     .append("\n");
        stringBuilder.append(stringBuilderMap.get("HMBC")
                                             .toString())
                     .append("\n");
        stringBuilder.append(stringBuilderMap.get("COSY")
                                             .toString())
                     .append("\n");
        stringBuilder.append(stringBuilderMap.get("SHIX")
                                             .toString())
                     .append("\n");
        stringBuilder.append(stringBuilderMap.get("SHIH")
                                             .toString())
                     .append("\n");

        // LIST PROP for certain limitations or properties of atoms in lists, e.g. hetero hetero bonds allowance
        stringBuilder.append(buildLISTsAndPROPs(molecularConnectivityMap, elementCounts,
                                                elucidationOptions.isAllowHeteroHeteroBonds()))
                     .append("\n");
        // DEFF and FEXP as filters (good/bad lists)
        stringBuilder.append(buildDEFFsAndFEXP(elucidationOptions))
                     .append("\n");

        return stringBuilder.toString();
    }
}
