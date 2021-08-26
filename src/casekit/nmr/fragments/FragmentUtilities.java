package casekit.nmr.fragments;

import casekit.nmr.model.Assignment;
import casekit.nmr.model.DataSet;
import casekit.nmr.model.Spectrum;
import casekit.nmr.similarity.Similarity;
import casekit.nmr.utils.Utils;
import com.google.gson.Gson;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.*;
import java.util.stream.Collectors;

public class FragmentUtilities {

    private final static Gson gson = new Gson();

    public static LinkedHashMap<String, List<DataSet>> sortByFrequencies(
            final Map<String, List<DataSet>> functionalGroupDataSetsMap) {
        final LinkedHashMap<String, List<DataSet>> sortedCollection = new LinkedHashMap<>();
        final List<Map.Entry<String, List<DataSet>>> sortedFrequencies = getSortedFrequencies(
                functionalGroupDataSetsMap);
        for (final Map.Entry<String, List<DataSet>> frequency : sortedFrequencies) {
            sortedCollection.put(frequency.getKey(), frequency.getValue());
        }

        return sortedCollection;
    }

    private static List<Map.Entry<String, List<DataSet>>> getSortedFrequencies(
            final Map<String, List<DataSet>> functionalGroupDataSets) {
        return functionalGroupDataSets.entrySet()
                                      .stream()
                                      .sorted(Map.Entry.comparingByValue((list1, list2) -> -1
                                              * Integer.compare(list1.size(), list2.size())))
                                      .collect(Collectors.toList());
    }

    public static Map<String, List<DataSet>> collectBySmiles(final List<DataSet> functionalGroupDataSets) {
        final Map<String, List<DataSet>> collection = new HashMap<>();
        String smiles;
        for (final DataSet functionalGroupDataSet : functionalGroupDataSets) {
            smiles = functionalGroupDataSet.getMeta()
                                           .get("smiles");
            if (smiles
                    != null) {
                collection.putIfAbsent(smiles, new ArrayList<>());
                collection.get(smiles)
                          .add(functionalGroupDataSet);
            }
        }

        return collection;
    }

    public static void filterByCommonSubstructures(final Map<String, List<DataSet>> smilesCollection1,
                                                   final Map<String, List<DataSet>> smilesCollection2) {
        // filter first by second collection
        Set<String> keysToRemove = new HashSet<>();
        for (final String keyCollection1 : smilesCollection1.keySet()) {
            if (!smilesCollection2.containsKey(keyCollection1)) {
                keysToRemove.add(keyCollection1);
            }
        }
        for (final String keyCollection1 : keysToRemove) {
            smilesCollection1.remove(keyCollection1);
        }
        // filter second by first collection
        keysToRemove = new HashSet<>();
        for (final String keyCollection2 : smilesCollection2.keySet()) {
            if (!smilesCollection1.containsKey(keyCollection2)) {
                keysToRemove.add(keyCollection2);
            }
        }
        for (final String keyCollection2 : keysToRemove) {
            smilesCollection2.remove(keyCollection2);
        }
    }

    public static Map<String, List<DataSet>> getGoodlistAndBadlist(final List<DataSet> dataSetList,
                                                                   final Spectrum querySpectrum, final String mf,
                                                                   final double shiftTol,
                                                                   final double maxAverageDeviation,
                                                                   final boolean checkMultiplicity,
                                                                   final List<List<String>> queryHybridizationList) {
        final List<DataSet> matches = new ArrayList<>();
        final List<DataSet> nonMatches = new ArrayList<>();
        Assignment matchAssignment;
        for (final DataSet dataSet : dataSetList) {
            matchAssignment = Similarity.matchSpectra(dataSet.getSpectrum()
                                                             .toSpectrum(), querySpectrum, 0, 0, shiftTol,
                                                      checkMultiplicity, true, true);
            if (isMatch(dataSet, querySpectrum, mf, matchAssignment, maxAverageDeviation, queryHybridizationList)) {
                matches.add(dataSet);
            } else if (isNonMatch(dataSet, querySpectrum, mf, matchAssignment)) {
                nonMatches.add(dataSet);
            }
        }
        final Map<String, List<DataSet>> lists = new HashMap<>();
        lists.put("goodlist", matches);
        lists.put("badlist", nonMatches);

        return lists;
    }

    private static boolean isMatch(final DataSet dataSet, final Spectrum querySpectrum, final String mf,
                                   final Assignment matchAssignment, final double maxAverageDeviation,
                                   final List<List<String>> queryHybridizationList) {
        // check for nuclei
        if (!dataSet.getSpectrum()
                    .getNuclei()[0].equals(querySpectrum.getNuclei()[0])) {
            return false;
        }
        // check for structural problems in fragment regarding the molecular formula
        if (!isStructuralMatch(dataSet, mf)) {
            return false;
        }
        final Spectrum spectrum = dataSet.getSpectrum()
                                         .toSpectrum();
        // check average deviation
        final Double averageDeviation = Similarity.calculateAverageDeviation(spectrum, querySpectrum, 0, 0,
                                                                             matchAssignment);
        if (averageDeviation
                == null
                || averageDeviation
                > maxAverageDeviation) {
            return false;
        }
        // check hybridazations after knowing that all signals in dataset have an assignment
        if (!checkHybridizations(dataSet, matchAssignment, queryHybridizationList)) {
            return false;
        }
        dataSet.addMetaInfo("matchAssignment", gson.toJson(matchAssignment, Assignment.class));
        final Double rmsd = Similarity.calculateRMSD(spectrum, querySpectrum, 0, 0, matchAssignment);
        dataSet.addMetaInfo("averageDeviation", Double.toString(averageDeviation));
        dataSet.addMetaInfo("rmsd", Double.toString(rmsd));


        return true;
    }

    private static boolean isNonMatch(final DataSet dataSet, final Spectrum querySpectrum, final String mf,
                                      final Assignment matchAssigment) {
        if (!isStructuralMatch(dataSet, mf)) {
            return false;
        }
        boolean isSpectralMatch = false;
        if (dataSet.getSpectrum()
                   .getNuclei()[0].equals(querySpectrum.getNuclei()[0])) {
            if (matchAssigment
                    != null
                    && matchAssigment.getSetAssignmentsCount(0)
                    > 0) {
                isSpectralMatch = true;
            }
        }

        return !isSpectralMatch;
    }

    private static boolean isStructuralMatch(final DataSet dataSet, final String mf) {
        final IAtomContainer fragment = dataSet.getStructure()
                                               .toAtomContainer();
        // check molecular formula with atom types in group
        // do not allow unsaturated fragments
        return Utils.compareWithMolecularFormulaLessOrEqual(fragment, mf)
                && !Utils.getUnsaturatedAtomIndices(fragment)
                         .isEmpty();
    }

    private static boolean checkHybridizations(final DataSet dataSet, final Assignment matchAssignment,
                                               final List<List<String>> queryHybridizationList) {
        if (queryHybridizationList.isEmpty()) {
            return true;
        }

        final IAtomContainer fragment = dataSet.getStructure()
                                               .toAtomContainer();
        final String atomType = Utils.getAtomTypeFromNucleus(dataSet.getSpectrum()
                                                                    .getNuclei()[0]);
        int signalIndexInDataSetSpectrum, signalIndexInQuerySpectrum;
        for (int i = 0; i
                < fragment.getAtomCount(); i++) {
            if (fragment.getAtom(i)
                        .getSymbol()
                        .equals(atomType)) {
                signalIndexInDataSetSpectrum = dataSet.getAssignment()
                                                      .getIndices(0, i)
                                                      .get(0);
                signalIndexInQuerySpectrum = matchAssignment.getAssignment(0, signalIndexInDataSetSpectrum, 0);
                if (!queryHybridizationList.get(signalIndexInQuerySpectrum)
                                           .contains(fragment.getAtom(i)
                                                             .getHybridization()
                                                             .name())) {
                    return false;
                }
            }
        }

        return true;
    }
}
