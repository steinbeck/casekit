package casekit.nmr.fragments;

import casekit.nmr.model.Assignment;
import casekit.nmr.model.DataSet;
import casekit.nmr.model.Spectrum;
import casekit.nmr.similarity.Similarity;
import casekit.nmr.utils.Utils;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.*;
import java.util.stream.Collectors;

public class FragmentUtilities {
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
                                                                   final boolean checkMultiplicity) {
        final List<DataSet> matches = new ArrayList<>();
        final List<DataSet> nonMatches = new ArrayList<>();
        for (final DataSet dataSet : dataSetList) {
            if (isMatch(dataSet, querySpectrum, mf, shiftTol, maxAverageDeviation, checkMultiplicity)) {
                matches.add(dataSet);
            } else if (isNonMatch(dataSet, querySpectrum, mf, shiftTol, checkMultiplicity)) {
                nonMatches.add(dataSet);
            }
        }
        final Map<String, List<DataSet>> lists = new HashMap<>();
        lists.put("goodlist", matches);
        lists.put("badlist", nonMatches);

        return lists;
    }

    public static boolean isMatch(final DataSet dataSet, final Spectrum querySpectrum, final String mf,
                                  final double shiftTol, final double maxAverageDeviation,
                                  final boolean checkMultiplicity) {
        // check for nuclei
        if (!dataSet.getSpectrum()
                    .getNuclei()[0].equals(querySpectrum.getNuclei()[0])) {
            return false;
        }
        // check for structural problems in fragment regarding the molecular formula
        if (!isStructuralMatch(dataSet, mf)) {
            return false;
        }
        // check average deviation
        final Double averageDeviation = Similarity.calculateAverageDeviation(dataSet.getSpectrum(), querySpectrum, 0, 0,
                                                                             shiftTol, checkMultiplicity, true, true);
        if (averageDeviation
                == null
                || averageDeviation
                > maxAverageDeviation) {
            return false;
        }
        final Double rmsd = Similarity.calculateRMSD(dataSet.getSpectrum(), querySpectrum, 0, 0, shiftTol,
                                                     checkMultiplicity, true, true);
        dataSet.getMeta()
               .put("avgDev", Double.toString(averageDeviation));
        dataSet.getMeta()
               .put("rmsd", Double.toString(rmsd));

        return true;
    }

    public static boolean isNonMatch(final DataSet dataSet, final Spectrum querySpectrum, final String mf,
                                     final double shiftTol, final boolean checkMultiplicity) {
        boolean isSpectralMatch = false;
        if (dataSet.getSpectrum()
                   .getNuclei()[0].equals(querySpectrum.getNuclei()[0])) {
            final Assignment matchAssigment = Similarity.matchSpectra(dataSet.getSpectrum(), querySpectrum, 0, 0,
                                                                      shiftTol, checkMultiplicity, true, true);
            if (matchAssigment
                    != null
                    && matchAssigment.getSetAssignmentsCount(0)
                    > 0) {
                isSpectralMatch = true;
            }
        }

        return !isSpectralMatch
                && !isStructuralMatch(dataSet, mf);
    }

    public static boolean isStructuralMatch(final DataSet dataSet, final String mf) {
        final IAtomContainer fragment = dataSet.getStructure()
                                               .toAtomContainer();
        if (Utils.getAtomTypeFromNucleus(dataSet.getSpectrum()
                                                .getNuclei()[0])
                 .equals("H")) {
            return AtomContainerManipulator.getImplicitHydrogenCount(fragment)
                    <= Utils.getAtomTypeCount(mf, "H");
        } else {
            // check molecular formula with atom types in group
            if (!Utils.compareWithMolecularFormulaLessOrEqual(fragment, mf)) {
                return false;
            }
            // do not allow unsaturated fragments with different size than given molecular formula
            return !Utils.getUnsaturatedAtomIndices(fragment)
                         .isEmpty()
                    || Utils.compareWithMolecularFormulaEqual(fragment, mf);
        }
    }
}
