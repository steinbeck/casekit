package casekit.nmr.fragments;

import casekit.nmr.model.DataSet;
import casekit.nmr.model.Spectrum;
import casekit.nmr.similarity.Similarity;
import casekit.nmr.utils.Utils;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

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

    public static List<DataSet> findMatches(final List<DataSet> dataSetList, final Spectrum querySpectrum,
                                            final String mf, final double shiftTol, final double maxAverageDeviation,
                                            final boolean checkMultiplicity) {
        final List<DataSet> matches = new ArrayList<>();
        for (final DataSet dataSet : dataSetList) {
            if (isValidMatch(dataSet, querySpectrum, mf, shiftTol, maxAverageDeviation, checkMultiplicity)) {
                matches.add(dataSet);
            }
        }

        return matches;
    }

    public static boolean isValidMatch(final DataSet dataSet, final Spectrum querySpectrum, final String mf,
                                       final double shiftTol, final double maxAverageDeviation,
                                       final boolean checkMultiplicity) {
        final IMolecularFormula iMolecularFormula = MolecularFormulaManipulator.getMolecularFormula(mf,
                                                                                                    SilentChemObjectBuilder.getInstance());
        final IAtomContainer group = dataSet.getStructure()
                                            .toAtomContainer();

        if (!dataSet.getSpectrum()
                    .getNuclei()[0].equals(querySpectrum.getNuclei()[0])) {
            return false;
        }
        final String atomTypeInSpectrum = Utils.getAtomTypeFromNucleus(dataSet.getSpectrum()
                                                                              .getNuclei()[0]);
        if (atomTypeInSpectrum.equals("H")) {
            if (AtomContainerManipulator.getImplicitHydrogenCount(dataSet.getStructure()
                                                                         .toAtomContainer())
                    > MolecularFormulaManipulator.getElementCount(iMolecularFormula, atomTypeInSpectrum)) {
                return false;
            }
        } else {
            // check molecular formula with atom types in group
            if (!Utils.compareWithMolecularFormulaLessOrEqual(group, mf)) {
                return false;
            }
            // do not allow unsaturated fragments with different size than given molecular formula
            if (Utils.getUnsaturatedAtomIndices(group)
                     .isEmpty()
                    && !Utils.compareWithMolecularFormulaEqual(group, mf)) {
                return false;
            }
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
}
