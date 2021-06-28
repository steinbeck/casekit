package casekit.nmr.fragments;

import casekit.nmr.model.DataSet;

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
}
