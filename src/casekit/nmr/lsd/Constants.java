package casekit.nmr.lsd;

import java.util.*;

public class Constants {

    // valid strings from LSD webpage: C N N5 O S S4 S6 F Cl Br I P P5 Si B X
    public static final Map<String, String> nucleiMap = createNucleiMap();
    public static final Map<String, int[]> defaultHybridizationMap = createDefaultHybridizationMap();
    public static final Map<String, int[]> defaultProtonsCountPerValencyMap = createDefaultProtonsCountPerValencyMap();
    public static final Map<String, String> defaultAtomLabelMap = createDefaultAtomLabelMap();
    public static final Map<String, Integer> hybridizationConversionMap = createHybridizationConversionMap();
    public static final Map<String, Map<Integer, Set<Integer>>> hybridizationsByProtonsCountMap = createHybridizationsByProtonsCountMap();

    private static Map<String, String> createNucleiMap() {
        final Map<String, String> nuclei = new HashMap<>();
        nuclei.put("C", "13C");
        nuclei.put("N", "15N");
        nuclei.put("H", "1H");

        return Collections.unmodifiableMap(nuclei);
    }

    private static Map<String, int[]> createDefaultHybridizationMap() {
        final Map<String, int[]> defaultHybridization = new HashMap<>();
        defaultHybridization.put("C", new int[]{1, 2, 3});
        defaultHybridization.put("N", new int[]{1, 2, 3});
        defaultHybridization.put("S", new int[]{1, 2, 3});
        defaultHybridization.put("O", new int[]{2, 3});
        defaultHybridization.put("I", new int[]{3});
        defaultHybridization.put("F", new int[]{3});

        return Collections.unmodifiableMap(defaultHybridization);
    }

    private static Map<String, int[]> createDefaultProtonsCountPerValencyMap() {
        final Map<String, int[]> defaultProtonsCountPerValency = new HashMap<>();
        defaultProtonsCountPerValency.put("C", new int[]{0, 1, 2, 3});
        defaultProtonsCountPerValency.put("N", new int[]{0, 1, 2});
        defaultProtonsCountPerValency.put("N5", new int[]{0, 1, 2, 3});
        defaultProtonsCountPerValency.put("N35", new int[]{0, 1, 2, 3});
        defaultProtonsCountPerValency.put("S", new int[]{0, 1});
        defaultProtonsCountPerValency.put("S4", new int[]{0, 1, 2, 3});
        defaultProtonsCountPerValency.put("S6", new int[]{0, 1, 2, 3});
        defaultProtonsCountPerValency.put("S246", new int[]{0, 1, 2, 3});
        defaultProtonsCountPerValency.put("O", new int[]{0, 1});
        defaultProtonsCountPerValency.put("I", new int[]{0});
        defaultProtonsCountPerValency.put("F", new int[]{0});

        return defaultProtonsCountPerValency;
    }

    private static Map<String, String> createDefaultAtomLabelMap() {
        final Map<String, String> defaultAtomLabel = new HashMap<>();
        defaultAtomLabel.put("C", "C");
        defaultAtomLabel.put("N", "N35");
        defaultAtomLabel.put("O", "O");
        defaultAtomLabel.put("S", "S246");
        defaultAtomLabel.put("I", "I");
        defaultAtomLabel.put("F", "F");

        return Collections.unmodifiableMap(defaultAtomLabel);
    }

    private static Map<String, Integer> createHybridizationConversionMap() {
        // @TODO access this information from MongoDB and store it instead of hard coding it
        // possible command in MongoDB: db.hybridizations.aggregate([{$match: {nucleus: "15N"}}, {$group: {_id: null, set: {$addToSet: "$hybridization"}}}])
        // nucleus -> hybridization string -> number
        final Map<String, Integer> hybridizationConversionMap = new HashMap<>();
        hybridizationConversionMap.put("PLANAR3", 3);
        hybridizationConversionMap.put("SP3", 3);
        hybridizationConversionMap.put("SP2", 2);
        hybridizationConversionMap.put("SP1", 1);

        return Collections.unmodifiableMap(hybridizationConversionMap);
    }

    private static Map<String, Map<Integer, Set<Integer>>> createHybridizationsByProtonsCountMap() {
        final Map<String, Map<Integer, Set<Integer>>> hybridizationsByProtonsCountMap = new HashMap<>();
        hybridizationsByProtonsCountMap.put("C", new HashMap<>());
        hybridizationsByProtonsCountMap.get("C")
                                       .put(4, new HashSet<>());
        hybridizationsByProtonsCountMap.get("C")
                                       .get(4)
                                       .add(3);
        hybridizationsByProtonsCountMap.get("C")
                                       .put(3, new HashSet<>());
        hybridizationsByProtonsCountMap.get("C")
                                       .get(3)
                                       .add(3);
        hybridizationsByProtonsCountMap.get("C")
                                       .put(2, new HashSet<>());
        hybridizationsByProtonsCountMap.get("C")
                                       .get(2)
                                       .add(3);
        hybridizationsByProtonsCountMap.get("C")
                                       .get(2)
                                       .add(2);
        hybridizationsByProtonsCountMap.get("C")
                                       .put(1, new HashSet<>());
        hybridizationsByProtonsCountMap.get("C")
                                       .get(1)
                                       .add(3);
        hybridizationsByProtonsCountMap.get("C")
                                       .get(1)
                                       .add(2);
        hybridizationsByProtonsCountMap.get("C")
                                       .get(1)
                                       .add(1);
        hybridizationsByProtonsCountMap.get("C")
                                       .put(0, new HashSet<>());
        hybridizationsByProtonsCountMap.get("C")
                                       .get(0)
                                       .add(3);
        hybridizationsByProtonsCountMap.get("C")
                                       .get(0)
                                       .add(2);
        hybridizationsByProtonsCountMap.get("C")
                                       .get(0)
                                       .add(1);
        // N (3)
        hybridizationsByProtonsCountMap.put("N", new HashMap<>());
        hybridizationsByProtonsCountMap.get("N")
                                       .put(3, new HashSet<>());
        hybridizationsByProtonsCountMap.get("C")
                                       .get(3)
                                       .add(3);
        hybridizationsByProtonsCountMap.get("N")
                                       .put(2, new HashSet<>());
        hybridizationsByProtonsCountMap.get("C")
                                       .get(2)
                                       .add(3);
        hybridizationsByProtonsCountMap.get("C")
                                       .get(2)
                                       .add(2);

        return Collections.unmodifiableMap(hybridizationsByProtonsCountMap);
    }
}
