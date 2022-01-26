package casekit.nmr.lsd;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

public class Constants {

    // valid strings from LSD webpage: C N N5 O S S4 S6 F Cl Br I P P5 Si B X
    public static final Map<String, String> nucleiMap = createNucleiMap();
    public static final Map<String, int[]> defaultHybridizationMap = createDefaultHybridizationMap();
    public static final Map<String, int[]> defaultProtonsCountPerValencyMap = createDefaultProtonsCountPerValencyMap();
    public static final Map<String, String> defaultAtomLabelMap = createDefaultAtomLabelMap();
    public static final Map<String, Integer> hybridizationConversionMap = createHybridizationConversionMap();

    private static Map<String, String> createNucleiMap() {
        final Map<String, String> nuclei = new HashMap<>();
        nuclei.put("C", "13C");
        nuclei.put("N", "15N");
        nuclei.put("H", "1H");
        nuclei.put("S", "33S");
        nuclei.put("F", "19F");
        nuclei.put("P", "31P");
        nuclei.put("Si", "29Si");

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
        defaultHybridization.put("Cl", new int[]{3});
        defaultHybridization.put("Br", new int[]{3});
        defaultHybridization.put("P", new int[]{1, 2, 3});
        defaultHybridization.put("Si", new int[]{1, 2, 3});

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
        defaultProtonsCountPerValency.put("Cl", new int[]{0});
        defaultProtonsCountPerValency.put("Br", new int[]{0});
        defaultProtonsCountPerValency.put("P", new int[]{0, 1, 2});
        defaultProtonsCountPerValency.put("P5", new int[]{0, 1, 2, 3});
        defaultProtonsCountPerValency.put("P35", new int[]{0, 1, 2, 3});
        defaultProtonsCountPerValency.put("Si", new int[]{0, 1, 2, 3});

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
        defaultAtomLabel.put("Cl", "Cl");
        defaultAtomLabel.put("Br", "Br");
        defaultAtomLabel.put("P", "P35");
        defaultAtomLabel.put("Si", "Si");

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
}
