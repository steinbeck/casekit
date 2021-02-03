package casekit.nmr.lsd;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

public class Constants {

    // valid strings from LSD webpage: C N N5 O S S4 S6 F Cl Br I P P5 Si B X
    public static final Map<String, String> nucleiMap = createNucleiMap();
    public static final Map<String, String> defaultHybridizationMap = createDefaultHybridizationMap();
    public static final Map<String, String> defaultProtonsCountPerValencyMap = createDefaultProtonsCountPerValencyMap();
    public static final Map<String, String> defaultAtomLabelMap = createDefaultAtomLabelMap();
    public static final Map<String, Map<String, Integer>> hybridizationConversionMap = createHybridizationConversionMap();

    private static Map<String, String> createNucleiMap() {
        final Map<String, String> nuclei = new HashMap<>();
        nuclei.put("C", "13C");
        nuclei.put("N", "15N");
        nuclei.put("H", "1H");

        return Collections.unmodifiableMap(nuclei);
    }

    private static Map<String, String> createDefaultHybridizationMap() {
        final Map<String, String> defaultHybridization = new HashMap<>();
        defaultHybridization.put("C", "(1 2 3)");
        defaultHybridization.put("N", "(1 2 3)");
        defaultHybridization.put("O", "(2 3)");
        defaultHybridization.put("S", "(1 2 3)");
        defaultHybridization.put("I", "3");

        return Collections.unmodifiableMap(defaultHybridization);
    }

    private static Map<String, String> createDefaultProtonsCountPerValencyMap() {
        final Map<String, String> defaultProtonsCountPerValency = new HashMap<>();
        defaultProtonsCountPerValency.put("C", "(0 1 2 3)");
        defaultProtonsCountPerValency.put("N", "(0 1 2)");
        defaultProtonsCountPerValency.put("N5", "(0 1 2 3)");
        defaultProtonsCountPerValency.put("N35", "(0 1 2 3)");
        defaultProtonsCountPerValency.put("S", "(0 1)");
        defaultProtonsCountPerValency.put("S4", "(0 1 2 3)");
        defaultProtonsCountPerValency.put("S6", "(0 1 2 3)");
        defaultProtonsCountPerValency.put("S246", "(0 1 2 3)");
        defaultProtonsCountPerValency.put("O", "(0 1)");
        defaultProtonsCountPerValency.put("I", "(0 1)");

        return defaultProtonsCountPerValency;
    }

    private static Map<String, String> createDefaultAtomLabelMap() {
        final Map<String, String> defaultAtomLabel = new HashMap<>();
        defaultAtomLabel.put("C", "C");
        defaultAtomLabel.put("N", "N35");
        defaultAtomLabel.put("O", "O");
        defaultAtomLabel.put("S", "S246");
        defaultAtomLabel.put("I", "I");

        return Collections.unmodifiableMap(defaultAtomLabel);
    }

    private static Map<String, Map<String, Integer>> createHybridizationConversionMap() {
        // @TODO access this information from MongoDB and store it instead of hard coding it
        // possible command in MongoDB: db.hybridizations.aggregate([{$match: {nucleus: "15N"}}, {$group: {_id: null, set: {$addToSet: "$hybridization"}}}])
        // nucleus -> hybridization string -> number
        final Map<String, Map<String, Integer>> hybridizationConversionMap = new HashMap<>();
        hybridizationConversionMap.put("13C", new HashMap<>());
        hybridizationConversionMap.get("13C").put("PLANAR3", 3);
        hybridizationConversionMap.get("13C").put("SP3", 3);
        hybridizationConversionMap.get("13C").put("SP2", 2);
        hybridizationConversionMap.get("13C").put("SP1", 1);
        hybridizationConversionMap.put("15N", new HashMap<>());
        hybridizationConversionMap.get("15N").put("PLANAR3", 3);
        hybridizationConversionMap.get("15N").put("SP3", 3);
        hybridizationConversionMap.get("15N").put("SP2", 2);
        hybridizationConversionMap.get("15N").put("SP1", 1);

        return Collections.unmodifiableMap(hybridizationConversionMap);
    }
}
