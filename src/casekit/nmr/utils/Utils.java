package casekit.nmr.utils;

import casekit.nmr.model.Spectrum;
import casekit.nmr.model.nmrdisplayer.Correlation;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Utils {

    public static String getMultiplicityFromProtonsCount(final Correlation correlation) {
        if (correlation.getAtomType().equals("C") && correlation.getProtonsCount().size() == 1) {
            switch (correlation.getProtonsCount().get(0)) {
                case 0:
                    return "s";
                case 1:
                    return "d";
                case 2:
                    return "t";
                case 3:
                    return "q";
                default:
                    return null;
            }
        }
        return null;
    }

    public static String getAtomTypeFromSpectrum(final Spectrum spectrum, final int dim) {
        if (spectrum.containsDim(dim)) {
            return getAtomTypeFromNucleus(spectrum.getNuclei()[dim]);
        }

        return null;
    }

    public static String getAtomTypeFromNucleus(final String nucleus) {
        final String[] nucleusSplit = nucleus.split("\\d");
        return nucleusSplit[nucleusSplit.length - 1];
    }

    public static Map<String, Integer> getMolecularFormulaElementCounts(final String mf) {
        final LinkedHashMap<String, Integer> counts = new LinkedHashMap<>();
        final IMolecularFormula iMolecularFormula = Utils.getMolecularFormulaFromString(mf);
        final List<String> elements = new ArrayList<>();
        final Matcher matcher = Pattern.compile("([A-Z][a-z]*)").matcher(mf);

        while (matcher.find()) {
            elements.add(matcher.group(1));
        }
        for (final String element : elements) {
            counts.put(element, MolecularFormulaManipulator.getElementCount(iMolecularFormula, element));
        }

        return counts;
    }

    public static IMolecularFormula getMolecularFormulaFromString(final String mf) {
        return MolecularFormulaManipulator.getMolecularFormula(mf, SilentChemObjectBuilder.getInstance());
    }
    

}
