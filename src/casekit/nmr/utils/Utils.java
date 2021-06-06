package casekit.nmr.utils;

import casekit.nmr.model.Spectrum;
import casekit.nmr.model.nmrdisplayer.Correlation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Utils {

    /**
     * Specified for carbons only -> not generic!!!
     *
     * @param protonsCount
     *
     * @return
     */
    public static String getMultiplicityFromProtonsCount(final int protonsCount) {
        switch (protonsCount) {
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

    public static String getMultiplicityFromProtonsCount(final Correlation correlation) {
        if (correlation.getAtomType()
                       .equals("C")
                && correlation.getProtonsCount()
                              .size()
                == 1) {
            return getMultiplicityFromProtonsCount(correlation.getProtonsCount()
                                                              .get(0));
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
        return nucleusSplit[nucleusSplit.length
                - 1];
    }

    public static IMolecularFormula getMolecularFormulaFromString(final String mf) {
        return MolecularFormulaManipulator.getMolecularFormula(mf, SilentChemObjectBuilder.getInstance());
    }

    public static String getSmilesFromAtomContainer(final IAtomContainer ac) throws CDKException {
        // SmiFlavor.Unique instead of SmiFlavor.Absolute because current errors with InChI generator
        final SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Unique);

        return smilesGenerator.create(ac);
    }

    /**
     * @param data
     *
     * @return
     */
    public static Double getRMSD(final Double[] data) {
        if (data
                == null
                || data.length
                == 0) {
            return null;
        }
        if (data.length
                == 1) {
            return data[0];
        }
        int nullCounter = 0;
        double qSum = 0;
        for (final Double d : data) {
            if (d
                    != null) {
                qSum += d
                        * d;
            } else {
                nullCounter++;
            }
        }

        return ((data.length
                - nullCounter)
                != 0)
               ? Math.sqrt(qSum
                                   / (data.length
                - nullCounter))
               : null;
    }

    public static Map<String, Integer> getMolecularFormulaElementCounts(final String mf) {
        final LinkedHashMap<String, Integer> counts = new LinkedHashMap<>();
        final List<String> elements = new ArrayList<>();
        Matcher matcher = Pattern.compile("([A-Z][a-z]{0,1})")
                                 .matcher(mf);
        while (matcher.find()) {
            elements.add(matcher.group(1));
        }
        int count;
        for (final String element : elements) {
            matcher = Pattern.compile("("
                                              + element
                                              + "\\d+)")
                             .matcher(mf);
            count = 1;
            if (matcher.find()) {
                count = Integer.parseInt(matcher.group(1)
                                                .split(element)[1]);
            }
            counts.put(element, count);
        }

        return counts;
    }
}
