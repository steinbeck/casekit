package casekit.nmr.elucidation.lsd;

import casekit.io.FileSystem;
import casekit.nmr.elucidation.Constants;
import casekit.nmr.model.DataSet;
import casekit.nmr.model.Signal;
import casekit.nmr.model.nmrium.Correlation;
import casekit.nmr.utils.Statistics;
import casekit.nmr.utils.Utils;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.*;
import java.util.stream.Collectors;

public class Utilities {
    public static String buildSSTR(final int sstrIndex, final String atomType, final List<Integer> hybridization,
                                   final List<Integer> protonsCount) {
        if (hybridization.isEmpty()) {
            hybridization.addAll(Arrays.stream(Constants.defaultHybridizationMap.get(atomType))
                                       .boxed()
                                       .collect(Collectors.toList()));
        }
        if (protonsCount.isEmpty()) {
            protonsCount.addAll(Arrays.stream(Constants.defaultProtonsCountPerValencyMap.get(atomType))
                                      .boxed()
                                      .collect(Collectors.toList()));
        }
        final StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append("SSTR S")
                     .append(sstrIndex)
                     .append(" ")
                     .append(atomType)
                     .append(" ");
        if (hybridization.size()
                == 1) {
            stringBuilder.append(hybridization.get(0))
                         .append(" ");
            if (protonsCount.size()
                    == 1) {
                stringBuilder.append(protonsCount.get(0));
            } else {
                stringBuilder.append(buildMultipleValuesString(protonsCount));
            }
        } else {
            stringBuilder.append(buildMultipleValuesString(hybridization));
            stringBuilder.append(" ");
            if (protonsCount.size()
                    == 1) {
                stringBuilder.append(protonsCount.get(0));
            } else {
                stringBuilder.append(buildMultipleValuesString(protonsCount));
            }
        }

        return stringBuilder.toString();
    }

    private static String buildMultipleValuesString(final List<Integer> values) {
        final StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append("(");
        for (int l = 0; l
                < values.size(); l++) {
            stringBuilder.append(values.get(l));
            if (l
                    < values.size()
                    - 1) {
                stringBuilder.append(" ");
            }
        }
        stringBuilder.append(")");

        return stringBuilder.toString();
    }

    public static boolean writeNeighborsFile(final String pathToNeighborsFile, final List<Correlation> correlationList,
                                             final Map<Integer, Object[]> indicesMap,
                                             final Map<Integer, Map<String, Map<Integer, Set<Integer>>>> neighbors) {
        final StringBuilder stringBuilder = new StringBuilder();
        Correlation correlation;
        Signal signal;
        String atomType;
        int indexInPyLSD;
        int sstrIndex = 1, sstrIndexCorrelation;
        Map<String, Map<Integer, Set<Integer>>> neighborsTemp;
        for (int i = 0; i
                < correlationList.size(); i++) {
            if (neighbors.containsKey(i)) {
                correlation = correlationList.get(i);
                signal = Utils.extractFirstSignalFromCorrelation(correlation);
                atomType = correlation.getAtomType();
                neighborsTemp = neighbors.get(i);

                // put in the extracted information per correlation and equivalent
                for (int k = 1; k
                        < indicesMap.get(i).length; k++) {
                    indexInPyLSD = (int) indicesMap.get(i)[k];
                    for (final String neighborAtomType : neighborsTemp.keySet()) {
                        for (final Map.Entry<Integer, Set<Integer>> entryPerHybridization : neighborsTemp.get(
                                                                                                                 neighborAtomType)
                                                                                                         .entrySet()) {
                            sstrIndexCorrelation = sstrIndex;
                            stringBuilder.append(
                                    buildSSTR(sstrIndexCorrelation, atomType, correlation.getHybridization(),
                                              correlation.getProtonsCount()));
                            stringBuilder.append("; ")
                                         .append(atomType)
                                         .append(" at ")
                                         .append(signal
                                                         != null
                                                 ? Statistics.roundDouble(signal.getShift(0), 2)
                                                 : "?")
                                         .append(" (")
                                         .append(indexInPyLSD)
                                         .append(")")
                                         .append("\n");
                            stringBuilder.append("ASGN S")
                                         .append(sstrIndexCorrelation)
                                         .append(" ")
                                         .append(indexInPyLSD)
                                         .append("\n");
                            sstrIndex++;

                            final List<Integer> tempList = new ArrayList<>();
                            if (entryPerHybridization.getKey()
                                    != -1) {
                                tempList.add(entryPerHybridization.getKey());
                            }
                            stringBuilder.append(buildSSTR(sstrIndex, neighborAtomType, tempList,
                                                           new ArrayList<>(entryPerHybridization.getValue())))
                                         .append("\n");
                            stringBuilder.append("LINK S")
                                         .append(sstrIndexCorrelation)
                                         .append(" S")
                                         .append(sstrIndex)
                                         .append("\n")
                                         .append("\n");
                            sstrIndex++;
                        }
                    }
                }
            }
        }

        System.out.println(stringBuilder);


        return !stringBuilder.toString()
                             .isEmpty()
                && FileSystem.writeFile(pathToNeighborsFile, stringBuilder.toString());
    }


    public static boolean writeFragmentFile(final String pathToFragmentsFile, final DataSet fragmentDataSet) {
        final StringBuilder stringBuilder = new StringBuilder();
        final IAtomContainer fragment = fragmentDataSet.getStructure()
                                                       .toAtomContainer();
        IAtom atom;
        IBond bond;
        final Map<IAtom, Integer> sstrMap = new HashMap<>();
        for (int i = 0; i
                < fragment.getAtomCount(); i++) {
            atom = fragment.getAtom(i);
            sstrMap.put(atom, sstrMap.size()
                    + 1);
            if (atom.getSymbol()
                    .equals("R")) {
                stringBuilder.append("SSTR S")
                             .append(sstrMap.size())
                             .append(" A (1 2 3) (0 1 2 3) ")
                             .append("\n");
            } else {
                stringBuilder.append("SSTR S")
                             .append(sstrMap.size())
                             .append(" ")
                             .append(atom.getSymbol())
                             .append(" ")
                             .append(atom.getSymbol()
                                         .equals("C")
                                     ? Constants.hybridizationConversionMap.get(atom.getHybridization()
                                                                                    .name())
                                     : "(1 2 3)")
                             .append(" ")
                             .append(atom.getImplicitHydrogenCount())
                             .append("\n");
            }
        }
        for (int i = 0; i
                < fragment.getBondCount(); i++) {
            bond = fragment.getBond(i);
            stringBuilder.append("LINK S")
                         .append(sstrMap.get(bond.getBegin()))
                         .append(" S")
                         .append(sstrMap.get(bond.getEnd()))
                         .append("\n");
        }
        System.out.println(stringBuilder);

        return !stringBuilder.toString()
                             .isEmpty()
                && FileSystem.writeFile(pathToFragmentsFile, stringBuilder.toString());
    }
}
