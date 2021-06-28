package casekit.nmr.lsd;

import casekit.nmr.model.Assignment;
import casekit.nmr.model.DataSet;
import casekit.nmr.model.Signal;
import casekit.nmr.model.Spectrum;
import casekit.nmr.utils.Match;
import casekit.nmr.utils.Statistics;
import casekit.nmr.utils.Utils;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;

public class RankedResultSDFParser {

    public static List<DataSet> parseRankedResultSDFile(final String pathToFile,
                                                        final String nucleus) throws CDKException, FileNotFoundException {
        return parseRankedResultSDFile(new FileReader(pathToFile), nucleus);
    }

    public static List<DataSet> parseRankedResultSDFileContent(final String fileContent,
                                                               final String nucleus) throws CDKException {
        final InputStream inputStream = new ByteArrayInputStream(fileContent.getBytes(StandardCharsets.UTF_8));
        return parseRankedResultSDFile(new InputStreamReader(inputStream), nucleus);
    }

    public static List<DataSet> parseRankedResultSDFile(final Reader fileReader,
                                                        final String nucleus) throws CDKException {
        final List<DataSet> dataSetList = new ArrayList<>();
        final IteratingSDFReader iterator = new IteratingSDFReader(fileReader, SilentChemObjectBuilder.getInstance());
        IAtomContainer structure;
        Spectrum experimentalSpectrum, predictedSpectrum;
        Assignment assignment;
        HashMap<String, String> meta;
        final CDKHydrogenAdder hydrogenAdder = CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance());
        IMolecularFormula mf;
        LinkedHashMap<String, String> shiftProperties1D;
        String[] split;
        Signal experimentalSignal;
        double experimentalShift, predictedShift;
        String multiplicity;
        Map<Integer, List<Double>> signalShiftList;
        DataSet dataSet;
        Double[] deviations;
        int signalCounter, matchedSignalIndex;
        List<Integer> closestSignalList;

        while (iterator.hasNext()) {
            structure = iterator.next();
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(structure);
            hydrogenAdder.addImplicitHydrogens(structure);
            Utils.setAromaticityAndKekulize(structure);
            meta = new HashMap<>();
            meta.put("title", structure.getTitle());
            meta.put("id", structure.getProperty("nmrshiftdb2 ID"));
            mf = Utils.getMolecularFormulaFromAtomContainer(structure);
            meta.put("mf", Utils.molecularFormularToString(mf));
            try {
                final String smiles = casekit.nmr.utils.Utils.getSmilesFromAtomContainer(structure);
                meta.put("smiles", smiles);
            } catch (final CDKException e) {
                e.printStackTrace();
            }
            shiftProperties1D = getShiftProperties1D(structure, Utils.getAtomTypeFromNucleus(nucleus));

            experimentalSpectrum = new Spectrum();
            experimentalSpectrum.setNuclei(new String[]{nucleus});
            experimentalSpectrum.setSignals(new ArrayList<>());


            for (final Map.Entry<String, String> shiftProperty1D : shiftProperties1D.entrySet()) {
                split = shiftProperty1D.getValue()
                                       .split("\\s");
                multiplicity = Utils.getMultiplicityFromProtonsCount(structure.getAtom(Integer.parseInt(split[0])
                                                                                               - 1)
                                                                              .getImplicitHydrogenCount());
                experimentalSignal = new Signal();
                experimentalSignal.setNuclei(new String[]{nucleus});
                experimentalSignal.setShifts(new Double[experimentalSignal.getNDim()]);
                experimentalShift = Double.parseDouble(split[1]); // exp. shift
                experimentalSignal.setShift(experimentalShift, 0);
                experimentalSignal.setEquivalencesCount(1);
                experimentalSignal.setMultiplicity(multiplicity);
                experimentalSpectrum.addSignal(experimentalSignal);
            }
            assignment = new Assignment();
            assignment.setNuclei(new String[]{nucleus});
            assignment.initAssignments(experimentalSpectrum.getSignalCount());

            predictedSpectrum = experimentalSpectrum.buildClone();

            deviations = new Double[predictedSpectrum.getSignalCountWithEquivalences()];
            signalCounter = 0;
            signalShiftList = new HashMap<>();
            for (final Map.Entry<String, String> shiftProperty1D : shiftProperties1D.entrySet()) {
                split = shiftProperty1D.getValue()
                                       .split("\\s");
                multiplicity = Utils.getMultiplicityFromProtonsCount(structure.getAtom(Integer.parseInt(split[0])
                                                                                               - 1)
                                                                              .getImplicitHydrogenCount());
                experimentalShift = Double.parseDouble(split[1]); // exp. shift
                predictedShift = Double.parseDouble(split[3]); // pred. shift

                // just to be sure that we take the right signal if equivalences are present
                closestSignalList = experimentalSpectrum.pickByClosestShift(experimentalShift, 0, 0.0);
                closestSignalList.retainAll(experimentalSpectrum.pickByMultiplicity(multiplicity));
                matchedSignalIndex = closestSignalList.get(0);

                deviations[signalCounter] = Math.abs(predictedShift
                                                             - experimentalShift);
                signalShiftList.putIfAbsent(matchedSignalIndex, new ArrayList<>());
                signalShiftList.get(matchedSignalIndex)
                               .add(predictedShift);
                assignment.addAssignmentEquivalence(0, matchedSignalIndex, Integer.parseInt(split[0])
                        - 1);
                signalCounter++;
            }
            for (final int signalIndex : signalShiftList.keySet()) {
                predictedSpectrum.getSignal(signalIndex)
                                 .setShift(Statistics.getMedian(signalShiftList.get(signalIndex)), 0);
                predictedSpectrum.getSignal(signalIndex)
                                 .setEquivalencesCount(signalShiftList.get(signalIndex)
                                                                      .size());
            }
            // if no spectrum could be built or the number of signals in spectrum is different than the atom number in molecule
            if (Utils.getDifferenceSpectrumSizeAndMolecularFormulaCount(predictedSpectrum, mf, 0)
                    != 0) {
                continue;
            }
            dataSet = new DataSet(structure, predictedSpectrum, assignment, meta);
            dataSet.addMetaInfo("rmsd", String.valueOf(Match.calculateRMSD(deviations)));
            dataSet.addMetaInfo("tanimoto", String.valueOf(
                    Match.calculateTanimotoCoefficient(dataSet.getSpectrum(), experimentalSpectrum, 0, 0)));

            dataSetList.add(dataSet);
        }
        // pre-sort by RMSD value
        dataSetList.sort((dataSet1, dataSet2) -> {
            if (Double.parseDouble(dataSet1.getMeta()
                                           .get("rmsd"))
                    < Double.parseDouble(dataSet2.getMeta()
                                                 .get("rmsd"))) {
                return -1;
            } else if (Double.parseDouble(dataSet1.getMeta()
                                                  .get("rmsd"))
                    > Double.parseDouble(dataSet2.getMeta()
                                                 .get("rmsd"))) {
                return 1;
            }
            return 0;
        });

        return dataSetList;
    }

    public static LinkedHashMap<String, String> getShiftProperties1D(final IAtomContainer ac, final String atomType) {
        final LinkedHashMap<String, String> shiftProperties1D = new LinkedHashMap<>();
        String[] split;
        for (final Object key : ac.getProperties()
                                  .keySet()) {
            if (key instanceof String
                    && ((String) key).startsWith("CS")) {
                split = ((String) key).split("CS");
                if (ac.getAtom(Integer.parseInt(split[1])
                                       - 1)
                      .getSymbol()
                      .equals(atomType)) {
                    shiftProperties1D.put(split[1], ac.getProperty(key));
                }
            }
        }

        return shiftProperties1D;
    }
}
