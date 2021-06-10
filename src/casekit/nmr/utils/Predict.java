/*
 * The MIT License
 *
 * Copyright 2019 Michael Wenk [https://github.com/michaelwenk].
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package casekit.nmr.utils;


import casekit.nmr.hose.HOSECodeBuilder;
import casekit.nmr.hose.model.ConnectionTree;
import casekit.nmr.hose.model.ConnectionTreeNode;
import casekit.nmr.model.Assignment;
import casekit.nmr.model.DataSet;
import casekit.nmr.model.Signal;
import casekit.nmr.model.Spectrum;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Predict {

    /**
     * Diastereotopic distinctions are not provided yet.
     *
     * @param hoseCodeShiftStatistics
     * @param structure
     * @param solvent
     * @param nucleus
     *
     * @return
     */
    public static DataSet predict1D(final Map<String, Map<String, Double[]>> hoseCodeShiftStatistics,
                                    final IAtomContainer structure, final String solvent, final String nucleus) {
        final int minMatchingSphere = 1;
        final Spectrum spectrum = new Spectrum();
        spectrum.setNuclei(new String[]{nucleus});
        spectrum.setSolvent(solvent);
        spectrum.setSignals(new ArrayList<>());
        final Assignment assignment = new Assignment();
        assignment.setNuclei(spectrum.getNuclei());
        assignment.initAssignments(0);

        final CDKHydrogenAdder hydrogenAdder = CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance());
        String hoseCode, atomTypeSpectrum;
        Signal signal;
        Double shift;
        Integer addedSignalIndex;
        ConnectionTree connectionTree;

        try {
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(structure);
            casekit.nmr.Utils.convertExplicitToImplicitHydrogens(structure);
            hydrogenAdder.addImplicitHydrogens(structure);
            casekit.nmr.Utils.convertImplicitToExplicitHydrogens(structure);
            casekit.nmr.Utils.setAromaticityAndKekulize(structure);

            for (int i = 0; i
                    < structure.getAtomCount(); i++) {
                atomTypeSpectrum = casekit.nmr.utils.Utils.getAtomTypeFromNucleus(nucleus);
                if (structure.getAtom(i)
                             .getSymbol()
                             .equals(atomTypeSpectrum)) {
                    connectionTree = HOSECodeBuilder.buildConnectionTree(structure, i, null);
                    shift = null;
                    for (int s = connectionTree.getMaxSphere(); s
                            >= minMatchingSphere; s--) {
                        hoseCode = HOSECodeBuilder.buildHOSECode(structure, i, s, false);
                        if (hoseCodeShiftStatistics.containsKey(hoseCode)
                                && hoseCodeShiftStatistics.get(hoseCode)
                                                          .containsKey(solvent)) {
                            shift = hoseCodeShiftStatistics.get(hoseCode)
                                                           .get(solvent)[3]; // take median value
                            break;
                        }
                    }
                    signal = new Signal();
                    signal.setNuclei(spectrum.getNuclei());
                    signal.setEquivalencesCount(1);
                    //                    signal.setMultiplicity();
                    signal.setKind("signal");
                    signal.setShifts(new Double[]{shift});
                    addedSignalIndex = spectrum.addSignal(signal);
                    if (addedSignalIndex
                            >= assignment.getSetAssignmentsCount(0)) {
                        assignment.addAssignment(0, new int[]{i});
                    } else {
                        assignment.addAssignmentEquivalence(0, addedSignalIndex, i);
                    }
                }
            }
        } catch (final CDKException e) {
            e.printStackTrace();
            return null;
        }

        return new DataSet(structure, spectrum, assignment, new HashMap<>());
    }

    /**
     * Predicts a 2D spectrum from two 1D spectra. Each 1D spectra needs to contain the same solvent information.
     * Diastereotopic distinctions are not provided yet.
     *
     * @param hoseCodeShiftStatistics HOSE code shift statistics
     * @param structure               structure to use for prediction
     * @param nuclei                  nuclei for 2D spectrum to predict
     * @param solvent                 solvent
     * @param minPathLength           minimal path length
     * @param maxPathLength           maximal path length
     *
     * @return
     */
    public static DataSet predict2D(final Map<String, Map<String, Double[]>> hoseCodeShiftStatistics,
                                    final IAtomContainer structure, final String[] nuclei, final String solvent,
                                    final int minPathLength, final int maxPathLength) {
        final DataSet predictionDim1 = predict1D(hoseCodeShiftStatistics, structure, solvent, nuclei[0]);
        final DataSet predictionDim2 = predict1D(hoseCodeShiftStatistics, structure, solvent, nuclei[1]);
        return Predict.predict2D(structure, predictionDim1.getSpectrum(), predictionDim2.getSpectrum(),
                                 predictionDim1.getAssignment(), predictionDim2.getAssignment(), minPathLength,
                                 maxPathLength);
    }

    /**
     * Predicts a 2D spectrum from two 1D spectra. Each 1D spectra needs to contain the same solvent information.
     * Diastereotopic distinctions are not provided yet.
     *
     * @param structure      structure to use for prediction
     * @param spectrumDim1   1D spectrum of first dimension
     * @param spectrumDim2   1D spectrum of second dimension
     * @param assignmentDim1 1D assignment of first dimension
     * @param assignmentDim2 1D assignment of second dimension
     * @param minPathLength  minimal path length
     * @param maxPathLength  maximal path length
     *
     * @return
     */
    public static DataSet predict2D(final IAtomContainer structure, final Spectrum spectrumDim1,
                                    final Spectrum spectrumDim2, final Assignment assignmentDim1,
                                    final Assignment assignmentDim2, final int minPathLength, final int maxPathLength) {
        if (!spectrumDim1.getSolvent()
                         .equals(spectrumDim2.getSolvent())) {
            return null;
        }
        final String[] nuclei2D = new String[]{spectrumDim1.getNuclei()[0], spectrumDim2.getNuclei()[0]};
        final String atomTypeDim1 = casekit.nmr.utils.Utils.getAtomTypeFromNucleus(spectrumDim1.getNuclei()[0]);
        final String atomTypeDim2 = casekit.nmr.utils.Utils.getAtomTypeFromNucleus(spectrumDim2.getNuclei()[0]);

        final Spectrum predictedSpectrum2D = new Spectrum();
        predictedSpectrum2D.setNuclei(nuclei2D);
        predictedSpectrum2D.setSignals(new ArrayList<>());
        predictedSpectrum2D.setSolvent(spectrumDim1.getSolvent());
        final Assignment assignment2D = new Assignment();
        assignment2D.setNuclei(predictedSpectrum2D.getNuclei());
        assignment2D.initAssignments(0);

        Signal signal2D;
        IAtom atom;
        Double shiftDim1, shiftDim2;
        int addedSignalIndex;
        ConnectionTree connectionTree;
        List<ConnectionTreeNode> nodesInSphere;
        for (int i = 0; i
                < structure.getAtomCount(); i++) {
            atom = structure.getAtom(i);
            if (atom.getSymbol()
                    .equals(atomTypeDim1)) {
                connectionTree = HOSECodeBuilder.buildConnectionTree(structure, i, maxPathLength);
                for (int s = minPathLength; s
                        <= connectionTree.getMaxSphere(); s++) {
                    nodesInSphere = connectionTree.getNodesInSphere(s, false);
                    for (final ConnectionTreeNode nodeInSphere : nodesInSphere) {
                        if (nodeInSphere.getAtom()
                                        .getSymbol()
                                        .equals(atomTypeDim2)) {
                            signal2D = new Signal();
                            signal2D.setNuclei(nuclei2D);
                            signal2D.setKind("signal");
                            signal2D.setEquivalencesCount(1);
                            shiftDim1 = spectrumDim1.getShift(assignmentDim1.getIndex(0, i), 0);
                            shiftDim2 = spectrumDim2.getShift(assignmentDim2.getIndex(0, nodeInSphere.getKey()), 0);
                            signal2D.setShifts(new Double[]{shiftDim1, shiftDim2});

                            addedSignalIndex = predictedSpectrum2D.addSignal(signal2D);
                            if (addedSignalIndex
                                    >= assignment2D.getSetAssignmentsCount(0)) {
                                assignment2D.addAssignment(0, new int[]{i});
                                assignment2D.addAssignment(1, new int[]{nodeInSphere.getKey()});
                            } else {
                                assignment2D.addAssignmentEquivalence(0, addedSignalIndex, i);
                                assignment2D.addAssignmentEquivalence(1, addedSignalIndex, nodeInSphere.getKey());
                            }
                        }
                    }
                }
            }
        }

        return new DataSet(structure, predictedSpectrum2D, assignment2D, new HashMap<>());
    }

    public static DataSet predictHSQC(final IAtomContainer structure, final Spectrum spectrumDim1,
                                      final Spectrum spectrumDim2, final Assignment assignmentDim1,
                                      final Assignment assignmentDim2) {
        return predict2D(structure, spectrumDim1, spectrumDim2, assignmentDim1, assignmentDim2, 1, 1);
    }

    public static DataSet predictHSQCEdited(final IAtomContainer structure, final Spectrum spectrumDim1,
                                            final Spectrum spectrumDim2, final Assignment assignmentDim1,
                                            final Assignment assignmentDim2) {
        final DataSet dataSet = predictHSQC(structure, spectrumDim1, spectrumDim2, assignmentDim1, assignmentDim2);

        final String atomTypeDim2 = Utils.getAtomTypeFromSpectrum(spectrumDim2, 0);
        IAtom atom;
        Integer explicitHydrogensCount;
        for (int i = 0; i
                < dataSet.getSpectrum()
                         .getSignalCount(); i++) {
            atom = structure.getAtom(dataSet.getAssignment()
                                            .getAssignment(1, i, 0));
            if (!atom.getSymbol()
                     .equals(atomTypeDim2)) {
                continue;
            }
            explicitHydrogensCount = AtomContainerManipulator.countExplicitHydrogens(structure, atom);
            if (explicitHydrogensCount
                    == 2) {
                dataSet.getSpectrum()
                       .getSignal(i)
                       .setPhase(-1);
            } else if (explicitHydrogensCount
                    == 1
                    || explicitHydrogensCount
                    == 3) {
                dataSet.getSpectrum()
                       .getSignal(i)
                       .setPhase(1);
            }
        }

        return dataSet;
    }
}
