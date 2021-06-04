package casekit.nmr.analysis;

import casekit.nmr.hose.HOSECodeBuilder;
import casekit.nmr.model.DataSet;
import casekit.nmr.model.Signal;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.*;

public class HOSECodeShiftStatisticsBuilder {

    public static Map<String, Map<String, Double[]>> buildHOSECodeShiftStatistics(final List<DataSet> dataSetList,
                                                                                  final int maxSphere) {
        final Map<String, Map<String, List<Double>>> hoseCodeShifts = new HashMap<>();
        IAtomContainer structure;
        Signal signal;
        String hoseCode;
        String solvent;
        for (final DataSet dataSet : dataSetList) {
            structure = dataSet.getStructure()
                               .toAtomContainer();
            solvent = dataSet.getSpectrum()
                             .getSolvent();
            for (int i = 0; i
                    < structure.getAtomCount(); i++) {
                signal = dataSet.getSpectrum()
                                .getSignal(dataSet.getAssignment()
                                                  .getIndex(0, i));
                if (signal
                        != null) {
                    try {
                        for (int sphere = 1; sphere
                                <= maxSphere; sphere++) {
                            hoseCode = HOSECodeBuilder.buildHOSECode(structure, i, sphere, false);
                            hoseCodeShifts.putIfAbsent(hoseCode, new HashMap<>());
                            hoseCodeShifts.get(hoseCode)
                                          .putIfAbsent(solvent, new ArrayList<>());
                            hoseCodeShifts.get(hoseCode)
                                          .get(solvent)
                                          .add(signal.getShift(0));
                        }
                    } catch (final CDKException e) {
                        e.printStackTrace();
                    }
                }
            }
        }
        final Map<String, Map<String, Double[]>> hoseCodeShiftStatistics = new HashMap<>();
        List<Double> values;
        for (final Map.Entry<String, Map<String, List<Double>>> hoseCodes : hoseCodeShifts.entrySet()) {
            hoseCodeShiftStatistics.put(hoseCodes.getKey(), new HashMap<>());
            for (final Map.Entry<String, List<Double>> solvents : hoseCodes.getValue()
                                                                           .entrySet()) {
                values = solvents.getValue(); //casekit.nmr.Utils.removeOutliers(solvents.getValue(), 1.5);
                hoseCodeShiftStatistics.get(hoseCodes.getKey())
                                       .put(solvents.getKey(),
                                            new Double[]{Double.valueOf(values.size()), Collections.min(values),
                                                         casekit.nmr.Utils.getMean(values),
                                                         // casekit.nmr.Utils.getRMS(values),
                                                         casekit.nmr.Utils.getMedian(values), Collections.max(values)});
            }
        }

        return hoseCodeShiftStatistics;
    }
}
