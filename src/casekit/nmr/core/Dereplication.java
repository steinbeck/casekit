package casekit.nmr.core;

import casekit.nmr.model.Assignment;
import casekit.nmr.model.DataSet;
import casekit.nmr.model.Spectrum;
import casekit.nmr.utils.Match;
import org.openscience.cdk.exception.CDKException;

import java.util.ArrayList;
import java.util.List;

public class Dereplication {

    public static List<DataSet> dereplicate1D(final Spectrum querySpectrum, final List<DataSet> compoundDataSets,
                                              final double shiftTol) {
        final List<DataSet> solutions = new ArrayList<>();

        for (final DataSet dataSet : compoundDataSets) {
            final Assignment matchAssignment = Match.matchSpectra(dataSet.getSpectrum(), querySpectrum, 0, 0, 1);
            if (matchAssignment.isFullyAssigned(0)) {
                try {
                    dataSet.addMetaInfo("tanimoto", String.valueOf(
                            Match.calculateTanimotoCoefficient(dataSet.getSpectrum(), querySpectrum, 0, 0)));
                } catch (CDKException e) {
                    e.printStackTrace();
                }
                dataSet.addMetaInfo("avgDev", String.valueOf(
                        Match.calculateAverageDeviation(dataSet.getSpectrum(), querySpectrum, 0, 0, shiftTol)));
                solutions.add(dataSet);
            }
        }

        return solutions;
    }
}
