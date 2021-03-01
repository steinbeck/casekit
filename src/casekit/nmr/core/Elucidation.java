package casekit.nmr.core;

import casekit.nmr.model.Assignment;
import casekit.nmr.model.DataSet;
import casekit.nmr.model.Spectrum;
import casekit.nmr.utils.Match;

import java.util.ArrayList;
import java.util.List;

public class Elucidation {

    public static List<DataSet> findFragments(final Spectrum querySpectrum, final List<DataSet> compoundDataSets,
                                              final double shiftTol) {
        final List<DataSet> fragments = new ArrayList<>();

        Assignment matchAssignment;
        for (final DataSet dataSet : compoundDataSets) {
            matchAssignment = Match.matchSpectra(dataSet.getSpectrum(), querySpectrum, 0, 0, shiftTol);

        }

        return fragments;
    }

    public static List<DataSet> elucidate() {
        final List<DataSet> solutions = new ArrayList<>();

        return solutions;
    }

}
