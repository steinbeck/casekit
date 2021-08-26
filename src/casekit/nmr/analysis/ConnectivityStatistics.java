package casekit.nmr.analysis;

import casekit.nmr.model.DataSet;
import casekit.nmr.model.Spectrum;
import casekit.nmr.utils.Utils;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class ConnectivityStatistics {


    /**
     * @param dataSetList
     * @param nucleus
     * @param connectivityStatistics multiplicity -> hybridization -> shift (int) -> connected atom symbol -> connected atom hybridization -> connected atom protons count -> occurrence
     */
    public static void buildConnectivityStatistics(final List<DataSet> dataSetList, final String nucleus,
                                                   final Map<String, Map<String, Map<Integer, Map<String, Map<String, Map<Integer, Integer>>>>>> connectivityStatistics) {
        IAtomContainer structure;
        Spectrum spectrum;
        IAtom atom;
        final String atomType = Utils.getAtomTypeFromNucleus(nucleus);
        String multiplicity;
        String hybridization;
        String connectedAtomType;
        String connectedAtomHybridization;
        int shift, atomIndex;
        for (final DataSet dataSet : dataSetList) {
            if (!dataSet.getSpectrum()
                        .getNuclei()[0].equals(nucleus)) {
                continue;
            }
            structure = dataSet.getStructure()
                               .toAtomContainer();
            spectrum = dataSet.getSpectrum()
                              .toSpectrum();
            for (int signalIndex = 0; signalIndex
                    < spectrum.getSignalCount(); signalIndex++) {
                shift = spectrum.getShift(signalIndex, 0)
                                .intValue();
                for (int equivalenceIndex = 0; equivalenceIndex
                        < dataSet.getAssignment()
                                 .getAssignment(0, signalIndex).length; equivalenceIndex++) {
                    atomIndex = dataSet.getAssignment()
                                       .getAssignment(0, signalIndex, equivalenceIndex);
                    atom = structure.getAtom(atomIndex);
                    if (atom.getSymbol()
                            .equals(atomType)) {
                        multiplicity = Utils.getMultiplicityFromProtonsCount(atom.getImplicitHydrogenCount());
                        if (multiplicity
                                == null) {
                            continue;
                        }
                        multiplicity = multiplicity.toLowerCase();
                        hybridization = atom.getHybridization()
                                            .name();
                        connectivityStatistics.putIfAbsent(multiplicity, new HashMap<>());
                        connectivityStatistics.get(multiplicity)
                                              .putIfAbsent(hybridization, new HashMap<>());
                        // check for connected hetero atoms
                        for (final IAtom connectedAtom : structure.getConnectedAtomsList(atom)) {
                            if (connectedAtom.getSymbol()
                                             .equals("H")) {
                                continue;
                            }
                            connectedAtomType = connectedAtom.getSymbol();
                            if (connectedAtom.getHybridization()
                                    == null) {
                                continue;
                            }
                            connectedAtomHybridization = connectedAtom.getHybridization()
                                                                      .name();
                            connectivityStatistics.get(multiplicity)
                                                  .get(hybridization)
                                                  .putIfAbsent(shift, new HashMap<>());
                            connectivityStatistics.get(multiplicity)
                                                  .get(hybridization)
                                                  .get(shift)
                                                  .putIfAbsent(connectedAtomType, new HashMap<>());
                            connectivityStatistics.get(multiplicity)
                                                  .get(hybridization)
                                                  .get(shift)
                                                  .get(connectedAtomType)
                                                  .putIfAbsent(connectedAtomHybridization, new HashMap<>());
                            connectivityStatistics.get(multiplicity)
                                                  .get(hybridization)
                                                  .get(shift)
                                                  .get(connectedAtomType)
                                                  .get(connectedAtomHybridization)
                                                  .putIfAbsent(connectedAtom.getImplicitHydrogenCount(), 0);
                            connectivityStatistics.get(multiplicity)
                                                  .get(hybridization)
                                                  .get(shift)
                                                  .get(connectedAtomType)
                                                  .get(connectedAtomHybridization)
                                                  .put(connectedAtom.getImplicitHydrogenCount(),
                                                       connectivityStatistics.get(multiplicity)
                                                                             .get(hybridization)
                                                                             .get(shift)
                                                                             .get(connectedAtomType)
                                                                             .get(connectedAtomHybridization)
                                                                             .get(connectedAtom.getImplicitHydrogenCount())
                                                               + 1);
                        }
                    }
                }
            }
        }
    }

    /**
     * @param connectivityStatistics   multiplicity -> hybridization -> shift (int) -> connected atom symbol -> connected atom hybridization -> connected atom protons count -> occurrence
     * @param multiplicity
     * @param hybridization
     * @param shift
     * @param molecularFormulaElements
     *
     * @return
     */
    public static Map<String, Map<String, Map<Integer, Integer>>> extractConnectivities(
            final Map<String, Map<String, Map<Integer, Map<String, Map<String, Map<Integer, Integer>>>>>> connectivityStatistics,
            final String multiplicity, final String hybridization, final int shift,
            final Set<String> molecularFormulaElements) {
        final Map<String, Map<String, Map<Integer, Integer>>> extractedConnectivities = new HashMap<>();
        if (connectivityStatistics.containsKey(multiplicity)
                && connectivityStatistics.get(multiplicity)
                                         .containsKey(hybridization)
                && connectivityStatistics.get(multiplicity)
                                         .get(hybridization)
                                         .containsKey(shift)) {
            for (final Map.Entry<String, Map<String, Map<Integer, Integer>>> entry : connectivityStatistics.get(
                    multiplicity)
                                                                                                           .get(hybridization)
                                                                                                           .get(shift)
                                                                                                           .entrySet()) {
                if (molecularFormulaElements.contains(entry.getKey())) {
                    extractedConnectivities.put(entry.getKey(), entry.getValue());
                }
            }
        }

        return extractedConnectivities;
    }
}
