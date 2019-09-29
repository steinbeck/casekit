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
package casekit.NMR.predict;


import casekit.NMR.Utils;
import casekit.NMR.model.Signal;
import casekit.NMR.model.Spectrum;
import hose.HOSECodeBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Predict {
    
    /**
     * Predicts a shift value for a central atom based on its HOSE code and a 
     * given HOSE code lookup table. The prediction is done by using the median
     * of all occurring shifts in lookup table for the given HOSE code. <br>
     * Specified for carbons (13C) only -> {@link casekit.NMR.Utils#getMultiplicityFromHydrogenCount(int)}.
     *
     * @param HOSECodeLookupTable HashMap containing HOSE codes as keys and a list of chemical shifts
     *                            of occurring central atoms as values
     * @param HOSECode specific HOSE code to use for shift prediction
     * @return null if HOSE code does not exist in lookup table
     *
     * @see casekit.NMR.Utils#getMedian(ArrayList)
     *
     */
    public static Double predictShift(final HashMap<String, ArrayList<Double>> HOSECodeLookupTable, final String HOSECode) {
        if (HOSECodeLookupTable.containsKey(HOSECode)) {
            return Utils.getMedian(HOSECodeLookupTable.get(HOSECode));
//            return Utils.getMean(HOSECodeLookupTable.get(HOSECode));
        }

        return null;
    }

    /**
     * Predicts a signal for a central atom based on its HOSE code and a
     * given HOSE code lookup table. The prediction is done by using the mean
     * of all occurring shifts in lookup table for the given HOSE code. <br>
     * Specified for carbons (13C) only -> {@link casekit.NMR.Utils#getMultiplicityFromHydrogenCount(int)}.
     *
     * @param HOSECodeLookupTable HashMap containing HOSE codes as keys and a list of chemical shifts
     *                            of occurring central atoms as values
     * @param ac structure to predict from
     * @param atomIndex index of central atom in structure for HOSE code generation
     * @param maxSphere maximum sphere to use for HOSE code generation or null for unlimited
     * @param nucleus nucleus (e.g. "13C") for signal creation
     *
     * @return null if HOSE code of selected atom does not exist in lookup table
     * 
     * @throws CDKException
     *
     * @see #predictShift(HashMap, String)
     *
     */
    public static Signal predictSignal(final HashMap<String, ArrayList<Double>> HOSECodeLookupTable, final IAtomContainer ac, final int atomIndex, final Integer maxSphere, final String nucleus) throws Exception {
        if (!Utils.checkIndexInAtomContainer(ac, atomIndex)) {
            return null;
        }
        final String HOSECode = HOSECodeBuilder.buildHOSECode(ac, atomIndex, maxSphere, false);
        final Double predictedShift = Predict.predictShift(HOSECodeLookupTable, HOSECode);
        if (predictedShift == null) {
            return null;
        }
        return new Signal(
                new String[]{nucleus},
                new Double[]{predictedShift},
                Utils.getMultiplicityFromHydrogenCount(ac.getAtom(atomIndex).getImplicitHydrogenCount()),
                null
        );
    }

    /**
     * Predicts a spectrum for a given structure based on HOSE code of atoms with specified nucleus and a
     * given HOSE code lookup table. <br>
     * Specified for carbons (13C) only -> {@link casekit.NMR.Utils#getMultiplicityFromHydrogenCount(int)}.
     *
     * @param HOSECodeLookupTable HashMap containing HOSE codes as keys and a list of chemical shifts
     *                            of occurring central atoms as values
     * @param ac structure to predict from
     * @param maxSphere maximum sphere to use for HOSE code generation or null for unlimited
     * @param nucleus nucleus (e.g. "13C") for signal creation
     * @return null if a HOSE code of one atom does not exist in lookup table
     * 
     * @throws org.openscience.cdk.exception.CDKException
     *
     * @see #predictSignal(HashMap, IAtomContainer, int, Integer, String)
     *
     */
    public static Spectrum predictSpectrum(final HashMap<String, ArrayList<Double>> HOSECodeLookupTable, final IAtomContainer ac, final Integer maxSphere, final String nucleus) throws Exception {
        final Spectrum predictedSpectrum = new Spectrum(new String[]{nucleus});
        Signal signal;
        for (final IAtom atom : ac.atoms()) {            
            if (atom.getSymbol().equals(Utils.getAtomTypeFromSpectrum(predictedSpectrum, 0))) {
                signal = Predict.predictSignal(HOSECodeLookupTable, ac, atom.getIndex(), maxSphere, nucleus);
                if(signal == null){
                    continue;
//                    return null;
                }
                predictedSpectrum.addSignal(signal);
            }
        }

        return predictedSpectrum;
    }
}
