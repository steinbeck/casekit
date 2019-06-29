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
import java.util.ArrayList;
import java.util.HashMap;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Predict {
    
    /**
     * Predicts a shift value for a central atom based on its HOSE code and a 
     * given HOSE code lookup table. The prediction is done by using the mean
     * of all occurring shifts in lookup table for the given HOSE code.
     *
     * @param HOSECodeLookupTable
     * @param HOSECode
     * @return null if HOSE code does not exist in lookup table
     *
     * @see casekit.NMR.Utils#getRMS(java.util.ArrayList)
     *
     */
    public static Double predictShift(final HashMap<String, ArrayList<Double>> HOSECodeLookupTable, final String HOSECode) {
        if (HOSECodeLookupTable.containsKey(HOSECode)) {
            return Utils.getMean(HOSECodeLookupTable.get(HOSECode));
        }

        return null;
    }

    /**
     * Specified for carbons (13C) only. Not generic at the moment because of
     * usage of {@link casekit.NMR.Utils#getMultiplicityFromHydrogenCount(int)}
     * with {@code hCount}.
     *
     * @param HOSECodeLookupTable
     * @param ac
     * @param atomIndex
     * @param maxSphere
     * @param nucleus
     * @param hCount
     * @return null if HOSE code of selected atom does not exist in lookup table
     * 
     * @throws CDKException
     *
     * @see #predictShift(java.util.HashMap, java.lang.String)
     *
     */
    public static Signal predictSignal(final HashMap<String, ArrayList<Double>> HOSECodeLookupTable, final IAtomContainer ac, final int atomIndex, final Integer maxSphere, final String nucleus, final Integer hCount) throws CDKException {
        if (!Utils.checkIndexInAtomContainer(ac, atomIndex) || (hCount == null)) {
            return null;
        }
        final Double predictedShift = Predict.predictShift(HOSECodeLookupTable, HOSECodeBuilder.buildHOSECode(ac, atomIndex, maxSphere, false));
        if (predictedShift == null) {
            return null;
        }
        return new Signal(
                new String[]{nucleus},
                new Double[]{predictedShift},
                Utils.getMultiplicityFromHydrogenCount(hCount),
                null
        );
    }

    /**
     * Specified for carbons (13C) only. Not generic at the moment because of
     * {@link casekit.NMR.Utils#getMultiplicityFromHydrogenCount(int)}.
     *
     * @param HOSECodeLookupTable
     * @param ac
     * @param maxSphere
     * @param nucleus
     * @return null if a HOSE code of one atom does not exist in lookup table
     * 
     * @throws org.openscience.cdk.exception.CDKException
     *
     * @see #predictSignal(java.util.HashMap,
     * org.openscience.cdk.interfaces.IAtomContainer, int, java.lang.Integer,
     * java.lang.String, java.lang.Integer)
     *
     */
    public static Spectrum predictSpectrum(final HashMap<String, ArrayList<Double>> HOSECodeLookupTable, final IAtomContainer ac, final Integer maxSphere, final String nucleus) throws CDKException {
        final Spectrum predictedSpectrum = new Spectrum(new String[]{nucleus});
        Signal signal;
        for (final IAtom atom : ac.atoms()) {            
            if (atom.getSymbol().equals(Utils.getAtomTypeFromSpectrum(predictedSpectrum, 0))) {
                signal = Predict.predictSignal(HOSECodeLookupTable, ac, atom.getIndex(), maxSphere, nucleus, atom.getImplicitHydrogenCount());
                if(signal == null){
                    return null;
                }
                predictedSpectrum.addSignal(signal);
            }
        }

        return predictedSpectrum;
    }
}
