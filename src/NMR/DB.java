/*
 * The MIT License
 *
 * Copyright 2018 Michael Wenk [https://github.com/michaelwenk].
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
package NMR;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.Fingerprinter;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.qsar.descriptors.atomic.AtomHybridizationDescriptor;
import org.openscience.cdk.silent.AtomContainerSet;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class DB {
    
    /**
     * Returns the molecules of a given NMRShiftDB file.
     * This function sets the molecule aromaticity (with allowed exocyclic pi
     * bonds) by using the
     * {@link #setAromaticitiesInAtomContainer(org.openscience.cdk.interfaces.IAtomContainer, int)}
     * function.
     *
     * @param pathToNMRShiftDB path to NMRShiftDB file
     * @param maxCycleSize maximum cycle size for setting the aromaticity in a
     * molecule
     * @return
     * @throws FileNotFoundException
     * @throws CDKException
     */
    public static IAtomContainerSet getStructuresFromNMRShiftDBFile(final String pathToNMRShiftDB, final int maxCycleSize) throws FileNotFoundException, CDKException {

        final IAtomContainerSet acSet = new AtomContainerSet();
        final IteratingSDFReader iterator = new IteratingSDFReader(
                new FileReader(pathToNMRShiftDB),
                SilentChemObjectBuilder.getInstance()
        );
        while (iterator.hasNext()) {
            acSet.addAtomContainer(NMR.Utils.setAromaticitiesInAtomContainer(iterator.next(), maxCycleSize));
        }

        return acSet;
    }

    public static Connection getDBConnection(final String server, final String options, final String user, final String pwd) throws SQLException {

        return DriverManager.getConnection(server + "?" + options, user, pwd);
    }

    public static int[][] countNeighborhoodBonds(final Connection DBConnection, final String[] bondsSet, final String elem, String[] neighborElems, final int minShift, final int maxShift, final int stepSize) throws FileNotFoundException, IOException, SQLException {

        if (DBConnection == null || stepSize < 1) {
            return null;
        }
        // creation of frequency counting matrix and shift indices holder
        final int[][] neighborhoodCountsMatrix = new int[stepSize * (maxShift - minShift + 1)][3 + 4 + neighborElems.length * bondsSet.length];
        HashMap<Integer, ArrayList<Integer>> signalAtomIndicesInNMRShiftDB = new HashMap<>(); // holding of all indices of each ac set (DB) entry [first value] and it's atom indices [second value] too 
        for (int i = 0; i < stepSize * maxShift; i++) {
            for (int k = 0; k < 3 + 4 + neighborElems.length * bondsSet.length; k++) {
                neighborhoodCountsMatrix[i][k] = 0;
            }
            signalAtomIndicesInNMRShiftDB.put(i, new ArrayList<>());
        }

        final Statement statement = DBConnection.createStatement();
        String multQuery = "SELECT (FLOOR(sh.VALUE*" + stepSize + ")/" + stepSize + ") AS shift, nmrsig.MULTIPLICITY AS mult, COUNT(FLOOR(sh.VALUE*" + stepSize + ")/" + stepSize + ") AS shiftCount \n"
                + "FROM         SHIFT AS sh, SIGNAL_ATOM AS sigatom, NMR_SIGNAL AS nmrsig, SPECTRUM AS spec, SPECTRUM_TYPE AS spectype \n"
                + "WHERE 	sh.SIGNAL_ID = sigatom.SIGNAL_ID AND \n"
                + "		sigatom.SIGNAL_ID = nmrsig.SIGNAL_ID AND \n"
                + "		nmrsig.SPECTRUM_ID = spec.SPECTRUM_ID AND \n"
                + "		spec.REVIEW_FLAG = \"true\" AND \n"
                + "		spectype.SPECTRUM_TYPE_ID = spec.SPECTRUM_TYPE_ID AND \n"
                + "		spectype.NAME = \"" + NMR.Utils.getNMRIsotopeIdentifier(elem) + "\" AND \n"
                + "		nmrsig.MULTIPLICITY IS NOT NULL AND \n"
                + "		nmrsig.MULTIPLICITY != \"\""
                + "GROUP BY     shift, mult \n"
                + "HAVING       shift >= " + minShift + " AND shift <= " + maxShift + ";";

        double shiftDouble;
        int shiftInt;
        System.out.println("\n\nneighborhoods:\nQUERY: " + multQuery);
        final ResultSet resultSet = statement.executeQuery(multQuery);
        while (resultSet.next()) {
            shiftDouble = Math.floor(resultSet.getDouble("shift") * stepSize) / (double) stepSize;
            if (shiftDouble < minShift || shiftDouble > maxShift - 1) {
                continue;
            }
            shiftInt = (int) (shiftDouble * stepSize);
            neighborhoodCountsMatrix[shiftInt - minShift][0] += resultSet.getInt("shiftCount");
            switch (resultSet.getString("mult")) {
                case "S": // for qC
                    neighborhoodCountsMatrix[shiftInt - minShift][3] += resultSet.getInt("shiftCount");
                    break;
                case "D": // for CH
                    neighborhoodCountsMatrix[shiftInt - minShift][4] += resultSet.getInt("shiftCount");
                    break;
                case "T": // for CH2
                    neighborhoodCountsMatrix[shiftInt - minShift][5] += resultSet.getInt("shiftCount");
                    break;
                case "Q": // for CH3
                    neighborhoodCountsMatrix[shiftInt - minShift][6] += resultSet.getInt("shiftCount");
                    break;
            }

        }

//        this.neighborhoodCountsMatrix[shiftDB - min][0] += 1; // increase number of this shift occurence
//        this.neighborhoodCountsMatrix[shiftDB - min][1] += (acDB.getAtom(atomIndexDB).isInRing()) ? 1 : 0; // increase if atom is a ring member
//        this.neighborhoodCountsMatrix[shiftDB - min][2] += ((acDB.getAtom(atomIndexDB).getImplicitHydrogenCount() != null) && (acDB.getAtom(atomIndexDB).getImplicitHydrogenCount() == 0)) ? 1 : 0; // qC count or equivalents, e.g. qN
//        this.neighborhoodCountsMatrix[shiftDB - min][3] += ((acDB.getAtom(atomIndexDB).getImplicitHydrogenCount() != null) && (acDB.getAtom(atomIndexDB).getImplicitHydrogenCount() == 1)) ? 1 : 0; // CH count or equivalents, e.g. NH
//        this.neighborhoodCountsMatrix[shiftDB - min][4] += ((acDB.getAtom(atomIndexDB).getImplicitHydrogenCount() != null) && (acDB.getAtom(atomIndexDB).getImplicitHydrogenCount() == 2)) ? 1 : 0; // CH2 count or equivalents, e.g. NH2
//        this.neighborhoodCountsMatrix[shiftDB - min][5] += ((acDB.getAtom(atomIndexDB).getImplicitHydrogenCount() != null) && (acDB.getAtom(atomIndexDB).getImplicitHydrogenCount() == 3)) ? 1 : 0; // CH3 count or equivalents, e.g. NH3
//        // add counts for a specific atom to matrix m
//        int[] counts = NMR.Utils.getNeighborhoodBondsCount(acDB, atomIndexDB, bondsSet, neighborElems);
//        for (int i = 0; i < counts.length; i++) {
//            this.neighborhoodCountsMatrix[shiftDB - min][2 + 4 + i] += counts[i];
//        }
//        // add this atom container index and atom index within it to belonging hash map
//        this.shiftIndicesInACSet.get(shiftDB).add(new Integer[]{k, atomIndexDB});
        return neighborhoodCountsMatrix;
    }
    
    
    public static ResultSet getResultSet(final Connection DBConnection, final String query) throws SQLException{
        
        if (DBConnection == null) {
            return null;
        }
        
        return DBConnection.createStatement().executeQuery(query);
    }
    

    /**
     *
     * @param DBConnection
     * @param minShift
     * @param maxShift
     * @param mult
     * @param minIntens
     * @param maxIntens
     * @param elem
     * @return
     * @throws SQLException
     */
    public static ArrayList<Integer> getSignalIDsFromNMRShiftDB(final Connection DBConnection, final double minShift, final double maxShift, final String mult, final Double minIntens, final Double maxIntens, final String elem) throws SQLException {

        final ArrayList<Integer> spectraIDs = new ArrayList<>();        
        String query = "SELECT nmrsig.SIGNAL_ID AS sigID"
                + " FROM SHIFT AS sh, NMR_SIGNAL AS nmrsig, SPECTRUM AS spec, SPECTRUM_TYPE AS spectype \n"
                + " WHERE 	sh.VALUE >= " + minShift + " AND sh.VALUE <= " + maxShift + " AND \n" // for filtering by means of shift values
                + "		sh.SIGNAL_ID = nmrsig.SIGNAL_ID AND \n"
                + "             nmrsig.SPECTRUM_ID = spec.SPECTRUM_ID AND \n"
                + "             spec.REVIEW_FLAG = \"true\" AND \n" // checks whether review flag is set to true
                + "		spectype.SPECTRUM_TYPE_ID = spec.SPECTRUM_TYPE_ID AND \n"
                + "		spectype.NAME = \"" + NMR.Utils.getNMRIsotopeIdentifier(elem) + "\" \n";
        if(mult != null && !mult.trim().isEmpty()){
            query += "		AND nmrsig.MULTIPLICITY = \"" + mult + "\" \n";
        } else {
            query += "          AND nmrsig.MULTIPLICITY IS NOT NULL AND nmrsig.MULTIPLICITY != \"\" \n";
        }
        if((minIntens != null && minIntens > 0.0) && (maxIntens != null && maxIntens > 0.0)){
            query += "		AND nmrsig.INTENSITY >= " + minIntens + " AND nmrsig.INTENSITY <= " + maxIntens + " \n";
        }
        query += " ;";
        System.out.println("\n\ngetSpectraIDs:\nQUERY: " + query);
        final ResultSet resultSet = NMR.DB.getResultSet(DBConnection, query);
        while (resultSet.next()) {
            spectraIDs.add(resultSet.getInt("sigID"));
        }

        return spectraIDs;
    }
    
    
    public static HashMap<Integer, ArrayList<Integer>> matchSpectrumAgainstDB(final Connection DBConnection, final ArrayList<NMR.Signal> spectrum, final double shiftDev, final Double intensDev, final int stepSize) throws SQLException{
       
        final HashMap<Integer, ArrayList<Integer>> hits = new HashMap<>(); double shift;
        for (int i = 0; i < spectrum.size(); i++) {
            hits.put(i, new ArrayList<>());
            shift = Math.floor(spectrum.get(i).getShift() * stepSize) / (double) stepSize;
            if(spectrum.get(i).getIntensity() != null){
                hits.get(i).addAll(NMR.DB.getSignalIDsFromNMRShiftDB(DBConnection, shift - shiftDev, shift + shiftDev, spectrum.get(i).getMultiplicity(), spectrum.get(i).getIntensity() - intensDev, spectrum.get(i).getIntensity() + intensDev, spectrum.get(i).getElement()));
            } else {
                hits.get(i).addAll(NMR.DB.getSignalIDsFromNMRShiftDB(DBConnection, shift - shiftDev, shift + shiftDev, spectrum.get(i).getMultiplicity(), spectrum.get(i).getIntensity(), spectrum.get(i).getIntensity(), spectrum.get(i).getElement()));
            }
        }
        
        return hits;
    }

    
    /**
     *
     * @param DBConnection
     * @param elem
     * @return
     * @throws SQLException
     */
    public static HashMap<String, ArrayList<Double>> getLookupTableFromNMRShiftDB(final Connection DBConnection, final String elem) throws SQLException {

        if (DBConnection == null) {
            return null;
        }
        final HashMap<String, ArrayList<Double>> lookup = new HashMap<>();
        final Statement statement = DBConnection.createStatement();
        final String query = "SELECT a.HOSE_CODE AS hose, sh.VALUE AS shift \n"
                + " FROM ATOM AS a, SHIFT AS sh, SIGNAL_ATOM AS sigatom, NMR_SIGNAL AS nmrsig, SPECTRUM AS spec, SPECTRUM_TYPE AS spectype \n"
                + " WHERE 	a.ATOM_ID = sigatom.ATOM_ID AND \n" // to get signals of each atom
                + "		sigatom.SIGNAL_ID = sh.SIGNAL_ID AND \n" // to get shift values
                + "		sigatom.SIGNAL_ID = nmrsig.SIGNAL_ID AND \n"
                + "             nmrsig.SPECTRUM_ID = spec.SPECTRUM_ID AND \n"
                + "             spec.REVIEW_FLAG = \"true\" AND \n" // checks whether review flag is set to true
                + "		spectype.SPECTRUM_TYPE_ID = spec.SPECTRUM_TYPE_ID AND \n"
                + "		spectype.NAME = \"" + NMR.Utils.getNMRIsotopeIdentifier(elem) + "\";";
        System.out.println("\n\ngetLookupTable:\nQUERY: " + query);
        final ResultSet resultSet = statement.executeQuery(query);
        while (resultSet.next()) {
            if (!lookup.containsKey(resultSet.getString("hose"))) {
                lookup.put(resultSet.getString("hose"), new ArrayList<>());
            }
            lookup.get(resultSet.getString("hose")).add(resultSet.getDouble("shift"));
        }

        return lookup;
    }
    
    
    /**
     *
     * @param DBConnection
     * @param minShift
     * @param maxShift
     * @param elem
     * @return
     * @throws SQLException
     */
    public static HashMap<String, Double> getRMS(final Connection DBConnection, final double minShift, final double maxShift, final String elem) throws SQLException {

        if (DBConnection == null) {
            return null;
        }
        final HashMap<String, Double> rms = new HashMap<>();
        final Statement statement = DBConnection.createStatement();
        final String query = "SELECT a.HOSE_CODE AS hose, COUNT(sh.VALUE) AS shiftCount, AVG(sh.VALUE) AS mean, SQRT(SUM(POW(sh.VALUE, 2))/COUNT(sh.VALUE)) AS rms \n"
                + " FROM ATOM AS a, SHIFT AS sh, SIGNAL_ATOM AS sigatom, NMR_SIGNAL AS nmrsig, SPECTRUM AS spec, SPECTRUM_TYPE AS spectype \n"
                + " WHERE 	sh.VALUE >= " + minShift + " AND sh.VALUE <= " + maxShift + " AND \n" // for filtering by means of shift values
                + "		sh.SIGNAL_ID = sigatom.SIGNAL_ID AND \n" // to get shift values
                + "             a.ATOM_ID = sigatom.ATOM_ID AND \n"
                + "		sigatom.SIGNAL_ID = nmrsig.SIGNAL_ID AND \n"
                + "             nmrsig.SPECTRUM_ID = spec.SPECTRUM_ID AND \n"
                + "             spec.REVIEW_FLAG = \"true\" AND \n" // checks whether review flag is set to true
                + "		spectype.SPECTRUM_TYPE_ID = spec.SPECTRUM_TYPE_ID AND \n"
                + "		spectype.NAME = \"" + NMR.Utils.getNMRIsotopeIdentifier(elem) + "\" AND \n"
                + "		nmrsig.MULTIPLICITY IS NOT NULL AND \n"
                + "		nmrsig.MULTIPLICITY != \"\" \n"
                + " GROUP BY hose;";
        System.out.println("\n\nRMS SQL:\nQUERY: " + query);
        final ResultSet resultSet = statement.executeQuery(query);
        while (resultSet.next()) {
            rms.put(resultSet.getString("hose"), resultSet.getDouble("rms"));
        }

        return rms;
    }
    
    
    public static ArrayList<String> getSpectraFromNMRShiftDBEntry(final IAtomContainer ac, final String elem) {

        ArrayList<String> props = (ArrayList<String>) (ArrayList<?>) (new ArrayList<>(ac.getProperties().keySet()));
        final ArrayList<String> spectra = new ArrayList<>();
        for (String prop : props) {
            if (prop.contains("Spectrum " + NMR.Utils.getNMRIsotopeIdentifier(elem))) {
                spectra.add(ac.getProperty(prop));
            }
        }

        return spectra;
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    // TRIALS
    
    
    public static void findSubstructuresInNMRShiftDB(final IAtomContainer acQ, final String pathToNMRShiftDB) throws CDKException, FileNotFoundException, CloneNotSupportedException {

        final IAtomContainer acQcopy = acQ.clone();
        AtomContainerManipulator.convertImplicitToExplicitHydrogens(acQcopy);

        final IFingerprinter fingerprinter = new Fingerprinter();
        final IBitFingerprint fingerprintQ = fingerprinter.getBitFingerprint(acQcopy);
//        System.out.println("Q: cardinality: " + fingerprintQ.cardinality() + ", bit set: " + Arrays.toString(fingerprintQ.getSetbits()));
        IBitFingerprint fingerprintDB;
        IAtomContainer acDB;
        final IteratingSDFReader iterator = new IteratingSDFReader(
                new FileReader(pathToNMRShiftDB),
                SilentChemObjectBuilder.getInstance()
        );
        final AtomHybridizationDescriptor hybridDesc = new AtomHybridizationDescriptor();
        int molCounter = 1;
        while (iterator.hasNext()) {
            acDB = iterator.next();
//            // skip structures which do not at least contain one carbon spectrum
//            if (!acDB.getProperties().containsKey("Spectrum 13C 0")) {
//                continue;
//            }
//            IAtomContainer acDBcopy = acDB.clone();
//            AtomContainerManipulator.convertImplicitToExplicitHydrogens(acDBcopy);
//            fingerprintDB = fingerprinter.getBitFingerprint(acDBcopy);
////            System.out.println("DB: cardinality: " + fingerprintDB.cardinality() + ", bit set: " + Arrays.toString(fingerprintDB.getSetbits()));
//            
////            fingerprintDB.and(fingerprintQ);
////            System.out.println("and: " + Arrays.toString(fingerprintDB.getSetbits()));
//            
//            if(Tanimoto.calculate(fingerprintQ, fingerprintDB) >= 0.25)
//                System.out.println("similarity: " + Tanimoto.calculate(fingerprintQ, fingerprintDB) + " at " + acDB.getProperties());

            int counter = 0;
            for (int i = 0; i < acDB.getAtomCount(); i++) {
                if (acDB.getAtom(i).getSymbol().equals("N")) {
                    for (IAtom neighbor : acDB.getConnectedAtomsList(acDB.getAtom(i))) {
                        if (neighbor.getSymbol().equals("C") && IAtomType.Hybridization.values()[Integer.parseInt(hybridDesc.calculate(neighbor, acDB).getValue().toString())].equals(IAtomType.Hybridization.SP2)) {
                            for (IAtom neighbor2 : acDB.getConnectedAtomsList(neighbor)) {
                                if (neighbor2.getSymbol().equals("C") && neighbor2.getImplicitHydrogenCount() == 3) {
                                    counter++;
                                }
                            }
                        }
                    }
                    if (counter >= 2) {
                        System.out.println(molCounter);
                    }
                    break;
                }
            }
            molCounter++;
        }

    }
}
