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
package casekit.NMR;

import casekit.NMR.model.Assignment;
import casekit.NMR.model.Signal;
import casekit.NMR.model.Spectrum;
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
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.AtomContainerSet;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class DB {
    
    /**
     * Returns the molecules of a given MOL/SDF file.
     * This function sets the molecule aromaticity (with allowed exocyclic pi
     * bonds) by using the
     * {@link Utils#setAromaticitiesInAtomContainer(org.openscience.cdk.interfaces.IAtomContainer)}
     * function.
     *
     * @param pathToNMRShiftDB path to NMRShiftDB file
     * @param setAromaticity whether to set aromaticities in structures or not
     * @return
     * @throws FileNotFoundException
     * @throws CDKException
     * @deprecated 
     */
    public static IAtomContainerSet getStructuresFromSDFile(final String pathToNMRShiftDB, final boolean setAromaticity) throws FileNotFoundException, CDKException {      
        final IAtomContainerSet acSet = new AtomContainerSet();
        final IteratingSDFReader iterator = new IteratingSDFReader(
                new FileReader(pathToNMRShiftDB),
                SilentChemObjectBuilder.getInstance()
        );
        IAtomContainer ac;
        while (iterator.hasNext()) {
            ac = iterator.next();
            if(setAromaticity){
                Utils.setAromaticitiesInAtomContainer(ac);
            }
            acSet.addAtomContainer(ac);
        }

        return acSet;
    }
    
    
    /**
     * Returns the spectra of a given MOL/SDF file containing NMRShiftDB properties.
     *
     * @param pathToNMRShiftDB path to NMRShiftDB file
     * @param NMRShiftDBSpectrumProperty spectrum property name to use
     * @param atomType atomType of requested spectra
     * @return
     * @throws FileNotFoundException
     * @throws CDKException
     * @deprecated 
     */
    public static ArrayList<Spectrum> getSpectraFromSDFile(final String pathToNMRShiftDB, final String NMRShiftDBSpectrumProperty , final String atomType) throws FileNotFoundException, CDKException {

        final ArrayList<Spectrum> spectrumSet = new ArrayList<>();
        final IteratingSDFReader iterator = new IteratingSDFReader(
                new FileReader(pathToNMRShiftDB),
                SilentChemObjectBuilder.getInstance()
        );
        IAtomContainer ac;
        while (iterator.hasNext()) {
            ac = iterator.next();
            if((ac == null) || (ac.getProperty(NMRShiftDBSpectrumProperty) == null)){
                spectrumSet.add(null);
            } else {                
                spectrumSet.add(DB.NMRShiftDBSpectrumToSpectrum(ac.getProperty(NMRShiftDBSpectrumProperty), atomType));
            }
        }

        return spectrumSet;
    }
    
    /**
     * Returns 3-tuples consisting of structure, spectrum and assignments 
     * for each molecule in the given NMRShiftDB file.
     *
     * @param pathToNMRShiftDB path to NMRShiftDB file
     * @param NMRShiftDBSpectrumProperty spectrum property string to use
     * @param atomType atomType of requested spectra
     * @return
     * @throws FileNotFoundException
     * @throws CDKException
     */
    public static HashMap<Integer, Object[]> getSSCComponentsFromNMRShiftDB(final String pathToNMRShiftDB, final String NMRShiftDBSpectrumProperty, final String atomType) throws FileNotFoundException, CDKException {
        final HashMap<Integer, Object[]> structureSetWithSpectra = new HashMap<>();
        final IteratingSDFReader iterator = new IteratingSDFReader(
                new FileReader(pathToNMRShiftDB),
                SilentChemObjectBuilder.getInstance()
        );
        IAtomContainer ac;
        Spectrum spectrum;
        Assignment assignment;
        while (iterator.hasNext()) {
            ac = iterator.next();
            Utils.setAromaticitiesInAtomContainer(ac);
            
            spectrum = DB.NMRShiftDBSpectrumToSpectrum(ac.getProperty(NMRShiftDBSpectrumProperty), atomType);
            assignment = DB.NMRShiftDBSpectrumToAssignment(ac.getProperty(NMRShiftDBSpectrumProperty), atomType);
            if ((ac != null) && (spectrum != null)) {                
                structureSetWithSpectra.put(structureSetWithSpectra.size(), new Object[]{ac, spectrum, assignment});
            }
        }

        return structureSetWithSpectra;
    }

    /**
     *
     * @param pathToDB
     * @return
     * @throws FileNotFoundException
     * @deprecated 
     */
    public static HashSet<String> getAtomTypesInDB(final String pathToDB) throws FileNotFoundException{
        final HashSet<String> atomTypes = new HashSet<>();
        final IteratingSDFReader iterator = new IteratingSDFReader(
                new FileReader(pathToDB),
                SilentChemObjectBuilder.getInstance()
        );
        while (iterator.hasNext()) {
            atomTypes.addAll(Utils.getAtomTypesInAtomContainer(iterator.next()));
        }
        
        return atomTypes;
    }
    
    /**
     *
     * @param server
     * @param options
     * @param user
     * @param pwd
     * @return
     * @throws SQLException
     * @deprecated 
     */
    public static Connection getDBConnection(final String server, final String options, final String user, final String pwd) throws SQLException {

        return DriverManager.getConnection(server + "?" + options, user, pwd);
    }

    /**
     *
     * @param DBConnection
     * @param bondsSet
     * @param elem
     * @param neighborElems
     * @param minShift
     * @param maxShift
     * @param stepSize
     * @return
     * @throws FileNotFoundException
     * @throws IOException
     * @throws SQLException
     * @deprecated 
     */
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
                + "		spectype.NAME = \"" + casekit.NMR.Utils.getIsotopeIdentifier(elem) + "\" AND \n"
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
    
    /**
     *
     * @param DBConnection
     * @param query
     * @return
     * @throws SQLException
     * @deprecated 
     */
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
     * @deprecated 
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
                + "		spectype.NAME = \"" + casekit.NMR.Utils.getIsotopeIdentifier(elem) + "\" \n";
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
        final ResultSet resultSet = casekit.NMR.DB.getResultSet(DBConnection, query);
        while (resultSet.next()) {
            spectraIDs.add(resultSet.getInt("sigID"));
        }

        return spectraIDs;
    }
    
    // currently only for 1D spectra

    /**
     *
     * @param DBConnection
     * @param spectrum
     * @param shiftDev
     * @param intensDev
     * @param stepSize
     * @param dim
     * @return
     * @throws SQLException
     * @deprecated 
     */
    public static HashMap<Integer, ArrayList<Integer>> matchSpectrumAgainstDB(final Connection DBConnection, final Spectrum spectrum, final double shiftDev, final Double intensDev, final int stepSize, final int dim) throws SQLException{
       
        final HashMap<Integer, ArrayList<Integer>> hits = new HashMap<>(); double shift;
        for (int i = 0; i < spectrum.getSignalCount(); i++) {
            hits.put(i, new ArrayList<>());
            shift = Math.floor(spectrum.getSignal(i).getShift(dim) * stepSize) / (double) stepSize;
            if(spectrum.getSignal(i).getIntensity() != null){
                hits.get(i).addAll(casekit.NMR.DB.getSignalIDsFromNMRShiftDB(DBConnection, shift - shiftDev, shift + shiftDev, spectrum.getSignal(i).getMultiplicity(), spectrum.getSignal(i).getIntensity() - intensDev, spectrum.getSignal(i).getIntensity() + intensDev, spectrum.getSignal(i).getNuclei()[dim]));
            } else {
                hits.get(i).addAll(casekit.NMR.DB.getSignalIDsFromNMRShiftDB(DBConnection, shift - shiftDev, shift + shiftDev, spectrum.getSignal(i).getMultiplicity(), spectrum.getSignal(i).getIntensity(), spectrum.getSignal(i).getIntensity(), spectrum.getSignal(i).getNuclei()[dim]));
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
     * @deprecated 
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
                + "		spectype.NAME = \"" + casekit.NMR.Utils.getIsotopeIdentifier(elem) + "\";";
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
     * @deprecated 
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
                + "		spectype.NAME = \"" + casekit.NMR.Utils.getIsotopeIdentifier(elem) + "\" AND \n"
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
    
    
    public static ArrayList<String> getNMRShiftDBSpectra(final IAtomContainer ac, final String elem) {

        ArrayList<String> props = (ArrayList<String>) (ArrayList<?>) (new ArrayList<>(ac.getProperties().keySet()));
        final ArrayList<String> spectra = new ArrayList<>();
        for (String prop : props) {
            if (prop.contains("Spectrum " + casekit.NMR.Utils.getIsotopeIdentifier(elem))) {
                spectra.add(ac.getProperty(prop));
            }
        }

        return spectra;
    }
    
    
    /**
     * Creates a two dimensional array of a given NMRShiftDB NMR entry 
     * with all signal shift values, intensities, multiplicities and atom indices.
     *
     * @param NMRShiftDBSpectrum
     * @return two dimensional array: 
     * 1. dimension: signal index (row); 
     * 2. dimension: signal shift value (column 1), signal intensity (column 2), 
     * signal multiplicity (column 3), atom index in structure (column 4)
     */
    public static String[][] parseNMRShiftDBSpectrum(final String NMRShiftDBSpectrum){
        if(NMRShiftDBSpectrum.trim().isEmpty()){
            return new String[][]{};
        }
        String[] signalSplit;
        final String[] shiftsSplit = NMRShiftDBSpectrum.split("\\|");
        final String[][] values = new String[shiftsSplit.length][4];
        for (int i = 0; i < shiftsSplit.length; i++) {
            signalSplit = shiftsSplit[i].split(";");
            values[i][0] = signalSplit[0]; // shift value
            values[i][1] = signalSplit[1].substring(0, signalSplit[1].length() - 1); // intensity
            values[i][2] = signalSplit[1].substring(signalSplit[1].length() - 1); // multiplicity
            values[i][3] = signalSplit[2]; // atom index
        }
        
        return values;
    }
    
    /**
     * Sets shifts, intensities and implicit hydrogen counts in atoms of an atom container 
     * by means of given spectrum property string.
     *
     * @param ac IAtomContainer to set
     * @param NMRShiftDBSpectrum Property string of spectrum in NMRShiftDB format.
     * @return 
     * 
     * @see DB#parseNMRShiftDBSpectrum(java.lang.String) 
     * @see Utils#getHydrogenCountFromMultiplicity(java.lang.String) 
     * @deprecated 
     */
    public static boolean setNMRShiftDBShiftsToAtomContainer(final IAtomContainer ac, final String NMRShiftDBSpectrum){
        if (ac.getProperty(NMRShiftDBSpectrum) == null) {
            return false;
        }
        final String[][] spectrumStringArray = DB.parseNMRShiftDBSpectrum(ac.getProperty(NMRShiftDBSpectrum));
        
        Integer atomIndexSpectrum;
//        String multiplicity;
        Double shift;
        
        for (int i = 0; i < spectrumStringArray.length; i++) {
            atomIndexSpectrum = Integer.parseInt(spectrumStringArray[i][3]);
            shift = Double.parseDouble(spectrumStringArray[i][0]);
//            multiplicity = spectrumStringArray[i][3];
            if(Utils.checkIndexInAtomContainer(ac, atomIndexSpectrum)){
                ac.getAtom(atomIndexSpectrum).setProperty(Utils.getNMRShiftConstant(ac.getAtom(atomIndexSpectrum).getSymbol()), shift);
//                ac.getAtom(atomIndexSpectrum).setImplicitHydrogenCount(Utils.getHydrogenCountFromMultiplicity(multiplicity));
            }    
        }            
        
        return true;
    }
    
    
    public static Spectrum NMRShiftDBSpectrumToSpectrum(final String NMRShiftDBSpectrum, final String atomType){
        if ((NMRShiftDBSpectrum == null) || NMRShiftDBSpectrum.trim().isEmpty()) {
            return null;
        }
        final String[][] spectrumStringArray = DB.parseNMRShiftDBSpectrum(NMRShiftDBSpectrum);
        final Spectrum spectrum = new Spectrum(new String[]{Utils.getIsotopeIdentifier(atomType)});
        String multiplicity;
        Double shift, intensity;
        for (int i = 0; i < spectrumStringArray.length; i++) {
            shift = Double.parseDouble(spectrumStringArray[i][0]);
            intensity = Double.parseDouble(spectrumStringArray[i][1]);
            multiplicity = spectrumStringArray[i][2];            
            spectrum.addSignal(new Signal(new String[]{Utils.getIsotopeIdentifier(atomType)}, new Double[]{shift}, intensity, multiplicity));
        }
        Utils.setSpectrumEquivalences(spectrum);
        
        return spectrum;
    }
    
    public static Assignment NMRShiftDBSpectrumToAssignment(final String NMRShiftDBSpectrum, final String atomType) {
        if ((NMRShiftDBSpectrum == null) || NMRShiftDBSpectrum.trim().isEmpty()) {
            return null;
        }
        final String[][] NMRShiftDBSpectrumStringArray = DB.parseNMRShiftDBSpectrum(NMRShiftDBSpectrum);
        final Spectrum spectrum = DB.NMRShiftDBSpectrumToSpectrum(NMRShiftDBSpectrum, atomType);
        final Assignment assignment = new Assignment(spectrum);
        for (int i = 0; i < NMRShiftDBSpectrumStringArray.length; i++) {
            assignment.setAssignment(0, i, new Integer(NMRShiftDBSpectrumStringArray[i][3]));
        }        

        return assignment;
    }
}
