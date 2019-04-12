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
import com.mongodb.MongoClient;
import com.mongodb.MongoClientOptions;
import com.mongodb.MongoCredential;
import com.mongodb.ServerAddress;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoDatabase;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import org.bson.Document;
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
     * Returns all spectra for each molecule and a given nucleus which exist as 
     * property in a NMRSHiftDB SDF.
     *
     * @param pathToNMRShiftDB path to NMRShiftDB file
     * @param nucleus nucleus of requested spectra
     * @return
     * @throws FileNotFoundException
     * @throws CDKException
     * 
     */
    public static ArrayList<ArrayList<Object[]>> getSpectraFromNMRShiftDB(final String pathToNMRShiftDB, final String nucleus) throws FileNotFoundException, CDKException {
        final ArrayList<ArrayList<Object[]>> spectraSet = new ArrayList<>();
        final IteratingSDFReader iterator = new IteratingSDFReader(
                new FileReader(pathToNMRShiftDB),
                SilentChemObjectBuilder.getInstance()
        );
        IAtomContainer ac;
        Spectrum spectrum;
        ArrayList<Object[]> spectra;
        HashMap<String, String> spectraStrings;
        String spectrumIndexInRecord, solvent;
        while (iterator.hasNext()) {
            ac = iterator.next();
            if(ac == null){     
                continue;
            }
            spectraStrings = DB.getSpectraStrings(ac, nucleus);
            if(spectraStrings.isEmpty() || (ac.getProperty("Solvent") == null)){
                continue;
            }
            spectra = new ArrayList<>();
            for (final String spectrumPropertyString : spectraStrings.keySet()) {  
                spectrum = DB.NMRShiftDBSpectrumToSpectrum(spectraStrings.get(spectrumPropertyString), nucleus);
                if(spectrum == null){
                    continue;
                }
                spectrumIndexInRecord = spectrumPropertyString.split("\\s")[spectrumPropertyString.split("\\s").length - 1]; 
                solvent = DB.getSolvent(ac.getProperty("Solvent"), spectrumIndexInRecord);
                if(solvent == null){
                    continue;
                }
                spectrum.setSolvent(solvent);
                
                if(Utils.getAtomTypeIndicesByElement(ac, nucleus.replaceAll("\\d", "")).size() != spectrum.getSignalCount()){
                    continue;
                }
                
                spectra.add(new Object[]{spectrum, DB.NMRShiftDBSpectrumToAssignment(spectraStrings.get(spectrumPropertyString), nucleus)});
            }
            spectraSet.add(spectra);
        }

        return spectraSet;
    }
    
    public static String getSolvent(final String solventPropertyString, final String spectrumIndexInRecord){
        final String[] solventPropertyStringSplit = solventPropertyString.split(":");
        String solvent;
        for (int i = 0; i < solventPropertyStringSplit.length; i++) {
            if (solventPropertyStringSplit[i].endsWith(spectrumIndexInRecord)) {
                solvent = solventPropertyStringSplit[i + 1];
                if(solvent.substring(solvent.length() - 1).matches("\\d")){
                    solvent = solvent.substring(0, solvent.length() - 1);                    
                }
                if(solvent.substring(solvent.length() - 1).matches("\\d")){
                    solvent = solvent.substring(0, solvent.length() - 1);                    
                }
                solvent = solvent.substring(0, solvent.length() - 1);
                
                return solvent;
            }
        }
        
        return null;
    }
    
    /**
     * Returns 3-tuples consisting of structure, spectrum and assignments 
     * for each valid molecule record in the given NMRShiftDB file. Valid means 
     * here that each molecule record has to contain the given spectrum 
     * property string as well as the number of signals in that spectrum has to 
     * be the same as atoms of that atom type in molecule.
     *
     * @param pathToNMRShiftDB path to NMRShiftDB file
     * @param NMRShiftDBSpectrumProperty spectrum property string to use
     * @return
     * @throws FileNotFoundException
     * @throws CDKException
     */
    public static HashMap<Integer, Object[]> getSSCComponentsFromNMRShiftDB(final String pathToNMRShiftDB, final String NMRShiftDBSpectrumProperty) throws FileNotFoundException, CDKException {
        final HashMap<Integer, Object[]> structureSetWithSpectra = new HashMap<>();
        final IteratingSDFReader iterator = new IteratingSDFReader(
                new FileReader(pathToNMRShiftDB),
                SilentChemObjectBuilder.getInstance()
        );
        IAtomContainer ac;
        Spectrum spectrum;
        Assignment assignment;
        final String nucleus = DB.getNucleusFromNMRShiftDBSpectrumProperty(NMRShiftDBSpectrumProperty);
        final String spectrumIndexInRecord = NMRShiftDBSpectrumProperty.split("\\s")[NMRShiftDBSpectrumProperty.split("\\s").length - 1];        
        while (iterator.hasNext()) {
            ac = iterator.next();
            // skip molecules which not contain any of requested spectrum information
            if(ac.getProperty(NMRShiftDBSpectrumProperty) == null){
                continue;
            }
            spectrum = DB.NMRShiftDBSpectrumToSpectrum(ac.getProperty(NMRShiftDBSpectrumProperty), nucleus);           
            // if no spectrum could be built or the number of signals in spectrum is different than the atom number in molecule
            if((spectrum == null) || Utils.getAtomTypeIndicesByElement(ac, nucleus.replaceAll("\\d", "")).size() != spectrum.getSignalCount()){
                continue;
            }                       
            if(ac.getProperty("Solvent") != null){
                spectrum.setSolvent(DB.getSolvent(ac.getProperty("Solvent"), spectrumIndexInRecord));
            }
            if(ac.getProperty("Field Strength [MHz]") != null){
                for (final String fieldStrength : ac.getProperty("Field Strength [MHz]").toString().split("\\s")) {
                    if (fieldStrength.startsWith(spectrumIndexInRecord + ":")) {
                        try {
                            spectrum.setSpectrometerFrequency(Double.parseDouble(fieldStrength.split(spectrumIndexInRecord + ":")[1]));
                        } catch (NumberFormatException e) {
//                            spectrum.setSpectrometerFrequency(null);
                        }
                        break;
                    }
                }
            }            
            
            assignment = DB.NMRShiftDBSpectrumToAssignment(ac.getProperty(NMRShiftDBSpectrumProperty), nucleus);
//            if ((ac != null) && (spectrum != null)) {                
                structureSetWithSpectra.put(structureSetWithSpectra.size(), new Object[]{ac, spectrum, assignment});
//            }
            
            Utils.setAromaticitiesInAtomContainer(ac);
        }

        return structureSetWithSpectra;
    }

    /**
     * Returns a hashmap containing combined keys (by "_") of solvents 
     * and lists of calculated deviations between all given spectra for a 
     * nucleus in molecule record as values. <br>
     * Here, only molecule records in NMRShiftDB file are considered which have
     * at least two different spectra for same nucleus. <br>
     * Example: "Spectrum 13C 0", "Spectrum 13C 1" will be used for given 
     * nucleus 13C.
     * 
     *
     * @param pathToNMRShiftDB
     * @param nucleus
     * @return
     * @throws FileNotFoundException
     * @throws CDKException
     */
    public static HashMap<String, ArrayList<Double>> getSolventDeviations(final String pathToNMRShiftDB, final String nucleus) throws FileNotFoundException, CDKException{
        int signalCount;
        Spectrum spectrum;
        Assignment assignment;
        final ArrayList<ArrayList<Object[]>> spectraSets = DB.getSpectraFromNMRShiftDB(pathToNMRShiftDB, nucleus);
        HashMap<Integer, ArrayList<Double>> shiftsPerAtom;
        HashMap<Integer, ArrayList<String>> solventsPerAtom;
        ArrayList<String> solvents;
        String[] solventsToSort;

        final HashMap<String, ArrayList<Double>> deviations = new HashMap<>();
        String combiKey;

        for (final ArrayList<Object[]> spectraSetInRecord : spectraSets) {
            shiftsPerAtom = new HashMap<>();
            solventsPerAtom = new HashMap<>();
            signalCount = -1;
            for (final Object[] spectrumAndAssignment : spectraSetInRecord) {
                spectrum = (Spectrum) spectrumAndAssignment[0];
                assignment = (Assignment) spectrumAndAssignment[1];
                if (signalCount == -1) {
                    signalCount = spectrum.getSignalCount();
                } else if (signalCount != spectrum.getSignalCount()) {
                    continue;
                }
                for (final int atomIndex : assignment.getAtomIndices(0)) {
                    if (!shiftsPerAtom.containsKey(atomIndex)) {
                        shiftsPerAtom.put(atomIndex, new ArrayList<>());
                        solventsPerAtom.put(atomIndex, new ArrayList<>());
                    }
                    shiftsPerAtom.get(atomIndex).add(spectrum.getSignal(assignment.getSignalIndex(0, atomIndex)).getShift(0));
                    solventsPerAtom.get(atomIndex).add(spectrum.getSolvent());
                }
            }            
            if (shiftsPerAtom.isEmpty() || (shiftsPerAtom.get(Collections.min(shiftsPerAtom.keySet())).size() < 2)) {
                continue;
            }
            solvents = new ArrayList<>(solventsPerAtom.get(Collections.min(solventsPerAtom.keySet())));
//                if(Collections.frequency(solvents, "Unreported") + Collections.frequency(solvents, "Unknown") > solvents.size() - 2){
//                    continue;
//                }

            for (final int atomIndex : shiftsPerAtom.keySet()) {
                for (int s1 = 0; s1 < solvents.size(); s1++) {
//                        if(solvents.get(s1).equals("Unreported") || solvents.get(s1).equals("Unknown")){
//                            continue;
//                        }
                    for (int s2 = s1 + 1; s2 < solvents.size(); s2++) {
//                            if (solvents.get(s2).equals("Unreported") || solvents.get(s2).equals("Unknown")) {
//                                continue;
//                            }
                        solventsToSort = new String[2];
                        solventsToSort[0] = solvents.get(s1);
                        solventsToSort[1] = solvents.get(s2);
                        Arrays.sort(solventsToSort);
                        combiKey = solventsToSort[0] + "_" + solventsToSort[1];
                        if (!deviations.containsKey(combiKey)) {
                            deviations.put(combiKey, new ArrayList<>());
                        }
                        deviations.get(combiKey).add(Math.abs(shiftsPerAtom.get(atomIndex).get(s1) - shiftsPerAtom.get(atomIndex).get(s2)));
                    }
                }
            }
        }

        return deviations;
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
    
    
    public static HashMap<String, String> getSpectraStrings(final IAtomContainer ac, final String nucleus) {
        final ArrayList<String> props = (ArrayList<String>) (ArrayList<?>) (new ArrayList<>(ac.getProperties().keySet()));
        final HashMap<String, String> spectra = new HashMap<>();
        for (final String prop : props) {
            if (prop.startsWith("Spectrum " + nucleus)) {
                spectra.put(prop, ac.getProperty(prop));
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
    
    public static String getNucleusFromNMRShiftDBSpectrumProperty(final String NMRShiftDBSpectrumProperty){        
        return NMRShiftDBSpectrumProperty.split(" ")[1];
    }
    
    public static Spectrum NMRShiftDBSpectrumToSpectrum(final String NMRShiftDBSpectrum, final String nucleus){
        if ((NMRShiftDBSpectrum == null) || NMRShiftDBSpectrum.trim().isEmpty()) {
            return null;
        }
        final String[][] spectrumStringArray = DB.parseNMRShiftDBSpectrum(NMRShiftDBSpectrum);
        final Spectrum spectrum = new Spectrum(new String[]{nucleus});
        String multiplicity;
        Double shift, intensity;
        try {
            for (int i = 0; i < spectrumStringArray.length; i++) {
                shift = Double.parseDouble(spectrumStringArray[i][0]);
                intensity = Double.parseDouble(spectrumStringArray[i][1]);
                multiplicity = spectrumStringArray[i][2];
                spectrum.addSignal(new Signal(new String[]{nucleus}, new Double[]{shift}, multiplicity, intensity));
            }
            Utils.setSpectrumEquivalences(spectrum);
        } catch (Exception e) {
            
            return null;
        }
        
        return spectrum;
    }
    
    public static Assignment NMRShiftDBSpectrumToAssignment(final String NMRShiftDBSpectrum, final String nucleus) {
        if ((NMRShiftDBSpectrum == null) || NMRShiftDBSpectrum.trim().isEmpty()) {
            return null;
        }
        final String[][] NMRShiftDBSpectrumStringArray = DB.parseNMRShiftDBSpectrum(NMRShiftDBSpectrum);
        final Spectrum spectrum = DB.NMRShiftDBSpectrumToSpectrum(NMRShiftDBSpectrum, nucleus);
        final Assignment assignment = new Assignment(spectrum);
        for (int i = 0; i < NMRShiftDBSpectrumStringArray.length; i++) {
            assignment.setAssignment(0, i, new Integer(NMRShiftDBSpectrumStringArray[i][3]));
        }        

        return assignment;
    }    
    
    public static MongoClient login(final String mongoUser, final String mongoPassword, final String mongoAuthDB) throws CDKException {
        MongoClient mongo;
        try {
            // Creating a Mongo client   
            mongo = new MongoClient(
                    new ServerAddress("127.0.0.1", 27017),
                    MongoCredential.createCredential(
                            mongoUser,
                            mongoAuthDB,
                            mongoPassword.toCharArray()),
                    MongoClientOptions.builder().build());
            System.out.println("Login to MongoDB was successfull");
            // Accessing the database             
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println(Thread.currentThread().getStackTrace()[1].getMethodName() + ": could not connect to MongoDB!");

            return null;
        }

        return mongo;
    }

    public static MongoDatabase getDatabase(final MongoClient mongo, final String mongoDBName){
        return mongo.getDatabase(mongoDBName);
    }
    
    public static MongoCollection<Document> getCollection(final MongoClient mongo, final String mongoDBName, final String mongoDBCollection) {
        final MongoDatabase database = DB.getDatabase(mongo, mongoDBName);
        if (database == null) {
            return null;
        }
        System.out.println("Access to database \"" + mongoDBName + "\" was successfull");
        // Retrieving a collection
        final MongoCollection<Document> collection = database.getCollection(mongoDBCollection);
        System.out.println("Retrieval of collection \"" + mongoDBCollection + "\" was successfull -> size: " + collection.countDocuments());

        return collection;
    }

    public static void logout(final MongoClient mongo) {
        mongo.close();
    }
}
