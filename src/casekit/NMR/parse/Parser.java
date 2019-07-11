/*
 * The MIT License
 *
 * Copyright 2019 Michael Wenk [https://github.com/michaelwenk]
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package casekit.NMR.parse;

import casekit.NMR.Utils;
import casekit.NMR.model.Signal;
import casekit.NMR.model.Spectrum;
import casekit.io.FileParser;
import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;
import javax.xml.parsers.ParserConfigurationException;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;


public class Parser {


    /**
     * Creates a Spectrum class object from given 1D NMR input file in CSV or XML format.
     * The extension of given file is used to determine the format.
     *
     * @param pathToFile path to peak table (Bruker's TopSpin csv or xml
     * file format)
     * @param nucleus nucleus to use for spectrum creation, e.g. "13C"
     * @return
     * @throws Exception
     */
    public static Spectrum parse1DNMR(final String pathToFile, final String nucleus) throws Exception {
        switch (Utils.getFileFormat(pathToFile)){
            case "csv":
                return CSVtoSpectrum(pathToFile, new int[]{4}, new String[]{nucleus}, 6);
            case "xml":
                return XMLtoSpectrum(pathToFile, 1, new int[]{1}, new String[]{nucleus});
            default:
                return null;
        }
    }

    /**
     * Creates a Spectrum class object from given 2D NMR input file in CSV or XML format.
     * The extension of given file is used to determine the format.
     *
     * @param pathToFile path to peak table (Bruker's TopSpin csv or xml
     * file format)
     * @param nuclei nuclei to use for spectrum creation, e.g. ["13C", "13C]
     *
     * @return
     * @throws Exception
     */
    public static Spectrum parse2DNMR(final String pathToFile, final String[] nuclei) throws Exception {
        switch (Utils.getFileFormat(pathToFile)){
            case "csv":
                return CSVtoSpectrum(pathToFile, new int[]{5, 6}, nuclei, 9);
            case "xml":
                return XMLtoSpectrum(pathToFile, 2, new int[]{2, 1}, nuclei);
            default:
                return null;
        }
    }

    /**
     * Reads a specific column of a NMR peak table and stores it into an
     * ArrayList object.
     *
     * @param pathToCSV path to NMR peak table in CSV file format
     * @param column column index to select in peak table
     * @return ArrayList of Double shift values
     * @throws IOException
     */
    private static ArrayList<Double> CSVtoPeakList(final String pathToCSV, final int column) throws IOException {
        final ArrayList<Double> shifts = new ArrayList<>();
        String line;
        String[] tokens;
        final BufferedReader fileReader = FileParser.parseText(pathToCSV);
        while ((line = fileReader.readLine()) != null) {
            tokens = line.split(",");
            // get shift value
            if (tokens[column].trim().matches("^[+|-]{0,1}\\d+\\.{0,1}\\d*")) {
                shifts.add(Double.parseDouble(tokens[column].trim()));
            }
        }
        fileReader.close();

        return shifts;
    }

    /**
     * Reads specific columns of one NMR peak table to obtain a Spectrum class
     * object and set intensitiy values.
     * The number of columns and atom types has to be the same and defines the
     * dimension of the returning spectrum.
     *
     * @param pathToCSV path to NMR peak table in CSV file format
     * @param columns column indices to select in peak table
     * @param nuclei nuclei for each dimension
     * @param intensityColumnIndex column index for intensity values
     * @return Spectrum class object containing the peak lists
     * @throws Exception
     */
    private static Spectrum CSVtoSpectrum(final String pathToCSV, final int[] columns, final String[] nuclei, final int intensityColumnIndex) throws Exception {
        // assumes the same number of selected columns (dimensions) and atom types
        if(columns.length != nuclei.length){
            return null;
        }
        final Spectrum spectrum = new Spectrum(nuclei);
        ArrayList<Double> shiftList;
        for (int col = 0; col < columns.length; col++) {
            shiftList = CSVtoPeakList(pathToCSV, columns[col]);
            if(col == 0){
                for (int i = 0; i < shiftList.size(); i++) {
                    spectrum.addSignal(new Signal(spectrum.getNuclei()));
                }
            }
            if(!spectrum.setShifts(shiftList, col)){
                return null;
            }
        }
        spectrum.setIntensities(CSVtoPeakList(pathToCSV, intensityColumnIndex));

        return spectrum;
     }

    /**
     * Reads a NMR peak XML file and returns one attribute of nodes (column) into an
     * ArrayList object.
     * The XML file must be in Bruker's TopSpin format.
     *
     * @param pathToXML Path to XML file
     * @param dim number of dimensions of given data 1 (1D) or 2 (2D)
     * @param attribute which attribute index in XML peak nodes should be used:
     * 1 (shift of 1st dimension), 2 (shift of 2nd dimension if 2D data,
     * intensity if 1D data) or 3 (intensity if 2D data)
     *
     * @return ArrayList of Double shift values
     * @throws IOException
     * @throws javax.xml.parsers.ParserConfigurationException
     * @throws org.xml.sax.SAXException
     */
    private static ArrayList<Double> XMLtoPeakList(final String pathToXML, final int dim, final int attribute) throws IOException, ParserConfigurationException, SAXException {
        // assumes a attribute value between 1 and 3
        if(attribute < 1 || attribute > 3){
            return null;
        }

        final ArrayList<Double> shifts = new ArrayList<>();
        final Document doc = FileParser.parseXML(pathToXML);

        final NodeList peakLists = doc.getElementsByTagName("Peak" + dim + "D");
        for (int i = 0; i < peakLists.getLength(); i++) {
            shifts.add(Double.parseDouble(peakLists.item(i).getAttributes().item(attribute - 1).getNodeValue()));
        }

        return shifts;
    }

    /**
     * Reads specific columns of NMR XML files to obtain a Spectrum class
     * object.
     * The XML file must be in Bruker's TopSpin format.
     *
     * @param pathToXML path to NMR XML file in Bruker's TopSpin XML file format
     * @param ndim number of dimensions: 1 (1D) or 2 (2D)
     * @param attributes which attribute indices in XML peak nodes should be used:
     * 1 (shift of 1st dimension), 2 (shift of 2nd dimension if 2D data)
     * @param nuclei nuclei for each dimension
     * @return Spectrum class object containing the selected peak lists
     * @throws Exception
     */
    private static Spectrum XMLtoSpectrum(final String pathToXML, final int ndim, final int[] attributes, final String[] nuclei) throws Exception {

        // assumes the same number of dims, attributes and atom types and a maximum number of dims of 2
        if((ndim != attributes.length) || (ndim != nuclei.length) || (attributes.length != nuclei.length)
                || (ndim < 1 || ndim > 2)){
            return null;
        }
        final Spectrum spectrum = new Spectrum(nuclei);
        ArrayList<Double> shiftList;
        for (int dim = 0; dim < ndim; dim++) {
            shiftList = XMLtoPeakList(pathToXML, ndim, attributes[dim]);
            if(dim == 0){
                for (int i = 0; i < (shiftList != null ? shiftList.size() : 0); i++) {
                    spectrum.addSignal(new Signal(spectrum.getNuclei()));
                }
            }
            if(!spectrum.setShifts(shiftList, dim)){
                return null;
            }
        }
        spectrum.setIntensities(XMLtoPeakList(pathToXML, ndim, ndim + 1));

        return spectrum;
     }
}
