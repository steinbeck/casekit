package casekit.nmr.analysis;

import casekit.nmr.hose.HOSECodeBuilder;
import casekit.nmr.model.DataSet;
import casekit.nmr.model.Signal;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.reflect.TypeToken;
import org.bson.Document;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.io.*;
import java.util.*;

public class HOSECodeShiftStatistics {

    private final static Gson GSON = new GsonBuilder().setLenient()
                                                      .create(); //.setPrettyPrinting()

    public static Map<String, Map<String, List<Double>>> collectHOSECodeShifts(final List<DataSet> dataSetList,
                                                                               final int maxSphere,
                                                                               final boolean use3D) {
        return collectHOSECodeShifts(dataSetList, maxSphere, use3D, new HashMap<>());
    }

    public static Map<String, Map<String, List<Double>>> collectHOSECodeShifts(final List<DataSet> dataSetList,
                                                                               final int maxSphere, final boolean use3D,
                                                                               final Map<String, Map<String, List<Double>>> hoseCodeShifts) {
        //        final ExtendedHOSECodeGenerator extendedHOSECodeGenerator = new ExtendedHOSECodeGenerator();
        IAtomContainer structure;
        Signal signal;
        String hoseCode;
        String solvent;

        for (final DataSet dataSet : dataSetList) {
            structure = dataSet.getStructure()
                               .toAtomContainer();
            //            if (use3D) {
            //                try {
            //                    /* !!! No explicit H in mol !!! */
            //                    Utils.convertExplicitToImplicitHydrogens(structure);
            //                    /* add explicit H atoms */
            //                    AtomUtils.addAndPlaceHydrogens(structure);
            //                    /* detect aromaticity */
            //                    Utils.setAromaticityAndKekulize(structure);
            //                } catch (final IOException | ClassNotFoundException | CDKException e) {
            //                    e.printStackTrace();
            //                    continue;
            //                }
            //            }
            solvent = dataSet.getSpectrum()
                             .getSolvent();
            if (solvent
                    == null
                    || solvent.equals("")) {
                solvent = "Unknown";
            }
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
                            //                            if (use3D) {
                            //                                hoseCode = extendedHOSECodeGenerator.getHOSECode(structure, structure.getAtom(i),
                            //                                                                                 maxSphere);
                            //                            } else {
                            hoseCode = HOSECodeBuilder.buildHOSECode(structure, i, sphere, false);
                            //                            }
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

        return hoseCodeShifts;
    }

    public static Map<String, Map<String, Double[]>> buildHOSECodeShiftStatistics(
            final Map<String, Map<String, List<Double>>> hoseCodeShifts) {

        final Map<String, Map<String, Double[]>> hoseCodeShiftStatistics = new HashMap<>();
        List<Double> values;
        for (final Map.Entry<String, Map<String, List<Double>>> hoseCodes : hoseCodeShifts.entrySet()) {
            hoseCodeShiftStatistics.put(hoseCodes.getKey(), new HashMap<>());
            for (final Map.Entry<String, List<Double>> solvents : hoseCodes.getValue()
                                                                           .entrySet()) {
                values = solvents.getValue(); //casekit.nmr.Utils.removeOutliers(solvents.getValue(), 1.5);
                hoseCodeShiftStatistics.get(hoseCodes.getKey())
                                       .put(solvents.getKey(),
                                            new Double[]{(double) values.size(), Collections.min(values),
                                                         casekit.nmr.Utils.getMean(values),
                                                         // casekit.nmr.Utils.getRMS(values),
                                                         casekit.nmr.Utils.getMedian(values), Collections.max(values)});
            }
        }

        return hoseCodeShiftStatistics;
    }

    public static boolean writeHOSECodeShiftStatistics(final Map<String, Map<String, Double[]>> hoseCodeShifts,
                                                       final String pathToJsonFile) {
        try {
            final BufferedWriter bw = new BufferedWriter(new FileWriter(pathToJsonFile));
            bw.append("{");
            bw.newLine();
            bw.flush();

            Document subDocument;
            String json;
            long counter = 0;
            for (final Map.Entry<String, Map<String, Double[]>> entry : hoseCodeShifts.entrySet()) {
                subDocument = new Document();
                subDocument.append("HOSECode", entry.getKey());
                subDocument.append("values", GSON.toJson(entry.getValue()));
                json = new Document(String.valueOf(counter), subDocument).toJson();
                bw.append(json, 1, json.length()
                        - 1);
                if (counter
                        < hoseCodeShifts.size()
                        - 1) {
                    bw.append(",");
                }
                bw.newLine();
                bw.flush();

                counter++;
            }

            bw.append("}");
            bw.flush();
            bw.close();

            return true;
        } catch (final IOException e) {
            e.printStackTrace();
        }

        return false;
    }

    public static Map<String, Map<String, Double[]>> readHOSECodeShiftStatistics(
            final String pathToJsonFile) throws FileNotFoundException {
        final BufferedReader br = new BufferedReader(new FileReader(pathToJsonFile));
        final Map<String, Map<String, Double[]>> hoseCodeShiftStatistics = new HashMap<>();
        // add all task to do
        br.lines()
          .forEach(line -> {
              if ((line.trim()
                       .length()
                      > 1)
                      || (!line.trim()
                               .startsWith("{")
                      && !line.trim()
                              .endsWith("}"))) {
                  final StringBuilder hoseCodeShiftsStatisticInJSON = new StringBuilder();
                  if (line.endsWith(",")) {
                      hoseCodeShiftsStatisticInJSON.append(line, 0, line.length()
                              - 1);
                  } else {
                      hoseCodeShiftsStatisticInJSON.append(line);
                  }
                  final JsonObject jsonObject = JsonParser.parseString(hoseCodeShiftsStatisticInJSON.substring(
                          hoseCodeShiftsStatisticInJSON.toString()
                                                       .indexOf("{")))
                                                          .getAsJsonObject();
                  hoseCodeShiftStatistics.put(jsonObject.get("HOSECode")
                                                        .getAsString(), GSON.fromJson(jsonObject.get("values")
                                                                                                .getAsString(),
                                                                                      new TypeToken<Map<String, Double[]>>() {
                                                                                      }.getType()));
              }
          });

        return hoseCodeShiftStatistics;
    }
}
