package casekit.nmr.utils;

import casekit.nmr.model.DataSet;
import casekit.nmr.model.ExtendedConnectionMatrix;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class SDFParser {

    public static List<DataSet> parseSDFile(final String pathToFile) throws CDKException, FileNotFoundException {
        return parseSDFile(new FileReader(pathToFile));
    }

    public static List<DataSet> parseSDFileContent(final String fileContent) throws CDKException {
        final InputStream inputStream = new ByteArrayInputStream(fileContent.getBytes(StandardCharsets.UTF_8));
        return parseSDFile(new InputStreamReader(inputStream));
    }

    public static List<DataSet> parseSDFile(final Reader fileReader) throws CDKException {
        final List<DataSet> dataSetList = new ArrayList<>();
        final IteratingSDFReader iterator = new IteratingSDFReader(fileReader, SilentChemObjectBuilder.getInstance());
        IAtomContainer structure;
        Map<String, String> meta;
        final CDKHydrogenAdder hydrogenAdder = CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance());
        IMolecularFormula mf;
        DataSet dataSet;

        while (iterator.hasNext()) {
            structure = iterator.next();
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(structure);
            hydrogenAdder.addImplicitHydrogens(structure);
            Utils.setAromaticityAndKekulize(structure);
            meta = new HashMap<>();
            meta.put("title", structure.getTitle());
            mf = Utils.getMolecularFormulaFromAtomContainer(structure);
            meta.put("mf", Utils.molecularFormularToString(mf));
            try {
                final String smiles = Utils.getSmilesFromAtomContainer(structure);
                meta.put("smiles", smiles);
            } catch (final CDKException e) {
                e.printStackTrace();
            }
            dataSet = new DataSet();
            dataSet.setStructure(new ExtendedConnectionMatrix(structure));
            dataSet.setMeta(meta);

            dataSetList.add(dataSet);
        }

        return dataSetList;
    }
}
