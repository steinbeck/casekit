package casekit.nmr.model;

import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.HashMap;
import java.util.Map;

public class DataSet {

    private ExtendedConnectionMatrix structure;
    private Spectrum spectrum;
    private Assignment assignment;
    private HashMap<String, String> meta;

    public DataSet() {
    }

    public DataSet(final ExtendedConnectionMatrix structure, final Spectrum spectrum, final Assignment assignment, Map<String, String> meta) {
        this.structure = structure;
        this.spectrum = spectrum;
        this.assignment = assignment;
        this.meta = new HashMap<>(meta);
    }

    public DataSet(final IAtomContainer structure, final Spectrum spectrum, final Assignment assignment, Map<String, String> meta) {
        this.structure = new ExtendedConnectionMatrix(structure);
        this.spectrum = spectrum;
        this.assignment = assignment;
        this.meta = new HashMap<>(meta);
    }

    public void addMetaInfo(final String key, final String value){
        this.meta.put(key, value);
    }

    public void removeMetaInfo(final String key){
        this.meta.remove(key);
    }

    public ExtendedConnectionMatrix getStructure() {
        return structure;
    }

    public void setStructure(final ExtendedConnectionMatrix structure) {
        this.structure = structure;
    }

    public Spectrum getSpectrum() {
        return spectrum;
    }

    public void setSpectrum(final Spectrum spectrum) {
        this.spectrum = spectrum;
    }

    public Assignment getAssignment() {
        return assignment;
    }

    public void setAssignment(final Assignment assignment) {
        this.assignment = assignment;
    }

    public Map<String, String> getMeta() {
        return meta;
    }

    public void setMeta(Map<String, String> meta) {
        this.meta = new HashMap<>(meta);
    }

    @Override
    public String toString() {
        return "DataSet{" +
                "structure=" + structure +
                ", spectrum=" + spectrum +
                ", assignment=" + assignment +
                ", meta=" + meta +
                '}';
    }
}
