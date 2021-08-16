package casekit.nmr.model;

import lombok.AllArgsConstructor;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.HashMap;
import java.util.Map;

@NoArgsConstructor
@AllArgsConstructor
@Getter
@Setter
public class DataSet {

    private ExtendedConnectionMatrix structure;
    private Spectrum spectrum;
    private Assignment assignment;
    private Map<String, String> meta;

    public DataSet(final IAtomContainer structure, final Spectrum spectrum, final Assignment assignment,
                   final Map<String, String> meta) {
        this.structure = new ExtendedConnectionMatrix(structure);
        this.spectrum = spectrum;
        this.assignment = assignment;
        this.meta = new HashMap<>(meta);
    }

    public void addMetaInfo(final String key, final String value) {
        if (this.meta
                == null) {
            this.meta = new HashMap<>();
        }
        this.meta.put(key, value);
    }

    public void removeMetaInfo(final String key) {
        this.meta.remove(key);
    }

    @Override
    public String toString() {
        return "DataSet{"
                + "structure="
                + this.structure
                + ", spectrum="
                + this.spectrum
                + ", assignment="
                + this.assignment
                + ", meta="
                + this.meta
                + '}';
    }
}
