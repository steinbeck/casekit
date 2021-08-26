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

    private StructureCompact structure;
    private SpectrumCompact spectrum;
    private Assignment assignment;
    private Map<String, String> meta;

    public DataSet(final IAtomContainer structure, final Spectrum spectrum, final Assignment assignment,
                   final Map<String, String> meta) {
        this.structure = new StructureCompact(structure);
        this.spectrum = new SpectrumCompact(spectrum);
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

    public DataSet buildClone() {
        final Map<String, String> metaTemp = this.meta
                                                     == null
                                             ? new HashMap<>()
                                             : new HashMap<>(this.meta);
        return new DataSet(this.structure.buildClone(), this.spectrum.buildClone(), this.assignment.buildClone(),
                           new HashMap<>(metaTemp));
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
