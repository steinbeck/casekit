package casekit.nmr.lsd.model;

import lombok.AllArgsConstructor;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

@NoArgsConstructor
@AllArgsConstructor
@Getter
@Setter
public class ElucidationOptions {

    // PyLSD options
    private String[] filterPaths;
    private boolean allowHeteroHeteroBonds;
    private boolean useElim;
    private int elimP1;
    private int elimP2;
    private int hmbcP3;
    private int hmbcP4;
    private int cosyP3;
    private int cosyP4;
}
