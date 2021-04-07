package casekit.nmr.lsd.model;

public class ElucidationOptions {

    // PyLSD options
    private String[] filterPaths;
    private boolean allowHeteroHeteroBonds;
    private boolean useElim;
    private int elimP1;
    private int elimP2;

    public ElucidationOptions() {
    }

    public String[] getFilterPaths() {
        return this.filterPaths;
    }

    public void setFilterPaths(final String[] filterPaths) {
        this.filterPaths = filterPaths;
    }

    public boolean isAllowHeteroHeteroBonds() {
        return this.allowHeteroHeteroBonds;
    }

    public void setAllowHeteroHeteroBonds(final boolean allowHeteroHeteroBonds) {
        this.allowHeteroHeteroBonds = allowHeteroHeteroBonds;
    }

    public boolean isUseElim() {
        return this.useElim;
    }

    public void setUseElim(final boolean useElim) {
        this.useElim = useElim;
    }

    public int getElimP1() {
        return this.elimP1;
    }

    public void setElimP1(final int elimP1) {
        this.elimP1 = elimP1;
    }

    public int getElimP2() {
        return this.elimP2;
    }

    public void setElimP2(final int elimP2) {
        this.elimP2 = elimP2;
    }
}
