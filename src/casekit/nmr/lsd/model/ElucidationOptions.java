package casekit.nmr.lsd.model;

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

    public int getHmbcP3() {
        return this.hmbcP3;
    }

    public void setHmbcP3(final int hmbcP3) {
        this.hmbcP3 = hmbcP3;
    }

    public int getHmbcP4() {
        return this.hmbcP4;
    }

    public void setHmbcP4(final int hmbcP4) {
        this.hmbcP4 = hmbcP4;
    }

    public int getCosyP3() {
        return this.cosyP3;
    }

    public void setCosyP3(final int cosyP3) {
        this.cosyP3 = cosyP3;
    }

    public int getCosyP4() {
        return this.cosyP4;
    }

    public void setCosyP4(final int cosyP4) {
        this.cosyP4 = cosyP4;
    }
}
