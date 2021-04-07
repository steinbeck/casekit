package casekit.nmr.lsd.model;

public class ElucidationOptions {

    // PyLSD options
    private String pathToLSDFilterList;
    private boolean allowHeteroHeteroBonds;
    private int elimP1;
    private int elimP2;

    public ElucidationOptions() {
    }

    public String getPathToLSDFilterList() {
        return this.pathToLSDFilterList;
    }

    public void setPathToLSDFilterList(final String pathToLSDFilterList) {
        this.pathToLSDFilterList = pathToLSDFilterList;
    }

    public boolean isAllowHeteroHeteroBonds() {
        return this.allowHeteroHeteroBonds;
    }

    public void setAllowHeteroHeteroBonds(final boolean allowHeteroHeteroBonds) {
        this.allowHeteroHeteroBonds = allowHeteroHeteroBonds;
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
