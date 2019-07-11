/*
 * The MIT License
 *
 * Copyright 2018 Michael Wenk [https://github.com/michaelwenk].
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package casekit.NMR.model;

import casekit.NMR.model.dimensional.DimensionalNMR;
import org.apache.commons.lang3.ArrayUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Assignment extends DimensionalNMR implements Cloneable {
    
    int[][] assignments;

    public Assignment(final Spectrum spectrum) {
        super(spectrum.getNuclei());
        this.assignments = this.initAssignments(this.getNDim(), spectrum.getSignalCount());
    }
    
    private int[][] initAssignments(final int nDim, final int nSignals){
        final int[][] temp = new int[nDim][nSignals];
        for (int i = 0; i < nDim; i++) {
            for (int j = 0; j < nSignals; j++) {
                temp[i][j] = -1;
            }
        }
        
        return temp;
    }
 
    /**
     * Sets an assignment as atom index for a signal position.
     *
     * @param dim
     * @param indexInSpectrum
     * @param indexInAtomContainer
     * @return
     */
    public boolean setAssignment(final int dim, final int indexInSpectrum, final int indexInAtomContainer){
        if(!this.containsDim(dim) || !this.checkSpectrumIndex(dim, indexInSpectrum)){
            return false;
        }
        this.assignments[dim][indexInSpectrum] = indexInAtomContainer;
        
        return true;
    }
    
    public boolean setAssignments(final int dim, final List<Integer> indicesInAtomContainer){
        if(!this.containsDim(dim) || !this.checkInputListSize(indicesInAtomContainer.size())){
            return false;
        }
        for (int i = 0; i < this.getAssignmentsCount(); i++) {
            this.setAssignment(dim, i, indicesInAtomContainer.get(i));
        }
        
        return true;
    }
    
    public Integer getAtomIndex(final int dim, final int indexInSpectrum){
        if(!this.containsDim(dim) || !this.checkSpectrumIndex(dim, indexInSpectrum)){
            return null;
        }

        return this.assignments[dim][indexInSpectrum];
    }
    
    public Integer getSignalIndex(final int dim, final int atomIndexInStructure){
        if(!this.containsDim(dim)){
            return null;
        }
        for (int signalIndex = 0; signalIndex < this.assignments[dim].length; signalIndex++) {
            if(this.getAtomIndex(dim, signalIndex) == atomIndexInStructure){
                return signalIndex;
            }
        }
        
        return -1;
    }
    
    public List<Integer> getAtomIndices(final int dim){
        if(!this.containsDim(dim)){
            return null;
        }

        return Arrays.asList(ArrayUtils.toObject(this.assignments[dim]));
    }
    
    public int getAssignmentsCount(){
        if(this.getNDim() > 0){
            return this.assignments[0].length;
        }
        return 0;
    }
    
    public int getSetAssignmentsCount(final int dim){
        int setAssignmentsCounter = 0;
        if(this.containsDim(dim)){
            for (int j = 0; j < this.assignments[dim].length; j++) {
                if(this.assignments[dim][j] != -1){
                    setAssignmentsCounter++;
                }
            }
        }
        return setAssignmentsCounter;
    }
    
    public Boolean isFullyAssigned(final int dim){
        if(!this.containsDim(dim)){
            return null;
        }
        
        return this.getSetAssignmentsCount(dim) == this.getAssignmentsCount();
    }
    
    /**
     * Adds a new assignment entry for a further signal. The given atom indices 
     * will be stored as atom index for each dimension of the signal/spectrum. 
     *
     * @param atomIndicesInStructure
     * @return
     */
    public boolean addAssignment(final int[] atomIndicesInStructure){
        if(atomIndicesInStructure.length != this.getNDim()){
            return false;
        }
        final int[][] extendedAssignments = new int[this.getNDim()][this.getAssignmentsCount()+1];
        for (int dim = 0; dim < this.getNDim(); dim++) {
            for (int i = 0; i < this.getAssignmentsCount(); i++) {
                extendedAssignments[dim][i] = this.getAtomIndex(dim, i);
            }
            extendedAssignments[dim][this.getAssignmentsCount()] = atomIndicesInStructure[dim];
        }
        this.assignments = extendedAssignments;
        
        return true;
    }

    private boolean checkSpectrumIndex(final int dim, final int indexInSpectrum){
       return (indexInSpectrum >= 0) && (indexInSpectrum < this.assignments[dim].length);
    }
    
    private boolean checkInputListSize(final int size){
       return (size == this.getAssignmentsCount());
   } 
    
    @Override
    public Assignment clone() throws CloneNotSupportedException{
        return (Assignment) super.clone();
    }
}
