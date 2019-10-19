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
     * Sets an assignment with value for an index position.
     *
     * @param dim
     * @param index
     * @param assignment
     * @return
     */
    public boolean setAssignment(final int dim, final int index, final int assignment){
        if(!this.containsDim(dim) || !this.checkIndex(dim, index)){
            return false;
        }
        this.assignments[dim][index] = assignment;
        
        return true;
    }
    
    public boolean setAssignments(final int dim, final List<Integer> assignments){
        if(!this.containsDim(dim) || !this.checkInputListSize(assignments.size())){
            return false;
        }
        for (int i = 0; i < this.getAssignmentsCount(); i++) {
            this.setAssignment(dim, i, assignments.get(i));
        }
        
        return true;
    }
    
    public Integer getAssignment(final int dim, final int index){
        if(!this.containsDim(dim) || !this.checkIndex(dim, index)){
            return null;
        }

        return this.assignments[dim][index];
    }
    
    public Integer getIndex(final int dim, final int assignment){
        if(!this.containsDim(dim)){
            return null;
        }
        for (int index = 0; index < this.assignments[dim].length; index++) {
            if(this.getAssignment(dim, index) == assignment){
                return index;
            }
        }
        
        return -1;
    }
    
    public List<Integer> getAssignments(final int dim){
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
     * Adds a new assignment entry (index), e.g. for a new signal. The given assignment indices
     * will be stored for each dimension of the new assignment entry (index).
     *
     * @param assignment assignment indices to store in each dimension of new assignment entry
     *
     * @return
     */
    public boolean addAssignment(final int[] assignment){
        if(!this.compareNDim(assignments.length)){
            return false;
        }
        final int[][] extendedAssignments = new int[this.getNDim()][this.getAssignmentsCount() + 1];
        for (int dim = 0; dim < this.getNDim(); dim++) {
            for (int i = 0; i < this.getAssignmentsCount(); i++) {
                extendedAssignments[dim][i] = this.getAssignment(dim, i);
            }
            extendedAssignments[dim][this.getAssignmentsCount()] = assignment[dim];
        }
        this.assignments = extendedAssignments;
        
        return true;
    }

    public boolean removeAssignment(final int index){
        if(!this.checkIndex(0, index)){
            return false;
        }
        final int[][] reducedAssignments = new int[this.getNDim()][this.getAssignmentsCount() - 1];
        int nextIndexToInsertCounter = 0;
        for (int i = 0; i < this.getAssignmentsCount(); i++) {
            if(i == index){
                continue;
            }
            for (int dim = 0; dim < this.getNDim(); dim++) {
                reducedAssignments[dim][nextIndexToInsertCounter] = this.getAssignment(dim, i);
            }
            nextIndexToInsertCounter++;
        }
        this.assignments = reducedAssignments;

        return true;
    }

    private boolean checkIndex(final int dim, final int index){
       return (index >= 0) && (index < this.assignments[dim].length);
    }
    
    private boolean checkInputListSize(final int size){
       return (size == this.getAssignmentsCount());
   } 
    
    @Override
    public Assignment clone() throws CloneNotSupportedException{
        return (Assignment) super.clone();
    }
}
