/*
* This class was copied and modified from NMRSpectrum class in casekit.model package (by Christoph Steinbeck)
*/


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

import java.util.ArrayList;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Spectrum {
                                                  
   /**
    * An arbitrary name or description that can be assigned to this spectrum for identification purposes.
    */
   private String description;
   /**
    * An arbitrary name to identify the type of this spectrum, like COSY, NOESY, HSQC, etc. I
    * decided not to provide static Strings with given experiment type since the there are
    * numerous experiments yielding basically identical information having different names
    */
   private String specType;
   
   /**
    * Declares how many axes are in involved in this spectrum.
    */
   private final int ndim;
   /**
    * The nuclei of the different axes.
    */
   private final String nuclei[];
   /**
    * The proton frequency of the spectrometer used to record this spectrum.
    */
   private Float spectrometerFrequency;
   private String solvent;
   private String standard;
   
   private final ArrayList<Signal> signals = new ArrayList<>();
  

   public Spectrum(final String[] nuclei, final ArrayList<Double>[] shiftLists, final ArrayList<Double> intensities) {
       this.nuclei = nuclei;
       this.ndim = this.nuclei.length;
       this.addSignals(shiftLists, intensities);
   }
   
   public Spectrum(final String[] nuclei, final ArrayList<Double>[] shiftLists) {
       this.nuclei = nuclei;
       this.ndim = this.nuclei.length;
       this.addSignals(shiftLists, null);
   }
   
   public Spectrum(final String[] nuclei) {
       this.nuclei = nuclei;
       this.ndim = this.nuclei.length;
   }
   
   public String[] getNuclei(){
       return this.nuclei;
   }
   
   public int getDimCount(){
       return this.ndim;
   }
   
   public void setSpecType(final String specType){
       
        this.specType = specType;
   }
   
   public String getSpecType(){
        return this.specType;
   }
   
   public void setSpecDescription(final String description){
       this.description = description;
   }
   
   public String getSpecDescription(){
        return this.description;
   }
   
   public final void addSignals(final ArrayList<Double>[] shiftLists, final ArrayList<Double> intensities){
       // assumes the same number of shift lists as set dimension number
       if(shiftLists.length != this.ndim){
           System.err.println("Unequal number of nuclei (dimension) and shift lists!!!");      
           return;
       }
       // assumes that the shift lists have the same number of entries
       int prevShiftListSize = shiftLists[0].size();
       for (int i = 0; i < shiftLists.length; i++) {
           if(shiftLists[i].size() != prevShiftListSize){
                System.err.println("Unequal number of shifts in " + (i+1) + " shift list!!!");      
                return;
           }
           if(intensities != null && shiftLists[i].size() != intensities.size()){
                System.err.println("Unequal number of shifts in shift list " + (i+1) + " and intensities number!!!");
                return;
           }
       }
       Double[] shifts;
       for (int row = 0; row < shiftLists[0].size(); row++) {
           shifts = new Double[this.ndim];
           for (int col = 0; col < this.ndim; col++) {
               shifts[col] = shiftLists[col].get(row);
           }
           if(intensities != null){
               this.addSignal(new casekit.NMR.model.Signal(this.nuclei, shifts, intensities.get(row)));
           } else {
               this.addSignal(new casekit.NMR.model.Signal(this.nuclei, shifts));
           }
       }
   }
   
   public int getSignalCount() {
       return this.signals.size();
   }

   /**
    * Adds a Signal ({@link casekit.NMR.model.Signal}) to this Spectrum class object.
     * @param signal
    */
   public void addSignal(final Signal signal) {
       this.signals.add(signal);
   }

   public void removeSignal(final Signal signal){
       this.signals.remove(this.getSignalIndex(signal));
   }
   
   public void removeSignal(final int signalIndex){
       this.signals.remove(signalIndex);
   }
   
   /**
    * Returns an NMRSignal at position number in the List
     * @param signalIndex
     * @return 
    */
   public Signal getSignal(final int signalIndex) {
       return this.signals.get(signalIndex);
   }
   
   public ArrayList<Double> getIntensities(){
       final ArrayList<Double> intensities = new ArrayList<>();
       for (Signal sig : this.signals) {
           intensities.add(sig.getIntensity());
       }
       
       return intensities;
   }
   
   public ArrayList<Double> getShiftsByDim(final int dim){
       final ArrayList<Double> shifts = new ArrayList<>();
       for (final Signal sig : this.signals) {
           shifts.add(sig.getShift(dim));
       }
       
       return shifts;
   }
   
   public int getAssignmentAtomIndexByDim(final int signalIndex, final int dim){
       return this.getSignal(signalIndex).getAssignedAtomIndex(dim);
   }
   
   public ArrayList<Integer> getAssignmentAtomIndicesByDim(final int dim){
       final ArrayList<Integer> indices = new ArrayList<>();
       for (final Signal sig : this.signals) {
           indices.add(sig.getAssignedAtomIndex(dim));
       }
       
       return indices;
   }

   /**
    * Returns the position of an NMRSignal the List
     * @param signal
     * @return 
    */
   public int getSignalIndex(final Signal signal) {
       for (int f = 0; f < this.signals.size(); f++) {
           if (this.signals.get(f) == signal) {
               return f;
           }
       }
       return -1;
   }

   public void setSpectrometerFrequency(final Float sf) {
       this.spectrometerFrequency = sf;
   }

   public float getSpectrometerFrequency() {
       return spectrometerFrequency;
   }

   public void setSolvent(final String solvent) {
       this.solvent = solvent;
   }

   public String getSolvent() {
       return solvent;
   }

   public void setStandard(final String standard) {
       this.standard = standard;
   }

   public String getStandard() {
       return standard;
   }
   
   public void setAssignedAtomIndicesByDim(final ArrayList<Integer> indices, final int dim){
       if(indices.size() != this.getShiftsByDim(0).size()){
           return;
       }
       for (int i = 0; i < indices.size(); i++) {
           this.getSignal(i).setAssignedAtomIndex(indices.get(i), dim);
       }
   }

   /**
    * Returns the signal closest to the shift sought. If no Signal is found within the interval
    * defined by pickprecision, null is returned.
     * @param shift
     * @param dim
     * @param pickPrecision
     * @return 
    */
   public Signal pickClosestSignal(final Double shift, final int dim, final double pickPrecision) {
       
       int thisPosition = -1;
       double diff = pickPrecision;
       /*
       * Now we search dimension dim for the chemical shift.
       */
       for (int f = 0; f < this.signals.size(); f++) {
           if (Math.abs((this.signals.get(f)).getShift(dim) - shift) < diff) {
               diff = Math.abs((this.signals.get(f)).getShift(dim) - shift);
               thisPosition = f;
           }
       }
       if(thisPosition >= 0){
           this.signals.get(thisPosition);
       }
       
       return null;
   }

   /**
    * Returns a List with signals within the interval defined by pickPrecision. If none is found
    * an empty ArrayList is returned.
     * @param shift
     * @param dim
     * @param pickPrecision
     * @return 
    */
   public ArrayList<Signal> pickSignals(final double shift, final int dim, final double pickPrecision) {
       final ArrayList<casekit.NMR.model.Signal> pickedSignals = new ArrayList<>();
       /*
        * Now we search dimension dim for the chemical shift.
        */
       for (final Signal sig : this.signals) {
           if (Math.abs(sig.getShift(dim) - shift) < pickPrecision) {
               pickedSignals.add(sig);
           }
       }
       return pickedSignals;
   }
   
   public ArrayList<Double>[] getShifts(){
       final ArrayList<Double>[] shiftsList = new ArrayList[this.ndim];
       for (int d = 0; d < this.ndim; d++) {
           shiftsList[d] = this.getShiftsByDim(d);
       }
       
       return shiftsList;
   }
   
   
}
