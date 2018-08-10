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
   private final int nDim;
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
   private final ArrayList<Integer> equivalences = new ArrayList<>();
  

   public Spectrum(final String[] nuclei) {
       this.nuclei = nuclei;
       this.nDim = this.nuclei.length;
   }
   
   
   public String[] getNuclei(){
       return this.nuclei;
   }
   
   public int getDimCount(){
       return this.nDim;
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
   
   public final boolean setShifts(final ArrayList<Double> shiftList, final int dim){
        if(!this.checkDimension(dim) || (!this.checkInputListSize(shiftList.size()))){
            return false;
        }
        for (int i = 0; i < shiftList.size(); i++) {
            this.setShift(shiftList.get(i), dim, i);
        }
        
        return true;
   }
   
   public final boolean setShift(final double shift, final int dim, final int signalIndex){
        if(!this.checkDimension(dim) || !this.checkSignalIndex(signalIndex)){
            return false;
        }    
        this.getSignal(signalIndex).setShift(shift, dim);
        
        return true;
   }
   
   public int getSignalCount() {
       return this.signals.size();
   }

   /**
    * Adds a Signal ({@link casekit.NMR.model.Signal}) to this Spectrum class object at the end.
     * @param signal
     * @return 
    */
   public boolean addSignal(final Signal signal) {
       if(!this.checkDimCount(signal.getDimCount()) || !this.checkNuclei(signal.getNuclei())){
           return false;
       }
       this.addSignal(signal, null);
       
       return true;
   }
   
   /**
    * Adds a Signal ({@link casekit.NMR.model.Signal}) to this Spectrum class object at given index.
     * @param signal
     * @param index index where to insert the signal, if null the signal will be added at the end of signal list
     * @return 
    */
   public boolean addSignal(final Signal signal, final Integer index) {
       if(!this.checkDimCount(signal.getDimCount()) || !this.checkNuclei(signal.getNuclei())){
           return false;
       }
       // is index valid? if yes then insert it there
       if(this.checkSignalIndex(index)){
           this.signals.add(index, signal);
           this.equivalences.add(index, -1);
       // if not then check for null value and add signal at the end    
       } else if(index == null){
           this.signals.add(signal);
           this.equivalences.add(-1);
       // no valid index value, nothing to insert or add in spectrum    
       } else {
           return false;
       }
       
       return true;
   }

   public boolean removeSignal(final Signal signal){
       return this.removeSignal(this.getSignalIndex(signal));
   }
   
   public boolean removeSignal(final int signalIndex){
       if(!this.checkSignalIndex(signalIndex)){
           return false;
       }
       this.signals.remove(signalIndex);
       this.equivalences.remove(signalIndex);
       
       return true;
   }
   
   private boolean checkSignalIndex(final Integer signalIndex){
       return (signalIndex != null) && (signalIndex >= 0) && (signalIndex < this.getSignalCount());
   }
   
   public boolean checkDimension(final int dim){
       return (dim >= 0) && (dim < this.nDim);
   }
   
   private boolean checkInputListSize(final int size){
       return (size == this.getSignalCount());
   }
   
   private boolean checkDimCount(final int ndim){
       return ndim == this.getDimCount();
   }
   
   private boolean checkNuclei(final String[] nuclei){
       return nuclei == this.getNuclei();
   }
   
   /**
    * Returns an NMRSignal at position number in the List
     * @param signalIndex
     * @return 
    */
   public Signal getSignal(final int signalIndex) {
       if(!this.checkSignalIndex(signalIndex)){
           return null;
       }
       
       return this.signals.get(signalIndex);
   }
   
   public ArrayList<Double> getIntensities(){
       final ArrayList<Double> intensities = new ArrayList<>();
       for (Signal sig : this.signals) {
           intensities.add(sig.getIntensity());
       }
       
       return intensities;
   }
   
   public Double getIntensity(final int signalIndex){
       if(!this.checkSignalIndex(signalIndex)){
           return null;
       }
       
       return this.getSignal(signalIndex).getIntensity();
   }
   
   public boolean setIntensities(final ArrayList<Double> intensities){
       if(!this.checkInputListSize(intensities.size())){
           return false;
       }
       for (int s = 0; s < this.getSignalCount(); s++) {
           this.setIntensity(intensities.get(s), s);
       }
       
       return true;
   }
   
   public boolean setIntensity(final double intensity, final int signalIndex){
       if(!this.checkSignalIndex(signalIndex)){
           return false;
       }
       this.getSignal(signalIndex).setIntensity(intensity);
       
       return true;
   }
   
   public ArrayList<Double> getShifts(final int dim){
       final ArrayList<Double> shifts = new ArrayList<>();
       if(!this.checkDimension(dim)){
           return shifts;
       }
       for (final Signal sig : this.signals) {
           shifts.add(sig.getShift(dim));
       }
       
       return shifts;
   }
   
   public Double getShift(final int SignalIndex, final int dim){
       if(!this.checkSignalIndex(SignalIndex)){
           return null;
       }
       
       return this.getSignal(SignalIndex).getShift(dim);
   }
   
   public boolean setMultiplicities(final ArrayList<String> multiplicities){
       if(!this.checkInputListSize(multiplicities.size())){
           return false;
       }
       for (int s = 0; s < this.getSignalCount(); s++) {
           this.setMultiplicity(multiplicities.get(s), s);
       }
       
       return true;
   }
   
   public boolean setMultiplicity(final String multiplicity, final int signalIndex){
       if(!this.checkSignalIndex(signalIndex)){
           return false;
       }
       this.getSignal(signalIndex).setMultiplicity(multiplicity);
       
       return true;
   }
   
   public ArrayList<Signal> getSignals(){
       return this.signals;
   }
   
   public ArrayList<Integer> getEquivalences(){
       return this.equivalences;
   }
   
   public Integer getEquivalence(final int signalIndex){
       if(!this.checkSignalIndex(signalIndex)){
           return null;
       }
       
       return this.equivalences.get(signalIndex);
   }
   
   public boolean setEquivalence(final int signalIndex, final int equivalentSignalIndex){
       if(!this.checkSignalIndex(signalIndex) || !this.checkSignalIndex(equivalentSignalIndex)){
           return false;
       }
       this.equivalences.set(signalIndex, equivalentSignalIndex);
       
       return true;
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
   

   /**
    * Returns the signal closest to the shift sought. If no Signal is found within the interval
    * defined by pickprecision, null is returned.
     * @param shift
     * @param dim
     * @param pickPrecision
     * @return 
    */
   public int pickClosestSignal(final double shift, final int dim, final double pickPrecision) {       
       int matchIndex = -1;
       if(!this.checkDimension(dim)){
           return matchIndex;
       }
       double diff = pickPrecision;
       for (int s = 0; s < this.getSignalCount(); s++) {
           if (Math.abs(this.getShift(s, dim) - shift) < diff) {
               diff = Math.abs(this.getShift(s, dim) - shift);
               matchIndex = s;
           }
       }       
       
       return matchIndex;
   }

   /**
    * Returns a List with signal indices within the interval defined by pickPrecision. If none is found
    * an empty ArrayList is returned.
     * @param shift
     * @param dim
     * @param pickPrecision
     * @return 
    */
   public ArrayList<Integer> pickSignals(final double shift, final int dim, final double pickPrecision) {
       final ArrayList<Integer> pickedSignals = new ArrayList<>();
       if(!this.checkDimension(dim)){
           return pickedSignals;
       }
       for (int s = 0; s < this.getSignalCount(); s++) {
           if (Math.abs(this.getShift(s, dim) - shift) < pickPrecision) {
               pickedSignals.add(s);
           }
       }
       
       return pickedSignals;
   }
   
   
   
}
