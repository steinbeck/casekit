/*
* This class was adopted and modified from an earlier version by Christoph Steinbeck
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
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;

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
   private Double spectrometerFrequency;
   private String solvent;
   private String standard;
   private final ArrayList<Signal> signals;
   private int signalCount;
   private final ArrayList<Integer> equivalences;
   private  ArrayList<Integer>[] equivalentSignals;
  

   public Spectrum(final String[] nuclei) {
       this.nuclei = nuclei;
       this.nDim = this.nuclei.length;
       this.signals = new ArrayList<>();
       this.signalCount = 0;
       this.equivalences = new ArrayList<>();
       this.equivalentSignals = new ArrayList[]{};
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
   
   public final boolean setShift(final Double shift, final int dim, final int signalIndex){
        if(!this.checkDimension(dim) || !this.checkSignalIndex(signalIndex)){
            return false;
        }    
        this.getSignal(signalIndex).setShift(shift, dim);
        
        return true;
   }
   
   public int getSignalCount() {
       return this.signalCount;
   }
   
    /**
     * Adds a list of signals to this spectrum.
     *
     * @param signals list of signals to add
     * @return
     */
   public boolean addSignals(final ArrayList<Signal> signals){
       for (final Signal signal : signals) {
           if (!this.checkDimCount(signal.getDimCount()) || !this.checkNuclei(signal.getNuclei())) {
               return false;
           }
       }
       for (final Signal signal : signals) {
           this.addSignal(signal);
       }
       
       return true;
   }
   
   /**
    * Adds a signal to this spectrum.
    * 
    * @param signal signal to add
    * @return 
    */
   public boolean addSignal(final Signal signal) {       
       return this.addSignal(signal, -1);
   }
   
   /**
    * Adds a signal to this spectrum and stores an equivalent signal index.
    * 
    * @param signal signal to add
    * @param equivalentSignalIndex index of equivalent signal in this spectrum
    * @return 
    */
   public boolean addSignal(final Signal signal, final int equivalentSignalIndex) {
       if((signal == null) || !this.checkDimCount(signal.getDimCount()) || !this.checkNuclei(signal.getNuclei())){
           return false;
       }
       // add signal at the end of signal list  
       if(this.signals.add(signal)){
           this.signalCount++;
           this.equivalences.add(equivalentSignalIndex);
           this.updateEquivalentSignalClasses();
           
           return true;
       }       
       
       return false;
   }

   public boolean removeSignal(final Signal signal){
       return this.removeSignal(this.getSignalIndex(signal));
   }
   
   public boolean removeSignal(final int signalIndex){
        if(!this.checkSignalIndex(signalIndex)){
           return false;
        }
        if(this.signals.remove(signalIndex) != null){
            this.signalCount--;
            this.equivalences.remove(signalIndex);
            this.updateEquivalentSignalClasses();   
            
            return true;
        }       
                  
        return false;
   }
   
   private boolean checkSignalIndex(final Integer signalIndex){       
       return (signalIndex != null) && (signalIndex >= 0) && (signalIndex < this.getSignalCount());
   }
   
    /**
     * Checks whether the input dimension exists in this spectrum or not.
     *
     * @param dim
     * @return
     */
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
       return Arrays.equals(nuclei, this.getNuclei());
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
       
       try {
           return this.signals.get(signalIndex);
       } catch (Exception e) {
           return null;
       }
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
   
    public ArrayList<String> getMultiplicities() {
        final ArrayList<String> multiplicities = new ArrayList<>();        
        for (final Signal sig : this.signals) {
            multiplicities.add(sig.getMultiplicity());
        }

        return multiplicities;
    }

    public String getMultiplicity(final int SignalIndex) {
        if (!this.checkSignalIndex(SignalIndex)) {
            return null;
        }

        return this.getSignal(SignalIndex).getMultiplicity();
    }
   
   public ArrayList<Signal> getSignals(){
       return this.signals;
   }
   
   public Boolean hasEquivalences(final int signalIndex){
       if(!this.checkSignalIndex(signalIndex)){
           return null;
       }
       
       return (this.getEquivalence(signalIndex) != -1) || (this.getEquivalences().contains(signalIndex));
   }
   
   private ArrayList<Integer> searchEquivalentSignals(final int signalIndex){       
       if(!this.checkSignalIndex(signalIndex)){
           return null;
       }
       final ArrayList<Integer> equivalentSignalIndices = new ArrayList<>();
       // case 1: signal was first input signal (root) of an equivalence class and is actually not knowing any of its equivalences; collect all equivalent signals
       if(this.getEquivalence(signalIndex) == -1){
           for (int i = 0; i < this.getEquivalences().size(); i++) {
               if((this.getEquivalences().get(i) != -1) && (this.getEquivalences().get(i) == signalIndex)) {
                   equivalentSignalIndices.add(i);
               }               
           }
       } else {
           // case 2: signal was not the first input signal of that equivalent class; store the class root signal
           equivalentSignalIndices.add(this.getEquivalences().get(signalIndex));
       }
       // check all stored signals for further equivalent signals (i.e. for the added root signal in case 2)
       for (int i = 0; i < equivalentSignalIndices.size(); i++) {
           for (int j = 0; j < this.getEquivalences().size(); j++) {
               // do not store the own signal index in own equ. signal class
               if(j == signalIndex){
                   continue;
               }
               if ((this.getEquivalences().get(j) != -1) 
                       && (Integer.compare(this.getEquivalences().get(j), equivalentSignalIndices.get(i)) == 0)
                       && !equivalentSignalIndices.contains(j)) {
                   equivalentSignalIndices.add(j);
               }
           }           
       }
       
       return equivalentSignalIndices;
   }
   
   private void updateEquivalentSignalClasses(){     
       this.equivalentSignals = new ArrayList[this.getSignalCount()];
       for(int i = 0; i < this.getSignalCount(); i++) {
           this.equivalentSignals[i] = this.searchEquivalentSignals(i);
       }
   }
   
    /**
     * Returns equivalent signals for requested signal. 
     *
     * @param signalIndex
     * @return
     */
    public ArrayList<Integer> getEquivalentSignals(final int signalIndex){
       if(!this.checkSignalIndex(signalIndex)){
           return null;
       }
       
       return this.equivalentSignals[signalIndex];
   }
    
    /**
     * Returns a hashmap of equivalent signal classes. 
     * The key set of that hashmap is just a numerical class index and is not 
     * belonging to any signal.
     *
     * @return
     */
    public HashMap<Integer, ArrayList<Integer>> getEquivalentSignalClasses(){
        this.updateEquivalentSignalClasses();
        // create a new HashMap object to return, containing the key signal index to have a full equivalent signal class
        final HashMap<Integer, ArrayList<Integer>> equivalentSignalClasses = new HashMap<>();
        final HashSet<Integer> storedSignalIndices = new HashSet<>();
        for (int i = 0; i < this.getSignalCount(); i++) {
            if (!storedSignalIndices.contains(i)) {               
                equivalentSignalClasses.put(equivalentSignalClasses.size(), new ArrayList<>(this.equivalentSignals[i]));
                equivalentSignalClasses.get(equivalentSignalClasses.size() - 1).add(i);
                storedSignalIndices.addAll(equivalentSignalClasses.get(equivalentSignalClasses.size() - 1));
            }
        }

        return equivalentSignalClasses;
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
   
   public boolean setEquivalence(final int signalIndex, final int isEquivalentToSignalIndex){
       if(!this.checkSignalIndex(signalIndex) || !this.checkSignalIndex(isEquivalentToSignalIndex)){
           return false;
       }
       this.equivalences.set(signalIndex, isEquivalentToSignalIndex);
       this.updateEquivalentSignalClasses();
       
       return true;
   }

   /**
    * Returns the position of an NMRSignal the List
    * @param signal
    * @return 
    */
   public int getSignalIndex(final Signal signal) {
       for (int s = 0; s < this.signals.size(); s++) {
           if (this.signals.get(s) == signal) {
               return s;
           }
       }
       return -1;
   }

   public void setSpectrometerFrequency(final Double sf) {
       this.spectrometerFrequency = sf;
   }

   public Double getSpectrometerFrequency() {
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
    * Returns the signal index closest to the given shift. If no Signal is found within the interval
    * defined by pickprecision, -1 is returned.
    * @param shift query shift
    * @param dim dimension in spectrum to look in
    * @param pickPrecision tolerance value for search window
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
    * Returns a list of signal indices within the interval defined by 
    * pickPrecision. That list is sorted by the distances to the query shift.    
    * If none is found an empty ArrayList is returned.
    * @param shift query shift
    * @param dim dimension in spectrum to look in
    * @param pickPrecision tolerance value for search window
    * @return 
    */
   public ArrayList<Integer> pickSignals(final Double shift, final int dim, final double pickPrecision) {
       final ArrayList<Integer> pickedSignals = new ArrayList<>();
       if(!this.checkDimension(dim)){
           return pickedSignals;
       }
       for (int s = 0; s < this.getSignalCount(); s++) {
           if (Math.abs(this.getShift(s, dim) - shift) <= pickPrecision) {
               pickedSignals.add(s);
           }
       }
       // sort signal indices by distance to query shift
       pickedSignals.sort(new Comparator<Integer>() {
           @Override
           public int compare(final Integer pickedSignalIndex1, final Integer pickedSignalIndex2) {
               return Double.compare(Math.abs(shift - getShift(pickedSignalIndex1, dim)), Math.abs(shift - getShift(pickedSignalIndex2, dim)));
           }
       });
       
       return pickedSignals;
   }
   
   public Spectrum getClone() {
       final Spectrum clone = new Spectrum(this.nuclei);       
       for (int i = 0; i < this.getSignalCount(); i++) {
           clone.addSignal(this.getSignal(i).getClone(), this.getEquivalence(i));
       }              
       clone.setSpecDescription(this.description);
       clone.setSolvent(this.solvent);
       clone.setSpecType(this.specType);
       clone.setSpectrometerFrequency(this.spectrometerFrequency);
       clone.setStandard(this.standard);
       
       return clone;
   }
   
}
