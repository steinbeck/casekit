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
package NMR;

import java.util.ArrayList;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Spectrum extends ArrayList<NMR.Signal>{

   /**
    * An arbitrary name that can be assigned to this spectrum for identification purposes.
    */
   private String name;
   /**
    * An arbitrary name to identify the type of this spectrum, like COSY, NOESY, HSQC, etc. I
    * decided not to provide static Strings with given experiment type since the there are
    * numerous experiments yielding basically identical information having different names
    */
   private String specType;
   
//   /**
//    * Not yet clear if this is needed.
//    */
//   private float[] pickPrecision;
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
//   protected transient EventListenerList changeListeners = new EventListenerList();

   public Spectrum(final String[] nuclei, final ArrayList<Double>[] shiftLists) {
       this.nuclei = nuclei;
       this.ndim = this.nuclei.length;
       this.setSignals(shiftLists);
   }
   
   public Spectrum(final String[] nuclei) {
       this.nuclei = nuclei;
       this.ndim = this.nuclei.length;
   }
   
   
   public void setSignals(final ArrayList<Double>[] shiftLists){
       if(shiftLists.length != this.ndim){
           System.err.println("Unequal number of nuclei and shift lists!!!");      
           return;
       }
       Double[] shifts;
       // assumes that the shift lists have the same number of entries
       for (int row = 0; row < shiftLists[0].size(); row++) {
           shifts = new Double[this.ndim];
           for (int col = 0; col < this.ndim; col++) {
               shifts[col] = shiftLists[col].get(row);
           }
           this.add(new NMR.Signal(this.nuclei, shifts, null, null, null));
       }
   }
   

   /**
    * Return the number of individual frequencies in the heteroatom shift list, which should be
    * equal or smaller than the number of respective atoms
     * @return 
    */
   public int getSignalNumber() {
       return this.size();
   }

   /**
    * Adds an NMRSignal to the NMRSpectrum.
     * @param signal
    */
   public void addSignal(Signal signal) {
       this.add(signal);
//       this.updateShiftLists();s
   }

//   /**
//    * Creates an empty signal with correct dimension
//    */
//   public void newSignal() {
//       System.out.println("nucleus: " + nucleus.length + nucleus[0]);
//       add(new NMRSignal(nucleus));
//       updateShiftLists();
//   }

   /**
    * Returns an NMRSignal at position number in the List
     * @param signalIndex
     * @return 
    */
   public Signal getSignal(final int signalIndex) {
       return this.get(signalIndex);
   }
   
   public ArrayList<Double> getShiftsByDim(final int dim){
       final ArrayList<Double> signals = new ArrayList<>();
       for (final Signal sig : this) {
           signals.add(sig.getShift(dim));
       }
       
       return signals;
   }

   /**
    * Returns the position of an NMRSignal the List
     * @param signal
     * @return 
    */
   public int getSignalIndex(final Signal signal) {
       for (int f = 0; f < this.size(); f++) {
           if (this.get(f) == signal) {
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
   public Signal pickClosestSignal(final Double shift, final int dim, final double pickPrecision) {
       
       int thisPosition = -1;
       double diff = pickPrecision;
       /*
       * Now we search dimension dim for the chemical shift.
       */
       for (int f = 0; f < this.size(); f++) {
           if (Math.abs((this.get(f)).getShift(dim) - shift) < diff) {
               diff = Math.abs((this.get(f)).getShift(dim) - shift);
               thisPosition = f;
           }
       }
       if(thisPosition >= 0){
           this.get(thisPosition);
       }
       
       return null;
   }

   /**
    * Returns a List with signals within the interval defined by pickPrecision. If none is found
    * an empty List is returned.
     * @param shift
     * @param dim
     * @param pickPrecision
     * @return 
    */
   public ArrayList<Signal> pickSignals(final double shift, final int dim, final double pickPrecision) {
       final ArrayList<NMR.Signal> pickedSignals = new ArrayList<>();
       /*
        * Now we search dimension dim for the chemical shift.
        */
       for (final Signal sig : this) {
           if (Math.abs(sig.getShift(dim) - shift) < pickPrecision) {
               pickedSignals.add(sig);
           }
       }
       return pickedSignals;
   }

//   /**
//    * Extracts a list of unique shifts from the list of cross signals and sorts them. This is to
//    * define the column and row headers for tables.
//    */
//   protected void updateShiftLists() {
//       Double shift; NMR.Signal nmrSignal;
//       for (int i = 0; i < this.size(); i++) {
//           nmrSignal = this.get(i);
//           for (int d = 0; d < nmrSignal.getDim(); d++) {
//               shift = nmrSignal.getShift(d);
//               if (!this.get(d).contains(shift)) {
//                   shiftList.get(d).add(shift);
//               }
//           }
//       }
//   }

   /**
    * Creates a 2D matrix of booleans, that models the set of crosspeaks in the 2D NMR spectrum.
    * The dimensions are taken from hetAtomShiftList and protonShiftList, which again are
    * produced by updateShiftLists based a collection of 2D nmrSignals
    * <p/>
    * private void createMatrix(){ boolean found; float het, prot; int hetPos, protPos;
    * hetCorMatrix = new boolean[hetAtomShiftList.length][protonShiftList.length]; for (int f =
    * 0; f < size(); f++){ HetCorNMRSignal hetCorSignal = (HetCorNMRSignal)elementAt(f); prot =
    * hetCorSignal.shift[NMRSignal.SHIFT_PROTON]; het =
    * hetCorSignal.shift[NMRSignal.SHIFT_HETERO]; found = false; hetPos =
    * isInShiftList(hetAtomShiftList, het, hetAtomShiftList.length); if (hetPos >= 0){ protPos =
    * isInShiftList(protonShiftList, prot, protonShiftList.length); if ( protPos >= 0){ found =
    * true; hetCorMatrix[hetPos][protPos] = true; } } } }
    */
   public void report() {
       String s = "";
       System.out.println("Report for nmr spectrum " + name + " of type "
               + specType + ": ");
       for (int i = 0; i < this.size(); i++) {
           System.out.println("ShiftList for dimension " + (i + 1) + ":");
           for (int d = 0; d < this.get(i).getDim(); d++) {
               s += this.get(i).getShift(d) + "; ";
           }
           System.out.println(s + "\n");
           s = "";
       }

   }
   
   public String[] getNuclei(){
       return this.nuclei;
   }
}
