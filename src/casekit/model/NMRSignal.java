package casekit.model;

/* NMRSignal.java
*
* Copyright (C)   Dr. Christoph Steinbeck
*
* Contact: christoph.steinbeck@uni-jena.de
*
* This software is published and distributed under MIT license
*/


/**
* A class to store the properties of a single N-dimensional NMR signal
*/

public class NMRSignal {

    int dim;

   /**
    * Am array of doubles to store the chemical shift of
    */
   public float shift[];
   public String[] nucleus;

   /* Signal intensity in arbitrary values */
   public float intensity;

   public int phase;
   public static int DIM_ONE = 1, DIM_TWO = 2, DIM_THREE = 3, DIM_FOUR = 4;
   public static int SHIFT_PROTON = 0, SHIFT_HETERO = 1;
   public static int PHASE_NEGATIVE = 2, PHASE_POSITIVE = 1, PHASE_NONE = 0;
   public static String[] PHASENAMES = {"NONE", "POSITIVE", "NEGATIVE"};

   public NMRSignal(String[] nucleus) {
       this.dim = nucleus.length;
       this.shift = new float[dim];
       this.nucleus = nucleus;
       for (int f = 0; f < dim; f++)
           shift[f] = 0;
       intensity = 1;
       phase = PHASE_POSITIVE;
   }

   public NMRSignal(String[] nucleus, float[] shift, float intensity, int phase) {
       this.dim = nucleus.length;
       this.shift = shift;
       this.nucleus = nucleus;
       this.intensity = intensity;
       this.phase = phase;
   }

   public void setShift(float sshift, String nnucleus) {
       for (int f = 0; f < nucleus.length; f++) {
           if (nucleus[f].equals(nnucleus)) {
               shift[f] = sshift;
               break;
           }
       }
   }

   public void setShift(float sshift, int dim) {
       shift[dim] = sshift;
   }

   public float getShift(String nnucleus) {
       for (int f = 0; f < nucleus.length; f++) {
           if (nucleus[f].equals(nnucleus)) {
               return shift[f];
           }
       }

       return Float.MAX_VALUE;

   }

   public float getShift(int dim) {
       return shift[dim];
   }

   
   public String toString() {
       String s = "";
       s += dim + " -dimensional NMRSignal for nuclei ";
       for (int f = 0; f < nucleus.length; f++)
           s += nucleus[f] + "; ";
       s += "\nShiftlist: ";
       for (int f = 0; f < shift.length; f++)
           s += shift[f] + "; ";
       s += "\n\n";
       return s;
   }

}