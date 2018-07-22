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

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Signal {
    
    private final String element;
    private final String multiplicity;
    private final double shift;
    private final Double intensity;
    
    public Signal(final String elem, final double shift, final String mult, final Double intens){
        this.element = elem;
        this.shift = shift;
        this.multiplicity = mult;
        this.intensity = intens;
    }
    
    public String getElement(){
        return this.element;
    }
    
    public double getShift(){
        return this.shift;
    }
    
    public String getMultiplicity(){
        return this.multiplicity;
    }
    
    public Double getIntensity(){
        return this.intensity;
    }
}
