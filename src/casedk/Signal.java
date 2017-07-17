/* 
* This Open Source Software is provided to you under the MIT License
 * Refer to doc/mit.license or https://opensource.org/licenses/MIT for more information
 * 
 * Copyright (c) 2017, Christoph Steinbeck 
 */

package casedk;

public class Signal {

	Double shift = null;
	Integer mult = null;
	
	public Signal() {
		// TODO Auto-generated constructor stub
	}
	
	Signal(double shift)
	{
		setShift(shift);
	}

	Signal(double shift, int mult)
	{
		setShift(shift);
		setMult(mult);
	}

	public Integer getMult() {
		return mult;
	}
	public void setMult(Integer mult) {
		this.mult = mult;
	}
	
	public Double getShift() {
		return shift;
	}

	public void setShift(Double shift) {
		this.shift = shift;
	}

}
