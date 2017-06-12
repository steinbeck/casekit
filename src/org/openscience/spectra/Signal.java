package org.openscience.spectra;

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
