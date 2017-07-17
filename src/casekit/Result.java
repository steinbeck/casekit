/* 
* This Open Source Software is provided to you under the MIT License
 * Refer to doc/mit.license or https://opensource.org/licenses/MIT for more information
 * 
 * Copyright (c) 2017, Christoph Steinbeck 
 */

package casekit;

import org.openscience.cdk.interfaces.IAtomContainer;

public class Result {
	public IAtomContainer ac;
	public double score;

	public Result(IAtomContainer ac, double score) {
		super();
		this.ac = ac;
		this.score = score;
	} 

		public IAtomContainer getAc() {
			return ac;
		}
		public void setAc(IAtomContainer ac) {
			this.ac = ac;
		}
		public double getScore() {
			return score;
		}
	
		public void setScore(double score) {
			this.score = score;
		}

		
}