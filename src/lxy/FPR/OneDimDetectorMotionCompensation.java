package lxy.FPR;

import edu.stanford.rsl.conrad.data.numeric.Grid1D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.utils.CONRAD;
import edu.stanford.rsl.jpop.GradientOptimizableFunction;
import edu.stanford.rsl.jpop.OptimizableFunction;
import edu.stanford.rsl.jpop.OptimizationOutputFunction;

public class OneDimDetectorMotionCompensation implements OptimizableFunction, OptimizationOutputFunction{
	
	boolean debug = true;

	int numberOfProcessingBlocks = 1;
	
	protected double[] oldParameterVector;
	
	Grid1D[] lineIntegrals = null;
	int nViews = 0;
	int nSamples = 0;
	
	public OneDimDetectorMotionCompensation(Grid1D[] integrals) {
		this.lineIntegrals = integrals;
		this.nViews = integrals.length;
		this.nSamples = integrals[0].getSize()[0];
	}

	@Override
	public void setNumberOfProcessingBlocks(int number) {
		numberOfProcessingBlocks = number;
	}

	@Override
	public int getNumberOfProcessingBlocks() {
		return numberOfProcessingBlocks;
	}
	
	@Override
	public double evaluate(double[] x, int block) {
		// shifts x[i] are given in mm!
		double val = 0;
		for(int i = 0; i < nViews-1; i++){
			Grid1D ref = (Grid1D)lineIntegrals[i].clone();
			ref.setOrigin(x[i]);
			
			for(int k = i+1; k < nViews-1; k++){
				Grid1D shifted = (Grid1D)lineIntegrals[k].clone();
				shifted.setOrigin(x[k]);
				for(int j = 0; j < nSamples; j++){
					float refVal = ref.getAtIndex(j);
					double worldC = ref.indexToPhysical(j);
					double shiftedC = shifted.physicalToIndex(worldC);
					float shiftedVal;
					if(shiftedC < 0 || shiftedC >  nSamples-1){
						shiftedVal = 0;
					}else{
						shiftedVal = InterpolationOperators.interpolateLinear(shifted, shiftedC);
					}
					float v1 = refVal-shiftedVal;//Math.min(refVal-shiftedVal, 25);
					val += v1*v1;
				}
			}
		}
		val /= (nViews*(nViews-1)*nSamples);

		if(debug){
			CONRAD.log("--------------------------------------------------------");
			CONRAD.log("Current Cost Fct. Value: " + val);
			CONRAD.log("--------------------------------------------------------");
		}
		return val;
	}

	
	@Override
	public void optimizerCallbackFunction(int currIterationNumber, double[] x,
			double currFctVal, double[] gradientAtX) {
		
		if(currIterationNumber%10==0){
			String line = "[ ";
			for(int i = 0; i < x.length; i++){
				line += x[i] + " ";
			}
			line += "] => "+ currFctVal;
			System.out.println(line);
		}
	}
}