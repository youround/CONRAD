package lxy.ex3;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;

public class reBinning {
	private double maxTheta, deltaTheta, maxS, deltaS;		
	int maxThetaIdx, maxSIdx;
	
	public reBinning(double maxTheta, double deltaTheta, double maxS, double deltaS)
	{
		this.maxTheta = maxTheta;
		this.deltaTheta = deltaTheta;
		this.maxS = maxS;
		this.deltaS = deltaS;
		maxThetaIdx = (int) (maxTheta / deltaTheta) + 1;
		maxSIdx = (int)(maxS / deltaS) + 1;
	}
	
	public Grid2D applyReBinning(Grid2D fanogram, double focalLength) 
	{
		Grid2D sinogram = new Grid2D(maxSIdx, maxThetaIdx);
		sinogram.setSpacing(deltaS, deltaTheta);
		
		int maxBetaIdx = fanogram.getSize()[1];
		int maxTIdx = fanogram.getSize()[0];
		double deltaT = fanogram.getSpacing()[0];
		double deltaBeta = fanogram.getSpacing()[1];
				
		double theta = 0;
		for (int thetaIdx = 0; thetaIdx < maxThetaIdx; ++thetaIdx)
		{
			double s = -maxS / 2;
			for(int sIdx = 0; sIdx < maxSIdx; ++sIdx)
			{
				double gamma = Math.asin(s / focalLength);
				double beta = theta - gamma;
				double t = s / Math.cos(gamma);	
				if (beta < 0)
				{
					beta += Math.PI + 2 * gamma;
					t = -t;				
				}
							
				double betaIdx  = beta / deltaBeta;	
				double tIdx = (t + maxTIdx * deltaT /2) / deltaT;
				
				float val = InterpolationOperators.interpolateLinear(fanogram, tIdx, betaIdx);
				sinogram.setAtIndex(sIdx, thetaIdx, val);				
				
				s += deltaS;
			}
			
			theta += deltaTheta;
			
		}
	
		return sinogram;
	}

}
