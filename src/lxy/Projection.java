package lxy;

import java.util.ArrayList;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Box;
import edu.stanford.rsl.conrad.geometry.shapes.simple.StraightLine;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.geometry.transforms.Transform;
import edu.stanford.rsl.conrad.geometry.transforms.Translation;

public class Projection{
	private final double samplingRate = 3.d;
	private double maxTheta, deltaTheta, maxS, deltaS;		
	int maxThetaIdx, maxSIdx;
	
	public Projection(double maxTheta, double deltaTheta, double maxS, double deltaS)
	{
		this.maxTheta = maxTheta;
		this.deltaTheta = deltaTheta;
		this.maxS = maxS;
		this.deltaS = deltaS;
		maxThetaIdx = (int) (maxTheta / deltaTheta) + 1;
		maxSIdx = (int)(maxS / deltaS) + 1;

	}
	
	public Grid2D projectRayDriven(Grid2D object)
	{
		Grid2D sinogram = new Grid2D(maxSIdx, maxThetaIdx);
		sinogram.setSpacing(deltaS, deltaTheta);
		
		Box b = new Box((object.getSize()[0] * object.getSpacing()[0]), (object.getSize()[1] * object.getSpacing()[1]), 2);
		Translation trans = new Translation(-(object.getSize()[0] * object.getSpacing()[0])/2, -(object.getSize()[1] * object.getSpacing()[1])/2, -1);
        b.applyTransform(trans);
        
		double theta = 0;
		for (int thetaIdx = 0; thetaIdx < maxThetaIdx; ++thetaIdx)
		{
			double cosTheta = Math.cos(theta);
			double sinTheta = Math.sin(theta);
			SimpleVector dir = new SimpleVector(-sinTheta, cosTheta, 0d);
			double s = -maxS / 2;
			for(int sIdx = 0; sIdx < maxSIdx; ++sIdx)
			{
				StraightLine ray = new StraightLine(new PointND(s * cosTheta, s * sinTheta, 0d), dir);
				ArrayList<PointND> points = b.intersect(ray);
				if(points.size() == 2)
				{
					PointND start = points.get(0);
					PointND end = points.get(1);
					SimpleVector stepsize = end.getAbstractVector();
					stepsize.subtract(start.getAbstractVector());
					double length = stepsize.normL2();
					int sampleNums = (int)(length * samplingRate);
					stepsize.divideBy(sampleNums);
								
					double sum = 0;
					SimpleVector sample = start.getAbstractVector();
					for(int i = 0; i < sampleNums; ++i)
					{
						double x = (sample.getElement(0) + object.getSize()[0] * object.getSpacing()[0] / 2 ) / object.getSpacing()[0],
								y = (sample.getElement(1) + object.getSize()[1] * object.getSpacing()[1] / 2 )  / object.getSpacing()[1];
						if(x >= 0 && x <= object.getSize()[0] - 1 && y >= 0 && y <= object.getSize()[1] - 1)
							sum += InterpolationOperators.interpolateLinear(object, x, y);
						sample.add(stepsize);
					}
					// normalize by the number of interpolation points
					sum /= samplingRate;
					// write integral value into the sinogram.
					sinogram.setAtIndex(sIdx, thetaIdx, (float)sum);
					
							
				}
				
				s += deltaS;
			}
			theta += deltaTheta;
		}
		
		return sinogram;
		
	}

}
