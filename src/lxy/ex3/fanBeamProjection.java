package lxy.ex3;

import java.util.ArrayList;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Box;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.geometry.shapes.simple.StraightLine;
import edu.stanford.rsl.conrad.geometry.transforms.Translation;
import edu.stanford.rsl.conrad.numerics.SimpleVector;

public class fanBeamProjection {
	
	final double samplingRate = 3.d; //samples per mm
	private double focalLength, maxBeta, deltaBeta, maxT, deltaT;
	private int maxTIdx, maxBetaIdx;
	
	public fanBeamProjection(double focalLength, double maxBeta, 
			double deltaBeta, double maxT, double deltaT) 
	{
		this.focalLength = focalLength;
		this.maxBeta = maxBeta;
		this.maxT = maxT;
		this.deltaBeta = deltaBeta;
		this.deltaT = deltaT;
		this.maxBetaIdx = (int) (maxBeta / deltaBeta ) +1;
		this.maxTIdx = (int) (maxT / deltaT) + 1 ;
	}
	
	public Grid2D projectRayDriven(Grid2D object)
	{
		Grid2D sinogram = new Grid2D(maxTIdx, maxBetaIdx);
		sinogram.setSpacing(deltaT, deltaBeta);
		
		Box b = new Box((object.getSize()[0] * object.getSpacing()[0]), (object.getSize()[1] * object.getSpacing()[1]), 2);
		Translation trans = new Translation(-(object.getSize()[0] * object.getSpacing()[0])/2, -(object.getSize()[1] * object.getSpacing()[1])/2, -1);
        b.applyTransform(trans);
        
		double beta = 0;
		for (int betaIdx = 0; betaIdx < maxBetaIdx; ++betaIdx)
		{
			double cosBeta = Math.cos(beta);
			double sinBeta = Math.sin(beta);
			PointND source = new PointND(-focalLength * sinBeta, focalLength * cosBeta, 0d);
		    PointND detector = new PointND(-maxT / 2 * cosBeta, -maxT / 2 * sinBeta, 0d);
			SimpleVector step = new SimpleVector(deltaT * cosBeta, deltaT * sinBeta, 0d);
			for(int tIdx = 0; tIdx < maxTIdx; ++tIdx)
			{
				StraightLine ray = new StraightLine(source, detector);
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
					sinogram.setAtIndex(tIdx, betaIdx, (float)sum);
											
				}
				
				detector.getAbstractVector().add(step);
			}
			beta += deltaBeta;
		}
		
		return sinogram;
	}

}
