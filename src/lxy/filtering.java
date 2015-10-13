package lxy;

import edu.stanford.rsl.conrad.data.numeric.Grid1D;
import edu.stanford.rsl.conrad.data.numeric.Grid1DComplex;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.NumericPointwiseOperators;

public class filtering extends Grid1DComplex {
	
	public filtering(final int size, double deltaS, int type) 
	{
		super(size);
		if (type == 0)
			buildRamp(deltaS);			
		else if (type == 1)	
			buildRamLak(deltaS);	
	}
	
	public void buildRamp(double deltaS) 
	{

		final int paddedSize = getSize()[0];
		float deltaF = 1.f/(float)(deltaS*paddedSize);
		float F = deltaF;
		for (int i = 1; i < paddedSize/2; ++i)
		{
			setAtIndex(i, F);
			setAtIndex(paddedSize - i - 1, F);
			F += deltaF;
		}
		
	}
	
	public void buildRamLak(double deltaS) 
	{

		final int paddedSize = getSize()[0];
		final float odd = -1.f / ((float) (Math.PI * Math.PI * deltaS));
		setAtIndex(0, (float) (0.25f / (deltaS)));
		for (int i = 1; i < paddedSize/2; ++i)
		{
			if (1 == (i % 2))
			{
				setAtIndex(i, odd / (i * i));
				setAtIndex(paddedSize - i, odd / (i * i));
			}
		}
		transformForward();

	}

	
	public void applyFiltering(Grid2D sino)
	{	
		for (int theta = 0; theta < sino.getSize()[1]; ++theta) 
		{
			Grid1DComplex subGrid = new Grid1DComplex(sino.getSubGrid(theta));
			subGrid.transformForward();			
			for (int idx = 0; idx < subGrid.getSize()[0]; ++idx) {
				subGrid.multiplyAtIndex(idx, getRealAtIndex(idx),
						getImagAtIndex(idx));
			}
			subGrid.transformInverse();
			Grid1D filteredSinoSub = subGrid.getRealSubGrid(0, sino.getSize()[0]);
			NumericPointwiseOperators.copy(sino.getSubGrid(theta), filteredSinoSub);
		}	
	}
	
	
		
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
