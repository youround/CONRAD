package lxy.ex3;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.tutorial.filters.RamLakKernel;
import ij.ImageJ;
import lxy.BackProjection;
import lxy.filtering;
import lxy.myphantoms;

public class fanBeamReconstruction {
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new ImageJ(); 		
		int x = 200;
		int y = 200;
		myphantoms phan = new myphantoms(x, y);
		phan.show("The Phantom");
		
		double gammaM = 20*Math.PI/180, 
				maxT = 500, 
				deltaT = 1.0, 
				focalLength = (maxT/2.0-0.5)*deltaT/Math.tan(gammaM),
				//maxBeta = 360*Math.PI/180,
				//Short scan
				maxBeta = Math.PI + 2 * gammaM,
				deltaBeta = maxBeta / 132;
		
		fanBeamProjection projector = new fanBeamProjection(
				focalLength, maxBeta, deltaBeta, maxT, deltaT);
		Grid2D fanogram = projector.projectRayDriven(phan);
		fanogram.show("The Fanogram");
		
	    //Rebinning Method - short scan
		reBinning rb = new reBinning(Math.PI, Math.PI/180.0, 400, 1);
		Grid2D sinogram = rb.applyReBinning(fanogram, focalLength);	
		sinogram.show("The Sinogram");
		
		Grid2D filteredSinogram = new Grid2D(sinogram);
		filtering RamLak = new filtering(400, 1, 1);
		RamLak.applyFiltering(filteredSinogram);
		filteredSinogram.show("RamLak Filtered Sinogram");
		
		BackProjection backproj = new BackProjection(200, 200, 1, 1);
		backproj.backprojectPixelDriven(filteredSinogram).show("The Reconstruction of Ramp filtered image");
	
        //Filtering Method
		
		filtering ramLak = new filtering((int) (maxT / deltaT), deltaT, 1);

		
		
	}

}
