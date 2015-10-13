package lxy;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.tutorial.parallel.ParallelProjector2D;
import ij.ImageJ;

public class parallelReconstruction {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new ImageJ(); 		
		int x = 500;
		int y = 500;
		myphantoms  phan = new myphantoms(x, y);
		phan.show("The Phantom");
		
		int detectorSize = 600;
		Projection projector = new Projection(Math.PI, Math.PI/180.0, detectorSize, 1);
		Grid2D sinogram = projector.projectRayDriven(phan);
		sinogram.show("The Sinogram");
		
		Grid2D filteredSinogram = new Grid2D(sinogram);
		//Ramp filtering
		filtering Ramp = new filtering(detectorSize, 1, 0);
		Ramp.applyFiltering(filteredSinogram);
		filteredSinogram.show("Ramp Filtered Sinogram");
		
		//RamLak filtering
		Grid2D filteredSinogram2 = new Grid2D(sinogram);
		filtering RamLak = new filtering(detectorSize, 1, 1);
		RamLak.applyFiltering(filteredSinogram2);
		filteredSinogram2.show("RamLak Filtered Sinogram");
		
		// Backproject and show
		BackProjection backproj = new BackProjection(x, y, 1, 1);
		backproj.backprojectPixelDriven(filteredSinogram).show("The Reconstruction of Ramp filtered image");
		backproj.backprojectPixelDriven(filteredSinogram2).show("The Reconstruction of RamLak filtered image");
		backproj.backprojectPixelDriven(sinogram).show("Direct Reconstruction");

	}

}
