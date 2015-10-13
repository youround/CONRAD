package lxy.ex4;

import ij.ImageJ;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import lxy.BackProjection;
import lxy.Projection;
import lxy.filtering;
import lxy.myphantoms;

public class OpenCLBackProjectionTest {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new ImageJ(); 		
		int x = 300;
		int y = 300;
		myphantoms phan = new myphantoms(x, y);
		phan.show("The Phantom");
		
		Projection projector = new Projection(Math.PI, Math.PI/180.0, 400, 1);
		Grid2D sinogram = projector.projectRayDriven(phan);
		sinogram.show("The Sinogram");
		
		Grid2D filteredSinogram = new Grid2D(sinogram);
		//Ramp filtering
		filtering Ramp = new filtering(400, 1, 0);
		Ramp.applyFiltering(filteredSinogram);
		filteredSinogram.show("Ramp Filtered Sinogram");
		
		//RamLak filtering
		Grid2D filteredSinogram2 = new Grid2D(sinogram);
		filtering RamLak = new filtering(400, 1, 1);
		RamLak.applyFiltering(filteredSinogram2);
		filteredSinogram2.show("RamLak Filtered Sinogram");
		
		// Backproject and show
		BackProjection backproj = new BackProjection(x, y, 1, 1);
		backproj.backprojectPixelDriven(filteredSinogram).show("The Reconstruction of Ramp filtered image");
		backproj.backprojectPixelDrivenCL(filteredSinogram2).show("The Reconstruction of RamLak filtered image");
		backproj.backprojectPixelDrivenCL(sinogram).show("Direct Reconstruction");

	}

}
