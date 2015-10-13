package lxy.FPR;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.NumericPointwiseOperators;
import edu.stanford.rsl.conrad.filtering.ExtremeValueTruncationFilter;
import edu.stanford.rsl.conrad.utils.ImageUtil;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class lungSegmentation {

	static public Grid2D thresholdFiltering(Grid2D img, float threshold)
	{
		ImageProcessor revan = new FloatProcessor(img.getWidth(), img.getHeight(), img.getBuffer()); 
		ImageProcessor bp = (ByteProcessor)revan.convertToByte(true);
		bp.autoThreshold();
		/*	for(int x = 0; x < img.getWidth();x++)
			{
				if(img.getAtIndex(x, y)< threshold)
					img.setAtIndex(x, y, 0);
				else
					img.setAtIndex(x, y, 1);
			}
		}*/
		
		return ImageUtil.wrapImageProcessor(bp);
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		new ImageJ();
		ImagePlus imp=IJ.openImage("E:\\study\\master\\mt\\3rd\\Flat  panel reconstruction\\Pat122\\Pat122_CTIJ_T07.nrrd");
		imp.show();
		
	    Grid2D phan=ImageUtil.wrapImagePlus(imp).getSubGrid(80);
	    phan.show("phantom");
	    lungSegmentation.thresholdFiltering(phan, 0.5f).show();
	   // lungSegmentation.thresholdFiltering(phan, 0.5f).show();

	}

}
