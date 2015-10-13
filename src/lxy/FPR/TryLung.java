package lxy.FPR;

import java.util.ArrayList;

import mmorpho.*;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.Grid3D;
import edu.stanford.rsl.conrad.utils.Configuration;
import edu.stanford.rsl.conrad.utils.ImageUtil;
import edu.stanford.rsl.conrad.utils.VisualizationUtil;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.OvalRoi;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import edu.stanford.rsl.conrad.data.numeric.NumericPointwiseOperators;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.geometry.splines.BSpline;
import edu.stanford.rsl.conrad.geometry.splines.BSplineCurveInterpolation;
import edu.stanford.rsl.conrad.geometry.splines.BSplineSurfaceInterpolation;
import edu.stanford.rsl.conrad.geometry.splines.BSplineVolumeRenderer;
import edu.stanford.rsl.conrad.geometry.splines.SurfaceBSpline;


public class TryLung {
	static String inputPath="E:\\study\\master\\mt\\3rd\\Flat  panel reconstruction\\motion compensation\\Pat122\\Pat122_CTIJ_T01.nrrd";
	public static void main(String[] args) {
		Configuration.loadConfiguration();
		new ImageJ();
		ImagePlus imp=IJ.openImage(inputPath);
	//	imp.show("input image");
		Grid3D phan=ImageUtil.wrapImagePlus(imp);
		Grid3D bg = backgroundDetection(phan);
		bg.show("background");
		Grid3D contour = LungSegemntation(bg);
		contour.show("contour");
		
		
	/*	Grid2D slice = contour.getSubGrid(20);
		ContourTracer tracer = new ContourTracer(slice, 54);
		tracer.mooreNeighborhoodTracer(0).show("test");
		tracer.resample();*/
	
        /*
         * 	Grid3D contour = LungSegemntation(bg);
        contour.show("contour");
		Grid2D slice = contour.getSubGrid(54);
		ContourTracer tracer = new ContourTracer(slice);
		tracer.mooreNeighborhoodTracer(0).show("test");
		System.out.println(tracer.getCountourPoints());
		int degree = 3;
		BSplineCurveInterpolation test = new BSplineCurveInterpolation(tracer.getCountourPoints(), degree);
		BSpline spline_closed = test.applyInterpolation(1, 0);
		System.out.println(spline_closed.getControlPoints());
		VisualizationUtil.createSplinePlot(spline_closed).show();*/

		
/*	Grid3D phan2=setThres(phan,-175);
	phan2.show("after threshold");
	setZThres(phan2,13,98);
	phan2.show("after Z thres");*/
/*	
	OvalRoi eclRoi=new OvalRoi(133, 177, 298, 167);
	getROIarea(phan2,eclRoi, 13, 52);
	eclRoi=new OvalRoi(139, 190, 288, 151);
	getROIarea(phan2,eclRoi, 53, 63);
	eclRoi=new OvalRoi(144, 193, 285, 144);
	getROIarea(phan2,eclRoi, 64, 75);
	eclRoi=new OvalRoi(170, 205, 221, 126);
	getROIarea(phan2,eclRoi, 76, 84);	
	eclRoi=new OvalRoi(170, 226, 221, 94);
	getROIarea(phan2,eclRoi, 85, 98);
	phan2.show("after ROI selection");*/
	}
	
	static public Grid3D backgroundDetection(Grid3D img)
	{
		Grid3D background = new Grid3D(img);
		int thres = - 200;
		int ystop = 373;
		for(int z = 0; z < img.getSize()[2]; z ++)
		{
			for(int i = 0; i < img.getSize()[0]; i ++)
			{
				for(int j = 0; j < img.getSize()[1]; j ++)
				{
					if(img.getAtIndex(i, j, z) > thres)
						break;
					else
						background.setAtIndex(i, j, z, -200);
				}
				for(int j = img.getSize()[1] - 1; j > ystop; j --)
				{
					background.setAtIndex(i, j, z, -200);
				}
				for(int j = ystop; j >= 0; j --)
				{
									
					if(img.getAtIndex(i, j, z) > thres)
						break;
					else
						background.setAtIndex(i, j, z, -200);
				}
	
			}
		}

		return background;
	}
	
	static public Grid3D LungSegemntation(Grid3D img)
	{
		int start = 25, end = 90;
		Grid3D left = new Grid3D(img.getSize()[0], img.getSize()[1], end - start);
		//Grid3D sample = new Grid3D(img.getSize()[0], img.getSize()[1], end - start);
		ArrayList<ArrayList<PointND>> list = new ArrayList<ArrayList<PointND>>();
		/*setZThres(res, start, end);
		res.show();*/
		for(int z = start; z < end; z ++)
		{
			ImageProcessor revan = new FloatProcessor(img.getSubGrid(z).getWidth(), img.getSubGrid(z).getHeight(), img.getSubGrid(z).getBuffer()); 
			ImageProcessor bp = (ByteProcessor)revan.convertToByte(true);
			bp.autoThreshold();
			int eltype = Constants.CIRCLE,  shift = 1,  radius = 5, offset[] = Constants.OFFSET0;
			StructureElement se = new StructureElement(eltype,  shift,  radius, offset);
			MorphoProcessor mp = new MorphoProcessor(se);
			mp.open(bp);
			
			left.setSubGrid(z - start ,ImageUtil.wrapImageProcessor(bp));
			Grid2D slice = ImageUtil.wrapImageProcessor(bp);
			ContourTracer tracer = new ContourTracer(slice, (z - start) * 2);
			left.setSubGrid(z - start ,tracer.mooreNeighborhoodTracer(0));
			
			//System.out.println(z + " slice has data points: "+ tracer.getCountourPoints().size());
			//sample.setSubGrid(z - start ,tracer.resample());
			list.add(tracer.getCountourPoints());		
		}
		//sample.show("sampling");
		
		resample(list, start);
		return left;
	}
	
	static void resample(ArrayList<ArrayList<PointND>> input, int startingSlice)
	{

		ArrayList<PointND> output = new ArrayList<PointND>();

		int size = input.size();
		int numU = 20;
		int numV = size / 5;
		
		Grid3D pointsMap = new Grid3D(512, 512, numV);
		for(int i = size - 1; i >= 0; i -= 5)
		{
			ArrayList<PointND>  sliceI = input.get(i);
			int slicePointnum = sliceI.size();
			if(startingSlice + i < 65)
			{
				int step = 18;
				for(int j = 0; j < numU / 2 ; j ++)
				{
					PointND cur = sliceI.get(j * step);
					output.add(cur);
					pointsMap.setAtIndex((int)cur.getCoordinates()[0], (int)cur.getCoordinates()[1],  (size - i)/5 ,1000);
				}
				for(int j = numU - numU / 2 ; j >= 0 ; j --)
				{
					PointND cur = sliceI.get(slicePointnum - j * step - 1);
					output.add(cur);
					pointsMap.setAtIndex((int)cur.getCoordinates()[0], (int)cur.getCoordinates()[1], (size - i)/5, 1000);
				}
				
			}
			else if(startingSlice + i < 80)
			{
				int step = slicePointnum/numU - 3;
				for(int j = 0; j < numU / 2; j ++)
				{
					PointND cur = sliceI.get(j * step);
					output.add(cur);
					pointsMap.setAtIndex((int)cur.getCoordinates()[0], (int)cur.getCoordinates()[1], (size - i)/5,1000);
				}
				for(int j = numU - numU / 2 ; j >= 0 ; j --)
				{
					PointND cur = sliceI.get(slicePointnum - j * step - 1);
					output.add(cur);
					pointsMap.setAtIndex((int)cur.getCoordinates()[0], (int)cur.getCoordinates()[1], (size - i)/5, 1000);
				}
				
			}
			else
			{
				int step = slicePointnum/numU - 2;
				for(int j = 0; j < numU / 2 + 2; j ++)
				{
					PointND cur = sliceI.get(j * step);
					output.add(cur);
					pointsMap.setAtIndex((int)cur.getCoordinates()[0], (int)cur.getCoordinates()[1], (size - i)/5,1000);
				}
				for(int j = numU - numU / 2 - 2 ; j >= 0 ; j --)
				{
					PointND cur = sliceI.get(slicePointnum - j * step - 1);
					output.add(cur);
					pointsMap.setAtIndex((int)cur.getCoordinates()[0], (int)cur.getCoordinates()[1], (size - i)/5, 1000);
				}
				
			}
		}
		pointsMap.show("sampling");
		System.out.println(output);
		BSplineCurveInterpolation slice = new BSplineCurveInterpolation(output, 3);
		BSpline spline = slice.applyInterpolation(1, 0);
		for(int j = 0; j < slice.parameters.length; j++)
		{
			PointND res = spline.evaluate(slice.parameters[j]);
			System.out.println("RES at "+slice.parameters[j]+" is "+res);
			
		}
		VisualizationUtil.createSplinePlot(spline).show();

		int degreeU = 3, degreeV = 3;
		BSplineSurfaceInterpolation test = new BSplineSurfaceInterpolation(output, degreeU, degreeV, numV, numU);
		SurfaceBSpline Sbsp = test.applyInterpolation(1, 0);
		BSplineVolumeRenderer test1 = new BSplineVolumeRenderer(Sbsp);

		
	}
	
	
}
