package lxy;

import java.io.IOException;
import java.nio.FloatBuffer;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLCommandQueue;
import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLDevice;
import com.jogamp.opencl.CLImage2d;
import com.jogamp.opencl.CLImageFormat;
import com.jogamp.opencl.CLKernel;
import com.jogamp.opencl.CLProgram;
import com.jogamp.opencl.CLImageFormat.ChannelOrder;
import com.jogamp.opencl.CLImageFormat.ChannelType;
import com.jogamp.opencl.CLMemory.Mem;

import edu.stanford.rsl.conrad.data.numeric.Grid1D;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.data.numeric.NumericPointwiseOperators;
import edu.stanford.rsl.conrad.numerics.SimpleOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.opencl.OpenCLUtil;

public class BackProjection {
	
	int imgSzX, imgSzY;       // image dimensions [GU]
	float pxSzX, pxSzY;       // pixel size [mm]
	
	public BackProjection(int imageSizeX, int imageSizeY, float pixelSizeX, float pixelSizeY)
	{
		this.imgSzX = imageSizeX;
		this.imgSzY = imageSizeY;
		this.pxSzX = pixelSizeX;
		this.pxSzY = pixelSizeY;
	}
	
	
	
	public Grid2D backprojectPixelDriven(Grid2D sino) 
	{
		Grid2D rec = new Grid2D(imgSzX, imgSzY);
		rec.setSpacing(pxSzX, pxSzY);
		rec.setOrigin(-(rec.getSize()[0]*rec.getSpacing()[0])/2, -(rec.getSize()[1]*rec.getSpacing()[1])/2);
		
		int maxThetaIdx = sino.getSize()[1];
		int maxSIdx = sino.getSize()[0];
		double deltaS = sino.getSpacing()[0];
		double deltaTheta = sino.getSpacing()[1];
		
		for(int thetaIdx = 0; thetaIdx < maxThetaIdx; ++thetaIdx)
		{
			double theta = thetaIdx * deltaTheta;
			double sinTheta = Math.sin(theta);
			double cosTheta = Math.cos(theta);
			Grid1D subgridTheta = sino.getSubGrid(thetaIdx);
			
			for(int i = 0; i < rec.getSize()[0]; ++i)
				for(int j = 0; j < rec.getSize()[1]; ++j)
				{
					double pos[] = rec.indexToPhysical(i, j);
					double s = pos[0]*cosTheta + pos[1]*sinTheta;
					double sIdx =  (s + deltaS * (maxSIdx - 1) / 2) / deltaS;
							//subgridTheta.physicalToIndex(s);
					if(sIdx > maxSIdx - 1 || sIdx < 0)
						continue;
					float val = InterpolationOperators.interpolateLinear(subgridTheta, sIdx);
					rec.addAtIndex(i, j, (float)(val / maxThetaIdx * Math.PI));
				}
		}
//		NumericPointwiseOperators.divideBy(rec, (float) (maxThetaIdx / Math.PI));
		return rec;
	}
	
	public Grid2D backprojectPixelDrivenCL(Grid2D sino) 
	{
		int maxThetaIdx = sino.getSize()[1];
		int maxSIdx = sino.getSize()[0];
		double deltaS = sino.getSpacing()[0];
		double deltaTheta = sino.getSpacing()[1];
		
		
		CLContext context = OpenCLUtil.createContext();
		CLDevice device = context.getMaxFlopsDevice();
		
		// load sources, create and build program
		CLProgram program = null;
		try {
			program = context.createProgram(this.getClass().getResourceAsStream("backprojectionPixelDrivenCL.cl"))
					.build();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.exit(-1);
		}
		
		int sinoSize = maxThetaIdx*maxSIdx;
		// Length of arrays to process
		int localWorkSize = Math.min(device.getMaxWorkGroupSize(), 8); // Local work size dimensions
		int globalWorkSizeS = OpenCLUtil.roundUp(localWorkSize, imgSzX); // rounded up to the nearest multiple of localWorkSize
		int globalWorkSizeTheta = OpenCLUtil.roundUp(localWorkSize, imgSzY); // rounded up to the nearest multiple of localWorkSize


		// create image from input grid
		CLBuffer<FloatBuffer> inputBuffer = context.createFloatBuffer(sinoSize, Mem.READ_ONLY);
		for (int i=0;i<sinoSize;++i){
				inputBuffer.getBuffer().put(sino.getBuffer()[i]);
		}
		inputBuffer.getBuffer().rewind();

		CLImageFormat format = new CLImageFormat(ChannelOrder.INTENSITY, ChannelType.FLOAT);
        CLImage2d<FloatBuffer> inputImage = context.createImage2d(
						inputBuffer.getBuffer(), maxSIdx, maxThetaIdx, format);
		inputBuffer.release();		

		// create memory for output image
		CLBuffer<FloatBuffer> outputBuffer = context.createFloatBuffer(imgSzX*imgSzY, Mem.WRITE_ONLY);
		
		// copy params
		CLKernel kernel = program.createCLKernel("backprojectPixelDriven2DCL");
		kernel.putArg(inputImage).putArg(outputBuffer)
		.putArg(imgSzX).putArg(imgSzY)
		.putArg(pxSzX).putArg(pxSzY)
		.putArg((float)maxSIdx).putArg((float)deltaS)
		.putArg((float)maxThetaIdx).putArg((float)deltaTheta); 

		// createCommandQueue
		CLCommandQueue queue = device.createCommandQueue();
		queue
		.putWriteImage(inputImage, true)
		.finish()
		.put2DRangeKernel(kernel, 0, 0, globalWorkSizeS, globalWorkSizeTheta,
				localWorkSize, localWorkSize).putBarrier()
		.putReadBuffer(outputBuffer, true)
		.finish();
		
	// write grid back to grid2D
		Grid2D img = new Grid2D(this.imgSzX, this.imgSzY);
		img.setSpacing(pxSzX, pxSzY);
		outputBuffer.getBuffer().rewind();
		for (int i = 0; i < imgSzX*imgSzY; ++i) {
			img.getBuffer()[i] = outputBuffer.getBuffer().get();
		}
		
		context.release();

		return img;


	}

}
