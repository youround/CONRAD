package lxy.ex4;

import ij.ImageJ;

import java.io.IOException;
import java.nio.FloatBuffer;

import lxy.myphantoms;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLCommandQueue;
import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLDevice;
import com.jogamp.opencl.CLKernel;
import com.jogamp.opencl.CLProgram;
import com.jogamp.opencl.CLMemory.Mem;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.opencl.OpenCLGrid2D;
import edu.stanford.rsl.conrad.data.numeric.opencl.OpenCLGridInterface;
import edu.stanford.rsl.conrad.opencl.OpenCLUtil;
import edu.stanford.rsl.tutorial.phantoms.DotsGrid2D;

public class OpenCLKernelTest {
	
	public static Grid2D addCL(OpenCLGrid2D gridA, OpenCLGrid2D gridB)
	{
		//OpenCLGridInterface clGrid = (OpenCLGridInterface)grid;
		int imgSizeX = gridA.getWidth();
		int imgSizeY = gridA.getHeight();
				
		OpenCLGridInterface clGridA = (OpenCLGridInterface)gridA;
		OpenCLGridInterface clGridB = (OpenCLGridInterface)gridB;

		clGridA.getDelegate().prepareForDeviceOperation();
		clGridB.getDelegate().prepareForDeviceOperation();

		// TODO check if both live on the same device.
		CLDevice device = clGridA.getDelegate().getCLDevice(); 
		
		//Input and Output Buffer
		CLBuffer<FloatBuffer> clmemA = clGridA.getDelegate().getCLBuffer();
		CLBuffer<FloatBuffer> clmemB = clGridB.getDelegate().getCLBuffer();
		
		//Context
		CLContext context = device.getContext();
		//Program
		CLProgram program = null;		
		try {
			program = context.createProgram(OpenCLKernelTest.class.getResourceAsStream("add.cl")).build();
		} catch (IOException e) {
			e.printStackTrace();
			program = null;
		}
		

		//workgroup (local) size
		int localSize = Math.min(device.getMaxWorkGroupSize(), 16);
		int elementCount = clmemA.getCLCapacity(); 
		int globalSize = OpenCLUtil.roundUp(localSize, elementCount);
	
		//output
		CLBuffer<FloatBuffer> resultBuffer = null;
		resultBuffer = context.createFloatBuffer(globalSize, Mem.READ_ONLY);
		
		//Kernel
		CLKernel kernel = program.createCLKernel("add");
		kernel.putArg(clmemA).putArg(clmemB).putArg(resultBuffer).putArg(elementCount);
		
		//CommandQueue
		CLCommandQueue queue = device.createCommandQueue();
		queue.putWriteBuffer(clmemA, false);
		queue.putWriteBuffer(clmemB, false);
		queue.put1DRangeKernel(kernel, 0, globalSize, localSize);
		queue.putReadBuffer(resultBuffer, true);
			
		Grid2D img = new Grid2D(imgSizeX, imgSizeY);
		//img.setSpacing(pxSzXMM, pxSzYMM);
		resultBuffer.getBuffer().rewind();
		for (int i = 0; i < imgSizeX*imgSizeY; ++i) {
				img.getBuffer()[i] = resultBuffer.getBuffer().get();
				
		}
		
		queue.release();
		resultBuffer.release();
		kernel.release();
		program.release();
		context.release();

		return img;


		
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		int x = 512;
		int y = 512;
		myphantoms phan1 = new myphantoms(x, y);
		phan1.show("Phantom1");	
		DotsGrid2D phan2 = new DotsGrid2D(x, y);
		phan2.show("Phantom2");
		OpenCLGrid2D gridCL1 = new OpenCLGrid2D(phan1);
		OpenCLGrid2D gridCL2 = new OpenCLGrid2D(phan2);
		Grid2D result = OpenCLKernelTest.addCL(gridCL1, gridCL2);
		new ImageJ(); 	
		result.show("Addition");
		


	}

}
