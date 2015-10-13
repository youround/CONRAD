package lxy.ex4;

import ij.ImageJ;
import edu.stanford.rsl.conrad.data.numeric.NumericPointwiseOperators;
import edu.stanford.rsl.conrad.data.numeric.opencl.OpenCLGrid2D;
import lxy.myphantoms;
import edu.stanford.rsl.tutorial.phantoms.Phantom;
import edu.stanford.rsl.tutorial.phantoms.UniformCircleGrid2D;

public class OpenCLGridTest {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		new ImageJ(); 		
		int x = 512;
		int y = 512;
		int N = 999;
		myphantoms phan = new myphantoms(x, y);
		phan.show("My Phantom");

		
		// OpenCL
		// Add the phantom to itself  for 1000 times
		OpenCLGrid2D gridCL = new OpenCLGrid2D(phan);
		OpenCLGrid2D resultCL = new OpenCLGrid2D(phan);
		long starttime= System.nanoTime();
		for (int i = 0; i < N; ++i)
		{
			NumericPointwiseOperators.addBy(resultCL, gridCL);
		}
		
		long endtime= System.nanoTime();
		long timecost= endtime - starttime;
		System.out.println("TIme cost on openCL: " + timecost);
		
		
		// CPU
		myphantoms result =  new myphantoms(x, y);
		starttime= System.nanoTime();
		for (int i = 0; i < N; ++i)
		{
			NumericPointwiseOperators.addBy(result, phan);
		}
        endtime= System.nanoTime();
		timecost= endtime - starttime;
		System.out.println("TIme cost on CPU: " + timecost);


	}

}
