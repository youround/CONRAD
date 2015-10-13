package lxy.FPR;

import java.util.Collection;
import java.util.LinkedList;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import ij.IJ;
import ij.ImageJ;
import edu.stanford.rsl.conrad.data.numeric.Grid1D;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.Grid3D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.utils.Configuration;
import edu.stanford.rsl.conrad.utils.ImageUtil;
import edu.stanford.rsl.conrad.utils.RegKeys;
import edu.stanford.rsl.conrad.utils.VisualizationUtil;
import edu.stanford.rsl.jpop.FunctionOptimizer;
import edu.stanford.rsl.jpop.FunctionOptimizer.OptimizationMode;
//import edu.stanford.rsl.science.unberath.coronary.structuring.Morphology;

/**
 * Implements horizontal integration following Blondel et al. 2006,
 * allowing image-based gating and possibly respiratory motion compensation
 * @author Mathias
 *
 */
public class HorizontalIntegrationDetectorV {

	private Grid3D projections = null;
	
	private double sigma = 0.1;
	
	/** This is related to update steps of the optimizer. */
	private int numTrustDigits = 2; // or 1;
	
	private static final double min = -50; // mm
	private static final double max = +50; // mm
	
	public static void main(String[] args){
		
		Configuration.loadConfiguration();
		//String filename = "D:/Data/Coronary/breathing/cavarev_COM.tif";
		String filename = "D:/Data/Coronary/inpaintTest/dsa.tif";
			
		new ImageJ();
		
		HorizontalIntegrationDetectorV integration = new HorizontalIntegrationDetectorV(filename);
		integration.run();
		
	}
	
	
	public HorizontalIntegrationDetectorV(String filename){
		this.projections = ImageUtil.wrapImagePlus(IJ.openImage(filename));
//		this.projections.setSpacing(0.32,0.32,1);
//		this.projections.show();
//		System.out.println();
	}
	
	public double[] run(){
		Grid1D[] integrals = integrateHorizontally(projections, false);
		
		for(int i = 0; i < 2; i++){
			integrals[i].show();
		}
		
		double[] oneDshifts = optimizeFunction(integrals);
		VisualizationUtil.createPlot("Detector v-Shifts for all projections", oneDshifts).show();
		Grid3D shifted = applyShift(oneDshifts);
		shifted.show();
		
		Grid1D[] integrals2 = integrateHorizontally(shifted, false);
		for(int i = 0; i < 2; i++){
			integrals2[i].show();
		}
		
		
		return oneDshifts;
	}

	
	
	private double[] optimizeFunction(Grid1D[] lineIntegrals){
		
		double[] initVal = new double[lineIntegrals.length];
		double[] minima = new double[lineIntegrals.length];
		double[] maxima = new double[lineIntegrals.length];
		for(int i = 0; i < lineIntegrals.length; i++){
			minima[i] = min;
			maxima[i] = max;
		}
		
		OneDimDetectorMotionCompensation optimFunc = new OneDimDetectorMotionCompensation(lineIntegrals);		
		
		FunctionOptimizer fo = new FunctionOptimizer();
		fo.setConsoleOutput(true);
		fo.setDimension(lineIntegrals.length);
		fo.setOptimizationMode(OptimizationMode.Function);
		fo.setNdigit(numTrustDigits);
		fo.setInitialX(initVal);
		fo.setMaxima(maxima);
		fo.setMinima(minima);

		double[] out = fo.optimizeFunction(optimFunc);
		return out;
	}

	private Grid1D[] integrateHorizontally(Grid3D projections, boolean weight) {
		int[] s = projections.getSize();
		double[] spacing = projections.getSpacing();
		
		Grid1D[] integrals = new Grid1D[s[2]];
		
		// precompute Gaussian-weighting Kernel
		float[] gaussianKern = new float[s[1]];
		if(weight){
			for(int i = 0; i < s[1]; i++){
				gaussianKern[i] = (float)(Math.exp(-Math.pow((i-s[1]/2) / (2*sigma*s[1]) , 2)) / (Math.sqrt(2*Math.PI)*sigma*s[1]));
			}
		}else{
			for(int i = 0; i < s[1]; i++){
				gaussianKern[i] = 1;
			}
		}
				
		ExecutorService executorService = Executors.newFixedThreadPool(
				Integer.valueOf(Configuration.getGlobalConfiguration().getRegistryEntry(RegKeys.MAX_THREADS)));
		Collection<Future<?>> futures = new LinkedList<Future<?>>();
						
		for(int count = 0; count < s[2]; count++){
			final int k = count;
			futures.add(
				executorService.submit(new Runnable() {
					@Override
					public void run() {						
						System.out.println("Integrating view "+String.valueOf(k+1)+" of "+ s[2]+".");
						Grid2D grid = projections.getSubGrid(k);
						Grid1D integ = new Grid1D(s[1]);
						integ.setSpacing(spacing[1]);
						for(int j = 0; j < s[1]; j++){
							float sum = 0;
							for(int i = 0; i < s[0]; i++){
								sum += grid.getAtIndex(i, j) * spacing[0];
							}
							integ.setAtIndex(j, sum*gaussianKern[j]);
						}
						integrals[k] = integ;
					}
				})
			);
		}
		for (Future<?> future : futures){
			   try{
			       future.get();
			   }catch (InterruptedException e){
			       throw new RuntimeException(e);
			   }catch (ExecutionException e){
			       throw new RuntimeException(e);
			   }
		}		
		return integrals;
	}
	
	/**
	 * Applies the detector v-coordinate projection shift
	 * @param shifts
	 * @return
	 */
	private Grid3D applyShift(double[] shifts){
		int[] s = projections.getSize();
		Grid3D shifted = new Grid3D(s[0],s[1],s[2]);
		double[] spacing = projections.getSpacing();
		shifted.setSpacing(spacing);
		
		ExecutorService executorService = Executors.newFixedThreadPool(
				Integer.valueOf(Configuration.getGlobalConfiguration().getRegistryEntry(RegKeys.MAX_THREADS)));
		Collection<Future<?>> futures = new LinkedList<Future<?>>();
						
		for(int count = 0; count < s[2]; count++){
			final int k = count;
			futures.add(
				executorService.submit(new Runnable() {
					@Override
					public void run() {						
						System.out.println("Applying shift to view "+String.valueOf(k+1)+" of "+ s[2]+".");
						Grid2D grid = projections.getSubGrid(k);
						Grid2D shi = new Grid2D(s[0],s[1]);
						shi.setSpacing(spacing[0], spacing[1]);
						shi.setOrigin(0, -shifts[k]);
						
						for(int i = 0; i < s[0]; i++){
							for(int j = 0; j < s[1]; j++){
								float val = grid.getAtIndex(i, j);
								double[] world = grid.indexToPhysical(i, j);
								double[] ind = shi.physicalToIndex(world[0], world[1]);
								if(ind[0] < 0 || ind[0] > s[0]-1 || ind[1] < 0 || ind[1] >= s[1]-1){
								}else{
									InterpolationOperators.addInterpolateLinear(shi, ind[0], ind[1], val);
								}
							}
							
						}
						shifted.setSubGrid(k, shi);
					}
				})
			);
		}
		for (Future<?> future : futures){
			   try{
			       future.get();
			   }catch (InterruptedException e){
			       throw new RuntimeException(e);
			   }catch (ExecutionException e){
			       throw new RuntimeException(e);
			   }
		}
		return shifted;
	}
	
	
}
