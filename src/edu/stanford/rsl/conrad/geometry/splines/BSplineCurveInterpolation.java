package edu.stanford.rsl.conrad.geometry.splines;

import java.util.ArrayList;
import java.util.Arrays;

import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.numerics.Solvers;
import edu.stanford.rsl.conrad.utils.VisualizationUtil;

/**
 * Class to fit a set of data points with a B-spline curve using global interpolation.
 * 
 * @author XinyunLi
 *
 */
public class BSplineCurveInterpolation {
	
	private ArrayList<PointND> dataPoints;
	
	protected int dimension;
	protected int degree;//degree
	protected int n; //number of data points
	protected int m; //number of control points
	public double[] parameters;
	
	public BSplineCurveInterpolation(ArrayList<PointND> dataPoints, int degree)
	{
		this.dataPoints = dataPoints;
		this.degree = degree;
		this.dimension =  dataPoints.get(0).getDimension();
		this.n = dataPoints.size();
		this.m = n;
		
	}
	/**
	 * 
	 * @param parametrization aa
	 * @param endpointCondition bb
	 * @return
	 */
	public BSpline applyInterpolation(int parametrization, int endpointCondition )
	{
		//compute the parameter value
		if(endpointCondition == 2) m = n + degree;
		double[] parameters = buildParameterVector(parametrization, endpointCondition);
		//System.out.println(Arrays.toString(parameters));
		//compute the knot vector
		double[] knots = buildKnotVector(parameters, endpointCondition);
		System.out.println("Knots" + Arrays.toString(knots));			
		// set up the linear equation system
		SimpleMatrix A = buildCoefficientMatrix(knots, parameters, endpointCondition);
		//System.out.println(A);
			
		SimpleVector[] b = buildDataVector(knots, endpointCondition);
		
		ArrayList<PointND> controlPoints = solveLinearEquations(A, b, endpointCondition);
		//System.out.println(controlPoints);
		
		return new BSpline(controlPoints, knots);		

	}
	
	public BSpline applyInterpolation(double[] parameters, double[] knots,  int endpointCondition)
	{
		SimpleMatrix A = buildCoefficientMatrix(knots, parameters, endpointCondition);
		SimpleVector[] b = buildDataVector(knots, endpointCondition);	
		ArrayList<PointND> controlPoints = solveLinearEquations(A, b, endpointCondition);
		return new BSpline(controlPoints, degree, knots);		
		
	}
	public BSpline applyInterpolation(int parametrization, int endpointCondition, ArrayList<PointND> endDerivatives)
	{
		double[] parameters = buildParameterVector(parametrization, endpointCondition);
		System.out.println(Arrays.toString(parameters));
		//compute the knot vector
		double[] knots = buildKnotVector(parameters, endpointCondition);
		System.out.println(Arrays.toString(knots));			
		// set up the linear equation system
		SimpleMatrix A = buildCoefficientMatrix(knots, parameters, endpointCondition);
		System.out.println(A);
    	SimpleVector[] b = buildDataVector(knots, endpointCondition);
		SimpleMatrix X = new SimpleMatrix(n, dimension);	
		for (int i = 0; i < dimension; i ++)
		{
            SimpleVector xi =  Solvers.solveLinearSysytemOfEquations(A,  b[i]);
            X.setColValue(i, xi);
		}
		System.out.println(X);
		ArrayList<PointND> controlPoints = new ArrayList<PointND>();
		for (int i = 0; i < n; i ++)
		{    
			controlPoints.add(new PointND(X.getRow(i)));
		}
		
		/*SimpleMatrix X = new SimpleMatrix(n + 2, dimension);	
		for (int i = 0; i < dimension; i ++)
		{
            SimpleVector xi =  Solvers.solveLinearSysytemOfEquations(A,  b[i]);
            X.setColValue(i, xi);
		}
		System.out.println(X);
		ArrayList<PointND> controlPoints = new ArrayList<PointND>();
		for (int i = 0; i < n + 2; i ++)
		{    
			controlPoints.add(new PointND(X.getRow(i)));
		}*/
		
		return new BSpline(controlPoints, degree, knots);	
	}
	
	// parametrization options: 0: uniform; 1:chordal; 2:centripetal
	protected double[] buildParameterVector(int parametrization, int endpointCondition)
	{
		
		parameters = new double[n];
		//closed curve
		if(endpointCondition == 2)
		{
			//very special case, the parameter corresponds to the first data point are both 0 and 1
			//So actually we have n+1 parameter, but we only take the first n for solving the linear equations
			if (parametrization == 0)
			{
				for(int i = 1; i < n - 1; i++)
					parameters[i] = i / n;
			}
			else if (parametrization == 1)
			{
				for(int i = 1; i < n; i++)
				{
					double dis = dataPoints.get(i).euclideanDistance(dataPoints.get(i - 1));
					parameters[i] = parameters[i - 1] + dis;
				}
				double total = parameters[n - 1] + dataPoints.get(0).euclideanDistance(dataPoints.get(n - 1));
				for(int i = 1; i < n; i++)
					parameters[i] /= total;
			}
			else if (parametrization == 2)
			{
				for(int i = 1; i < n; i++)
				{
					double dis = dataPoints.get(i).euclideanDistance(dataPoints.get(i - 1));
					parameters[i] = parameters[i - 1] + Math.sqrt(dis);
				}
				double total = parameters[n - 1] + dataPoints.get(0).euclideanDistance(dataPoints.get(n - 1));
				for(int i = 1; i < n; i++)
					parameters[i] /= total;			
			}
		}
		else
		{
			//clamped, multiplicity n+1 at endpoints
			//compute parameters uk
			
			if (parametrization == 0)
			{
				for(int i = 1; i < n - 1; i++)
					parameters[i] = i / (n - 1);
			}
			else if (parametrization == 1)
			{
				for(int i = 1; i < n; i++)
				{
					double dis = dataPoints.get(i).euclideanDistance(dataPoints.get(i - 1));
					parameters[i] = parameters[i - 1] + dis;
				}
				for(int i = 1; i < n - 1; i++)
					parameters[i] /= parameters[n - 1];
			}
			else if (parametrization == 2)
			{
				for(int i = 1; i < n; i++)
				{
					double dis = dataPoints.get(i).euclideanDistance(dataPoints.get(i - 1));
					parameters[i] = parameters[i - 1] + Math.sqrt(dis);
				}
				for(int i = 1; i < n - 1; i++)
					parameters[i] /= parameters[n - 1];			
			}
			parameters[n - 1] = 1d; 
		}
		return parameters;
	}
	
	// endpoints conditions: 0: clamped; 1: open; 2:closed
	protected double[] buildKnotVector(double[] parameters, int endpointCondition)
	{
		//with endderivative
/*		if (endpointCondition == 3)
		{
			m = n + 2;
			int k = m + degree + 1; //number of knots
			double[] knots = new double[k];
			// compute the knot vector
			for(int i = 0; i <= degree; i++)
			{
				knots[i] = 0;
				knots[k-i-1] = 1;
			}
			
			for(int j = 0; j <= m - degree - 1; j++)
			{
				double sum = 0;
				for(int i = j; i < j + degree; i++)
					sum += parameters[i];
				knots[j + degree + 1] =  sum / degree;
			}
			return knots;
			
		}
		else if(endpointCondition == 1)
		{
			int k = m + degree + 1; //number of knots
			double[] knots = new double[k];
			knots[m] = 1;
			for(int j = 1; j <= m - degree - 1; j++)
			{
				double sum = 0;
				for(int i = j; i < j + degree; i++)
					sum += parameters[i];
				knots[j + degree] =  sum / degree;
			}
			for(int i = 0; i < degree; i++)
			{
				knots[degree - i - 1] = knots[degree - i] - (knots[m - i] - knots[m - i - 1]);
				knots[m + i + 1] = knots[m + i] + (knots[degree + i + 1] - knots[degree + i]);
						
			}
			return knots;
		}
		else if(endpointCondition == 2)
		{
			//uniform knots
			m = n + degree;
			int k = m + degree + 1; //number of knots
			double[] knots = new double[k];
			for(int j = 0; j < k; j++)
				knots[j] = ((double)(j - degree)) / (m - degree);
			return knots;
			
			
		}
		else 
		{
			int k = m + degree + 1; //number of knots
			double[] knots = new double[k];
			// compute the knot vector
			for(int i = 0; i <= degree; i++)
			{
				knots[i] = 0;
				knots[k-i-1] = 1;
			}
			
			for(int j = 1; j <= m - degree - 1; j++)
			{
				double sum = 0;
				for(int i = j; i < j + degree; i++)
					sum += parameters[i];
				knots[j + degree] =  sum / degree;
			}
			return knots;
			
		}*/
		
		int num = parameters.length;
		if(endpointCondition == 1)
		{
			int k = num + degree + 1; //number of knots
			double[] uKnots = new double[k];
			uKnots[num] = 1;
			for(int j = 1; j <= num - degree - 1; j++)
			{
				double sum = 0;
				for(int i = j; i < j + degree; i++)
					sum += parameters[i];
				uKnots[j + degree] =  sum / degree;
			}
			for(int i = 0; i < degree; i++)
			{
				uKnots[degree - i - 1] = uKnots[degree - i] - (uKnots[num - i] - uKnots[num - i - 1]);
				uKnots[num + i + 1] = uKnots[num + i] + (uKnots[degree + i + 1] - uKnots[degree + i]);
						
			}
			return uKnots;
		}
		else if(endpointCondition == 2)
		{
			//uniform knots
			int k = num + degree * 2 + 1; //number of knots
			double[] knots = new double[k];
			for(int j = 0; j < k; j++)
				knots[j] = ((double)(j - degree)) / num;
			return knots;
			
			
		}
		else 
		{
			int k = num + degree + 1; //number of knots
			double[] uKnots = new double[k];
			// compute the knot vector
			for(int i = 0; i <= degree; i++)
			{
				uKnots[i] = 0;
				uKnots[k-i-1] = 1;
			}
			
			for(int j = 1; j <= num - degree - 1; j++)
			{
				double sum = 0;
				for(int i = j; i < j + degree; i++)
					sum += parameters[i];
				uKnots[j + degree] =  sum / degree;
			}
			return uKnots;
			
		}
		

	}

	protected SimpleMatrix buildCoefficientMatrix(double[] knots, double[] parameters, int endpointCondition )
	{
		//default initialization
		if (endpointCondition == 3)
		{
			double[][] A = new double[n + 2][n + 2];
			A[0][0] = 1;
			A[1][0] = -1; A[1][1] = 1;
			for (int i = 1; i < n - 1; i ++)
			{
				int index = Arrays.binarySearch(knots, parameters[i]);
				if (index < 0) index = -1 * index - 2; 
				if (index < degree) index = degree;
				if (index >= n) index = n-1;
				double[] Ai = BasisFuns(knots, parameters[i], index);
				System.arraycopy(Ai, 0, A[i + 1], index - degree, Ai.length );
			}
			A[n][n] = -1; A[n][n + 1] = 1;
			A[n + 1][n + 1] = 1;
			return new SimpleMatrix(A);
		}
		else if(endpointCondition == 2)
		{
			double[][] A = new double[n][n];
			for (int i = 0; i < n; i ++)
			{
				int index = Arrays.binarySearch(knots, parameters[i]);
				if (index < 0) index = -1 * index - 2; 
				if (index < degree) index = degree;
				if (index >= m) index = m-1;
				double[] Ai = BasisFuns(knots, parameters[i], index);
				if(index >= n)
				{
					System.arraycopy(Ai, 0, A[i], index - degree, degree + n -index);
					System.arraycopy(Ai, degree + n -index, A[i], 0, index - n + 1);
				}
				else
					System.arraycopy(Ai, 0, A[i], index - degree, Ai.length );
			}
			return new SimpleMatrix(A);
		}
			
		else
		{
			double[][] A = new double[n][n];
			for (int i = 0; i < n; i ++)
			{
				int index = Arrays.binarySearch(knots, parameters[i]);
				if (index < 0) index = -1 * index - 2; 
				if (index < degree) index = degree;
				if (index >= n) index = n-1;
				double[] Ai = BasisFuns(knots, parameters[i], index);
				System.arraycopy(Ai, 0, A[i], index - degree, Ai.length );
			}

			return new SimpleMatrix(A);
		}
		
	}
	
	protected SimpleVector[] buildDataVector(double[] knots, int endpointCondition)
	{
		SimpleVector[] b = new SimpleVector[dimension];

		double[][] dataPointsArray = new double[dimension][n];
 		for (int i = 0; i < dimension; i++)
		{
			for(int j = 0; j < n; j++)
				dataPointsArray[i][j] = dataPoints.get(j).get(i);
		}
		for (int i =0; i < dimension; i++)
		{
			b[i] = new SimpleVector(dataPointsArray[i]);
		}		
		return b;
		
	}
	
/*	protected SimpleVector[] buildDataVector(double[] knots, int endpointCondition, ArrayList<PointND> endDerivatives)
	{
		SimpleVector[] b = new SimpleVector[dimension];
		if(endpointCondition == 1)
		{
			double[][] dataPointsArray = new double[dimension][n + 2];
			for (int i =0; i < dimension; i++)
			{
				dataPointsArray[i][0] = dataPoints.get(0).get(i);
				dataPointsArray[i][n + 1] = dataPoints.get(n - 1).get(i);
				dataPointsArray[i][1] = knots[degree + 1] / degree * endDerivatives.get(0).get(i);
				dataPointsArray[i][n] = (1 - knots[n + 1]) / degree * endDerivatives.get(1).get(i);
			}
			for (int i =0; i < dimension; i++)
			{
				for(int j = 1; j < n - 1; j++)
					dataPointsArray[i][j + 1] = dataPoints.get(j).get(i);
			}
			for (int i =0; i < dimension; i++)
			{
				b[i] = new SimpleVector(dataPointsArray[i]);
			}
		}
		
		
		return b;
	}*/
	
	protected ArrayList<PointND> solveLinearEquations(SimpleMatrix A, SimpleVector[] b, int endpointCondition)
	{
		SimpleMatrix X = new SimpleMatrix(n, dimension);	
		for (int i = 0; i < dimension; i ++)
		{
            SimpleVector xi =  Solvers.solveLinearSysytemOfEquations(A,  b[i]);
            X.setColValue(i, xi);
		}

		ArrayList<PointND> controlPoints = new ArrayList<PointND>();
		for (int i = 0; i < n; i ++)
		{    
			controlPoints.add(new PointND(X.getRow(i)));
		}
		if(endpointCondition == 2)
		{
			for (int i = 0; i < degree; i ++)
			{    
				controlPoints.add(new PointND(X.getRow(i)));
			}
		}
		
		return controlPoints;		
	}
	
	protected double[] BasisFuns(double[] knots, double internalCoordinate, int index) 
	{
		//Nurbs Book Algorithms A2.2 
		//Compute the nonvanishing basis functions
		double[] N = new double[degree + 1];
		double[] left = new double[degree + 1];
		double[] right = new double[degree + 1];
		N[0] = 1.0d;
		for (int j = 1; j <= degree; j++)
		{
			left[j] = internalCoordinate - knots[index + 1 - j];
			right[j] = knots[index + j] - internalCoordinate;
			double saved = 0.0d;
			for (int r = 0; r < j; r++)
			{
				double temp = N[r] / (right[r+1] + left[j - r]);
				N[r] = saved + right[r + 1] * temp;
				saved  = left[j - r]  * temp;
			}
			N[j] = saved;
		}
		return N;
	}

	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		ArrayList<PointND> list = new ArrayList<PointND>();
		list.add(new PointND (-2, 5));
		list.add(new PointND (-5, 0));
		list.add(new PointND (-1, -3));
		list.add(new PointND (2, -3));
		list.add(new PointND (5, 2));
		list.add(new PointND (1, 3));

		
		ArrayList<PointND> endder = new ArrayList<PointND>();
		endder.add(new PointND (-5, 5));
		endder.add(new PointND (0, 10));
		
		int degree = 3;
		BSplineCurveInterpolation test = new BSplineCurveInterpolation(list, degree);
		BSpline spline_clampped = test.applyInterpolation(1, 0);
		BSpline spline_open = test.applyInterpolation(1, 1);
		BSpline spline_closed = test.applyInterpolation(1, 2);
		//BSpline spline = test.applyInterpolation(1, 1, endder);
		VisualizationUtil.createSplinePlot(spline_clampped).show();
		VisualizationUtil.createSplinePlot(spline_open).show();
		VisualizationUtil.createSplinePlot(spline_closed).show();
		

	}

}
