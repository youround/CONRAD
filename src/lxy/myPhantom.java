package lxy;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.numerics.SimpleMatrix;
import edu.stanford.rsl.conrad.numerics.SimpleVector;

public class myPhantom extends Grid2D {
	
	private String title;
	final private SimpleMatrix Ellipses; 
	
	public myPhantom(int x, int y) {
		// TODO Auto-generated constructor stub
		super(new float[x*y], x, y);
		this.setSpacing(1,1);
		this.setOrigin(x/2.d,y/2.d);
		this.title = "My";
		Ellipses = CreateEllipses();
		CreatePhantom();

	}
	
	private SimpleMatrix CreateEllipses()
	{
		SimpleMatrix M1 = new SimpleMatrix(10,6);	
		M1.setRowValue(0, new SimpleVector(new double[] {2.0, 0.69, 0.92, 0, 0, 0}));
		M1.setRowValue(1, new SimpleVector(new double[] {-0.98, 0.6624, 0.8740, 0, -0.0184, 0}));
		M1.setRowValue(2, new SimpleVector(new double[] {-0.02, 0.1100, 0.3100, 0.22, 0.0, -18.0}));
		M1.setRowValue(3, new SimpleVector(new double[] {-0.02, 0.1600, 0.4100, -0.22, 0.0, 18.0}));
		M1.setRowValue(4, new SimpleVector(new double[] {0.01, 0.2100, 0.2500, 0, 0.35, 0}));
		M1.setRowValue(5, new SimpleVector(new double[] {0.01, 0.0460, 0.0460, 0, 0.1, 0}));
		M1.setRowValue(6, new SimpleVector(new double[] {0.01, 0.0460, 0.0460, 0, -0.1, 0}));
		M1.setRowValue(7, new SimpleVector(new double[] {0.01, 0.0460, 0.0230, -0.08, -0.605, 0}));
		M1.setRowValue(8, new SimpleVector(new double[] {0.01, 0.0230, 0.0230, 0, -0.606, 0}));
		M1.setRowValue(9, new SimpleVector(new double[] {0.01, 0.0230, 0.0460,  0.06, -0.605, 0}));
		return M1;
	}
	
	private void CreatePhantom()
	{
		double sizeX = (double)super.getSize()[0];
		double sizeY = (double)super.getSize()[1];
		
		for (int i=0; i < super.getSize()[0]; ++i)
		{
			double x = ((double)i-(sizeX-1)/2.0) / ((sizeX-1)/2.0);
			for (int j=0; j < super.getSize()[1]; ++j)
			{
				double y = ((double)j-(sizeY-1)/2.0) / ((sizeY-1)/2.0);
				super.setAtIndex(i, super.getSize()[1]-j-1, 0.f);
				for (int k=0; k < Ellipses.getRows(); ++k)
				{
					// Extract the ellipse properties here
					double xc = x - Ellipses.getElement(k, 3);
					double yc = y - Ellipses.getElement(k, 4);
					double phi = Ellipses.getElement(k, 5)*Math.PI/180.0;
					double cos = Math.cos(phi);
					double sin = Math.sin(phi);
					double asq = Ellipses.getElement(k, 1)*Ellipses.getElement(k, 1);
					double bsq = Ellipses.getElement(k, 2)*Ellipses.getElement(k, 2);
					double Val = Ellipses.getElement(k, 0);
					
					// Check if this pixel is part of the ellipse, if yes, add the given intensity value to it
					double help = Math.pow((xc*cos + yc*sin),2.0);
					double help2 = Math.pow((yc*cos - xc*sin),2.0);
					if ( help/asq + help2/bsq <= 1.0 )
						super.setAtIndex(i, super.getSize()[1]-j-1, super.getAtIndex(i, super.getSize()[1]-j-1) + (float)Val);
				}
				
			}
		}
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		myPhantom test = new myPhantom(256, 256);
		test.show("My Phantom");
		
		

	}

}



