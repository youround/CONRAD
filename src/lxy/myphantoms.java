package lxy;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.NumericGrid;
import edu.stanford.rsl.conrad.data.numeric.NumericPointwiseOperators;
import edu.stanford.rsl.tutorial.phantoms.MickeyMouseGrid2D;
import ij.ImageJ;

public class myphantoms extends Grid2D {
	private String title;

	public myphantoms(int x, int y) {
	this(x, y, 4);
}

	public myphantoms(int x, int y, int phase) {
	this(x,y,"myphantoms");
	
	int val1 = 1;
	int val2 = 3;
	int val3 = 5;

	double r1 = 0.4*x;
	double r2 = 0.1*x;
	double w1 = 0.2*x;
	double h1 = 0.06*y;
	
	int xCenter1 = x/2;
	int yCenter1 = 6*y/10;
	
	int xCenter2 = (int) (0.25*x);
	int yCenter2 = (int) (0.55*y);
	int xCenter3 = (int) (0.75*x);
	int yCenter3 = (int) (0.55*y);
	
	int xCenter4 = (int) (0.5*x);
	int yCenter4 = (int) (0.8*y);
	
	for(int i = 0; i < x; i++){
		for(int j = 0; j < y; j++) {
		
			if (phase > 3) if( Math.pow(i - xCenter1, 2)  + Math.pow(j - yCenter1, 2) <= (r1*r1) ) {
				super.setAtIndex(i, j, val1);
			}
			if (phase > 2) if( Math.pow(i - xCenter2, 2)  + Math.pow(j - yCenter2, 2) <= (r2*r2) ) {
				super.setAtIndex(i, j, val2);
			}
			if (phase > 1) if( Math.pow(i - xCenter3, 2)  + Math.pow(j - yCenter3, 2) <= (r2*r2) ) {
				super.setAtIndex(i, j, val2);
			}	
			if (phase > 0) if( Math.abs(i - xCenter4) <= w1  && Math.abs(j - yCenter4) <= h1 ) {
				super.setAtIndex(i, j, val3);			
			}
		}
	}
}
	public myphantoms(int x, int y, String title) {
		// TODO Auto-generated constructor stub
		super(new float[x*y], x, y);
		this.setSpacing(1,1);
		this.setOrigin(x/2.d,y/2.d);
		this.title = title;
	}

	public void setTitle(String t) {
		this.title = t;
	}
	
	public String getTitle() {
		return this.title;
	}
	
	public static void main (String [] args){
		new ImageJ();
		myphantoms test = new myphantoms(256,256);
		MickeyMouseGrid2D test2 = new MickeyMouseGrid2D(256, 256, 3);

		test.show("myphantoms");
		test2.show("micky");
		NumericGrid res = NumericPointwiseOperators.addedBy(test, test2);
		float m = NumericPointwiseOperators.mean(test);
		res.show("result");
		System.out.println("mean is " + m);
		
		
		}

}
