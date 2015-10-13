package lxy.FPR;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.numerics.SimpleOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.conrad.utils.CONRAD;

public class ContourTracer {
	
	protected int z;
	protected Grid2D img;
	protected Grid2D contour;
	private ArrayList<Point> contourPoints;
	private static final Map<Point, Point> clockwiseOffset;
	static
	{
		Map<Point, Point> aMap = new HashMap<Point, Point>();
		aMap.put(new Point(1,0), new Point(1,-1));     // right        => down-right
		aMap.put(new Point(1,-1), new Point(0,-1));    // down-right   => down
		aMap.put(new Point(0,-1), new Point(-1,-1));   // down         => down-left
		aMap.put(new Point(-1,-1), new Point(-1,0));   // down-left    => left
		aMap.put(new Point(-1,0), new Point(-1,1));    // left         => top-left
		aMap.put(new Point(-1,1), new Point(0,1));     // top-left     => top
		aMap.put(new Point(0,1), new Point(1,1));      // top          => top-right
		aMap.put(new Point(1,1), new Point(1,0));      // top-right    => right

		clockwiseOffset = Collections.unmodifiableMap(aMap);	
	}

	
	public ContourTracer(Grid2D img, int z)
	{
		this.img = img;
		contour = new Grid2D(img.getHeight(), img.getWidth());
		contourPoints = new ArrayList<Point>();
		this.z = z;
	
	}
	
	public ArrayList<PointND> getCountourPoints()
	{
		ArrayList<PointND> pND = new ArrayList<PointND>();;
		for (Point p : contourPoints) 
		{
			pND.add(new PointND(p.getX(),p.getY(), z));
		}
   
		    
		return pND;
	}
	
	public Grid2D mooreNeighborhoodTracer(int direc)
	{
		
		Point prev = new Point(0, 0);
		Point curr = new Point(0, 0);
		Point start =  new Point(0, 0);
		Point boundary =  new Point(0, 0);
					
		//fi the first black point
		boolean allwhite = true;
		if (direc == 0)
		{
			findStart:
			for (int i = 0; i < img.getWidth(); i ++)
			{
				for(int j = 0; j < img.getHeight(); j ++)
				{
					if(img.getAtIndex(i, j) == 0)
					{
						start = new Point(i, j);
						prev = new Point(i, j - 1);
					    contourPoints.add(start);
					    contour.setAtIndex(i, j, 1000);
					    allwhite = false;
					    break findStart;				    
					}			
				}
			}
		}
		else if (direc == 1)
		{
		    findStart:
			for (int i = img.getWidth() - 1; i >= 0; i --)
			{
				for(int j = 0; j < img.getHeight(); j ++)
				{
					if(img.getAtIndex(i, j) == 0)
					{
						start = new Point(i, j);
						prev = new Point(i, j - 1);
					    contourPoints.add(start);
					    contour.setAtIndex(i, j, 1000);
					    allwhite = false;
					    break findStart;				    
					}			
				}
			}
	}
		
		if(!allwhite)
		{
			boundary = start;
			curr = Clockwise(boundary, prev);
			while(!curr.equals(start))
			{
				int x = curr.getX();
				int y = curr.getY();
				if(x < img.getWidth() && x >=0 && y < img.getHeight() && y >= 0 && img.getAtIndex(x, y) == 0)
				{
					contourPoints.add(curr);
					contour.setAtIndex(x, y, 1000);
					prev = boundary;
					boundary = curr;
		            curr = Clockwise(boundary, prev);
				}
				else
				{
					prev = curr;
	                curr = Clockwise(boundary, prev);
				}
			}
		}
		return contour;
     
	}
	
	public Grid2D resample()
	{
		Grid2D pointsMap = new Grid2D(img.getHeight(), img.getWidth());
		ArrayList<Point> output = new ArrayList<Point>();
		int pointNum = 20;
	    int slicePointnum = contourPoints.size();
		if(slicePointnum > 300)
		{
			int step = 18;
			for(int j = 0; j < pointNum / 2 ; j ++)
			{
				Point cur = contourPoints.get(j * step);
				output.add(cur);
				pointsMap.setAtIndex(cur.getX(), cur.getY(), 1000);
			}
			for(int j = pointNum - pointNum / 2 ; j >= 0 ; j --)
			{
				Point cur = contourPoints.get(slicePointnum - j * step - 1);
				output.add(cur);
				pointsMap.setAtIndex(cur.getX(), cur.getY(), 1000);
			}
				
			}
			else
			{
				int step = slicePointnum/pointNum;
				for(int j = 0; j < pointNum * step; j += step)
				{
					Point cur = contourPoints.get(j);
					output.add(cur);
					pointsMap.setAtIndex(cur.getX(), cur.getY(), 1000);
				}
			}
		   // pointsMap.show("Resample");
		    return pointsMap;
	}


	private static Point Clockwise(Point target, Point prev)
	{
		Point next = clockwiseOffset.get(prev.subtract(target)).add(target);
	    return next;
	}
	
	private static Point key(Point x)
	{
		Point next = clockwiseOffset.get(x);
		System.out.println(next);
	    return x;
	}

	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Map<Point, Point> aMap = new HashMap<Point, Point>();
		aMap.put(new Point(1,0), new Point(1,-1)); 
		System.out.println(aMap.size());
		System.out.println(aMap.get(new Point(1,0)));

	}
	

}


class Point
{
	private int x, y;
	public Point(int x, int y)
	{
		this.x = x;
		this.y = y;
	}
	public void setCoordinates(int x, int y)
	{
		this.x = x;
		this.y = y;
	}
	public int getX()
	{
		return x;
	}
	public int getY()
	{
		return y;
	}
	public Point add(Point p)
	{
		return new Point(x + p.x, y + p.y);
	}
	public Point subtract(Point p)
	{
		return new Point(x - p.x, y - p.y);
	}
	@Override
	public boolean equals(Object o) 
	{
		if (o instanceof Point) 
		{
			Point p = (Point) o;
			return (x == p.x && y == p.y);
		} 
		else 
		{
			return false;
		}
	}
	public int hashCode() 
	{
	    return (x * 31) ^ y;
	}
	public String toString() { 
	    return "(" + x + ", " + y + ")";
	} 
	
	
	
}


