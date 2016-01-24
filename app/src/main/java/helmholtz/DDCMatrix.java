package helmholtz;

import java.util.Vector;
import java.util.Random;

/*
 * Dense Double Complex Matrix implementation
*/
public class DDCMatrix extends Object{
	private Complex[] values;
	private int mrows, ncols;	// matrix dimensions

	public DDCMatrix(int m, int n){
		mrows = m;
		ncols = n;

		values = new Complex[mrows*ncols];
	}

	public DDCMatrix(int m, int n, Complex z){
		mrows = m;
		ncols = n;

		values = new Complex[mrows*ncols];

		for (int i=0; i<mrows*ncols; i++){
			values[i] = z.copy();
		}
	}

	public DDCMatrix(DDCMatrix mat){
		mrows = mat.rows();
		ncols = mat.cols();

		values = new Complex[mrows*ncols];

		for (int i=0; i<mrows; i++){
			for (int j=0; j<ncols; j++){
				values[j*ncols+i] = mat.at(i, j).copy();
			}
		}
	}

	public DDCMatrix copy(){
		return new DDCMatrix(this);
	}

	public static DDCMatrix eye(int m, int n){
		DDCMatrix out = new DDCMatrix(m,n,new Complex(0,0));
		int mindim = Math.min(m,n);
		for (int i=0; i<mindim; i++){
			out.put(i,i,new Complex(1,0));
		}
		return out;
	}

	public static DDCMatrix rand(int m, int n){
		Random rgen = new Random();
		DDCMatrix out = new DDCMatrix(m,n);

		for (int i=0; i<m; i++){
			for (int j=0; j<n; j++){
				out.put(i,j, new Complex(rgen.nextDouble(), rgen.nextDouble()));
			}
		}
		return out;
	}

	public DCVector MatVec(DCVector v){
		DCVector out = new DCVector(mrows);
		Complex s;

		// add up each row
		for (int i=0; i<mrows; i++){
			s = new Complex(0,0);
			for (int j=0; j<ncols; j++){
				s = s.plus(at(i,j).times(v.at(j)));
			}
			out.put(i,s);
		}
		return out;
	}

	public DCVector forward_solve(DCVector v){
		System.out.println("FORWARD!!!!!!!!!!!!!!");
		DCVector out = new DCVector(mrows);
		Complex s;

		for (int j=0; j<mrows; j++){
			s = new Complex(0,0);
			for (int k=0; k<j; k++){
				s = s.plus(out.at(k).times(at(j,k)));
			}
			out.put(j, (v.at(j).minus(s)).div(at(j,j)));
		}
		return out;
	}

	public DCVector back_solve(DCVector v){
		System.out.println("BACKWARD!!!!!!!!!!!!!!");
		DCVector out = new DCVector(mrows);
		Complex s;

		for (int j=mrows-1; j>=0; j--){
			s = new Complex(0,0);
			for (int k=j+1; k<mrows; k++){
				s = s.plus(out.at(k).times(at(j,k)));
			}
			out.put(j, (v.at(j).minus(s)).div(at(j,j)));
		}
		return out;
	}

	public int rows(){return mrows;};
	public int cols(){return ncols;};

	public void put(int i, int j, Complex z){
		values[j*ncols +i] = z;
	}

	public void put(IndexPair idx, Complex z){
		values[idx.j*ncols+idx.i] = z;
	}

	public Complex at(int i, int j){
		return values[j*ncols+i];
	}

}

