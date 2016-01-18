package helmholtz;

import java.util.Vector;

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

	public int rows(){return mrows;};
	public int cols(){return ncols;};

	public void put(int i, int j, Complex z){
		values[j*ncols +i] = z;
	}

	public Complex at(int i, int j){
		return values[j*ncols+i];
	}

}

