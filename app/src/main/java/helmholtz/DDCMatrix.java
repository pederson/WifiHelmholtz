package helmholtz;

import java.util.Vector;

/*
 * Dense Double Complex Matrix implementation
*/
public class DDCMatrix extends Object{
	private Vector<Complex> values;
	private int mrows, ncols;	// matrix dimensions

	public DDCMatrix(int m, int n){
		mrows = m;
		ncols = n;

		values = new Vector<Complex>(mrows*ncols);
	}

	public DDCMatrix(int m, int n, Complex z){
		mrows = m;
		ncols = n;

		values = new Vector<Complex>(mrows*ncols);

		for (int i=0; i<mrows*ncols; i++){
			values.set(i, z.copy());
		}
	}

	public DDCMatrix(DDCMatrix mat){
		mrows = mat.rows();
		ncols = mat.cols();

		values = new Vector<Complex>(mrows*ncols);

		for (int i=0; i<mrows; i++){
			for (int j=0; j<ncols; j++){
				put(i, j, mat.at(i, j).copy());
			}
		}
	}

	public DDCMatrix copy(){
		return new DDCMatrix(this);
	}

	public int rows(){return mrows;};
	public int cols(){return ncols;};

	public void put(int i, int j, Complex z){
		values.set(j*ncols +i, z);
	}

	public Complex at(int i, int j){
		return values.get(j*ncols+i);
	}

}

