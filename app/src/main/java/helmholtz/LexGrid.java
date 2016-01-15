package helmholtz;

public class LexGrid implements MultiGrid{
	protected LexGrid supergrid;	// finer grid
	protected LexGrid subgrid;		// coarser grid
	protected int level;			// the grid coarseness depth level
	protected int mrows, ncols;

	public LexGrid(int rows, int cols){
		mrows = rows;
		ncols = cols;

	}

	public LexGrid coarsen(){
		int newrows, newcols;

		newrows = (mrows-1)/2;
		newcols = (ncols-1)/2;

		subgrid = new LexGrid(newrows, newcols);
		return subgrid;
	}

	public int numpoints(){
		return mrows*ncols;
	}

	public DCVector restrict(DCVector longvec){

	}

	public DCVector interpolate(DCVector shortvec){

	}

	public DCVector solve(DCVector vec){

	}

	public DCVector smooth(DCVector coarsevec){
		
	}

}