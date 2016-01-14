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

	}

}