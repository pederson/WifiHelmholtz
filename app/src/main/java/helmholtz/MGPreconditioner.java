package helmholtz;

public class MGPreconditioner implements Preconditioner{
	private LexGrid finegrid;
	private FCycle fc;

	public MGPreconditioner(int rows, int cols, double dx){
		finegrid = new LexGrid(rows, cols, dx);

		// how many levels
		int shortdim = Math.min(rows, cols);
		int dim = shortdim;
		int nlevels = 0;
		while (dim >1){
			nlevels++;
			dim = (dim-1)/2+1;
		}

		fc = new FCycle(finegrid, nlevels, 1, 1);
	}
	
	public DCVector multiplyLeft(DCVector vec){
		return multiply(vec);
	}

	public DCVector multiplyRight(DCVector vec){
		return multiply(vec);
	}

	public DCVector multiply(DCVector vec){
		return vec;
	}

	public DCVector solveLeft(DCVector vec){
		return solve(vec);
	}

	public DCVector solveRight(DCVector vec){
		return solve(vec);
	}

	public DCVector solve(DCVector vec){
		return fc.solve(vec);
	}
}