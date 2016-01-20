package helmholtz;

public class MGPreconditioner implements Preconditioner{
	private LexGrid finegrid;
	private FCycle fc;

	public MGPreconditioner(FDOperator op){
		finegrid = op.getGrid();

		// how many levels
		int shortdim = Math.min(finegrid.rows(), finegrid.cols());
		int dim = shortdim;
		int nlevels = 0;
		while (dim >3 && dim%2 != 0){
			nlevels++;
			dim = (dim+1)/2;
		}

		fc = new FCycle(op, nlevels, 1, 1);
	}
	
	public DCVector multiplyLeft(DCVector vec){
		return multiply(vec);
	}

	public DCVector multiplyRight(DCVector vec){
		return vec.copy();
	}

	public DCVector multiply(DCVector vec){
		return vec.copy();
	}

	public DCVector solveLeft(DCVector vec){
		return solve(vec);
	}

	public DCVector solveRight(DCVector vec){
		return vec.copy();
	}

	public DCVector solve(DCVector vec){
		return fc.solve(vec);
	}
}