package helmholtz;

public class JacobiPreconditioner implements Preconditioner{
	SDCMatrix K;

	public JacobiPreconditioner(SDCMatrix A){
		K = A;
	}


	/*
	 * Multiply left PC by a vector
	 *
	*/
	public DCVector multiplyLeft(DCVector vec){
		return multiply(vec);
	}

	/*
	 * Multiply right PC by a vector
	 *
	*/
	public DCVector multiplyRight(DCVector vec){
		return vec;
	}


	/*
	 * Multiply total PC by a vector
	 *
	*/
	public DCVector multiply(DCVector vec){
		return K.MatVec(vec);
	}

	
	/*
	 * Diagonal solve using K
	 *
	*/
	public DCVector solveLeft(DCVector vec){
		return solve(vec);
	}

	/*
	 * does nothing
	 *
	*/
	public DCVector solveRight(DCVector vec){
		return vec.copy();
	}

	/*
	 * Solve a general problem
	 *
	*/
	public DCVector solve(DCVector vec){
		return K.solveDiagonal(vec);
	}
}