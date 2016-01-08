package helmholtz;

public class ILUPreconditioner implements Preconditioner{
	SDCMatrix L, U;

	public ILUPreconditioner(SDCMatrix A){

	}

	protected void incompleteLU(SDCMatrix A, SDCMatrix La, SDCMatrix Ua){

	}

	/*
	 * Multiply left PC by a vector
	 *
	*/
	public DCVector multiplyLeft(DCVector vec){
		return L.MatVec(vec);
	}

	/*
	 * Multiply right PC by a vector
	 *
	*/
	public DCVector multiplyRight(DCVector vec){
		return U.MatVec(vec);
	}


	/*
	 * Multiply total PC by a vector
	 *
	*/
	public DCVector multiply(DCVector vec){
		return L.MatVec(U.MatVec(vec));
	}

	/*
	 * Solve a general problem
	 *
	*/
	public DCVector solve(DCVector vec){
		DCVector result = new DCVector(vec.size());

		return result;
	}

	/*
	 * Performs forward substituion using L
	 *
	*/
	public DCVector solveLeft(DCVector vec){
		DCVector result = new DCVector(vec.size());

		return result;
	}

	/*
	 * Performs back substituion using U
	 *
	*/
	public DCVector solveRight(DCVector vec){
		DCVector result = new DCVector(vec.size());

		return result;
	}
}