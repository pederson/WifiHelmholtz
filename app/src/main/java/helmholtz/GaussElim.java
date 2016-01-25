package helmholtz;


public class GaussElim{
	
	public GaussElim(){

	}

	public DCVector solve(Operator op, DCVector b){
		return solve(op.getDenseMatrix(), b);
	}

	public DCVector solve(DDCMatrix A, DCVector b){
		int m = A.rows();
		int n = A.cols();

		DDCMatrix L = DDCMatrix.eye(m,n);
		DDCMatrix U = A.copy();

		// construct L and U
		for (int k=0; k<m-1; k++){
			for (int i=k+1; i<m; i++){

				L.put(i,k,U.at(i,k).div(U.at(k,k)));

				for (int j=k; j<m; j++){
					U.put(i,j, U.at(i,j).minus(L.at(i,k).times(U.at(k,j))));
				}
			}
		}

		// forward substitution
		DCVector y = L.forward_solve(b);

		// back substitution
		DCVector out = U.back_solve(y);


		return out;
	}
}