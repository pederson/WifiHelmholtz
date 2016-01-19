package helmholtz;


public class GaussElim{
	
	public GaussElim(){

	}

	public DCVector solve(Operator op, DCVector b){
		return solve(op.getDenseMatrix(), b);
	}

	public DCVector solve(DDCMatrix A, DCVector b){
		int r = A.rows();
		int c = A.cols();

		DCVector bp = b.copy();
		DDCMatrix U = A.copy();

		// for (int i=0; i<r; i++){
		// 	System.out.print("mat: ");
		// 	for (int j=0; j<c; j++){
		// 		System.out.print(""+U.at(i,j)+",");
		// 	}
		// 	System.out.println(" ");

		// }

		// for (int i=0; i<b.size(); i++){
		// 	System.out.println("rhs: "+b.at(i));
		// }

		// forward elimination
		for (int k=0; k<r-1; k++){
			for (int i=k+1; i<r; i++){
				U.put(i, k, U.at(i,k).div(U.at(k,k)));
				for (int j=k+1; j<c; j++){
					U.put(i,j, U.at(i,j).minus(U.at(i,k).times(U.at(k,j))));
				}
				// update rhs
				bp.put(i, bp.at(i).minus(U.at(i,k).times(bp.at(k))));
			}
		}


		// back substitution
		for (int i=r-1; i>=0; i--){
			for (int j=i+1; j<r; j++){
				bp.put(i, bp.at(i).minus(U.at(i,j).times(bp.at(j))));
			}
			bp.put(i, bp.at(i).div(U.at(i,i)));
		}

		return bp;
	}
}