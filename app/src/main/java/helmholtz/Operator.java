package helmholtz;


public interface Operator{

	public DCVector multiply(DCVector v);

	public DDCMatrix getDenseMatrix();

	public SDCMatrix getSparseMatrix();
}