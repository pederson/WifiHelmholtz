package helmholtz;


public interface FDOperator extends Operator{

	public DCVector multiply(DCVector v);

	public DDCMatrix getDenseMatrix();

	public SDCMatrix getSparseMatrix();

	public FDOperator coarsen();

	public LexGrid getGrid();
}