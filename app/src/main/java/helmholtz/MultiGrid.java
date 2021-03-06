package helmholtz;

public interface MultiGrid{

	
	public DCVector restrict(DCVector longvec);

	public DCVector interpolate(DCVector shortvec);

	public DCVector solve(DCVector vec);

	public DCVector smooth(DCVector roughvec, DCVector rhs, int ntimes);

	public int numpoints();


}