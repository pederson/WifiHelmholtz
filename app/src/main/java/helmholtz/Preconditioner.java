package helmholtz;

public interface Preconditioner{

	public DCVector multiplyLeft(DCVector vec);

	public DCVector multiplyRight(DCVector vec);

	public DCVector multiply(DCVector vec);

	public DCVector solveLeft(DCVector vec);

	public DCVector solveRight(DCVector vec);

	public DCVector solve(DCVector vec);
}