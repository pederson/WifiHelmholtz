package helmholtz;

import cern.colt.matrix.tdcomplex.impl.SparseDComplexMatrix2D;
import java.util.concurrent.*;

public class SDCMatrix2D extends SparseDComplexMatrix2D{

	public SDCMatrix2D(int m, int n){
		super(m,n);
	}

	public SDCMatrix2D(int m, int n, ConcurrentHashMap<Long, double[]> map, int i, int d, int not, int know){
		super(m, n, map, i, d, not, know);
	}
}