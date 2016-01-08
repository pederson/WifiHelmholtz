package helmholtz;


public class BiCGSTAB extends Object{
	private int max_iters;
	private double tol;

	public BiCGSTAB(){
		max_iters = 300;
		tol = 1e-4;
	}

	public void set_max_iters(int n){max_iters = n;};

	public void set_tolerance(double t){tol = t;};

	public DCVector solve(SDCMatrix A, DCVector b, DCVector x){

		x = new DCVector(b.size());
		x.assign(new Complex(0.0, 0.0));
		DCVector r = b.minus(A.MatVec(x));
		DCVector rhat = r.copy();
		DCVector v = x.copy();
		DCVector p = x.copy();
		DCVector s, t;

		Complex rho = new Complex(1.0, 0.0);
		Complex alpha = rho.copy();
		Complex omega = rho.copy();
		Complex beta;
		Complex rholast = rho.copy();

		// calculate initial residual magnitude
		Complex residm = r.dot(r);
		double resid0 = residm.mod();
		double resid = resid0;
		
		int it=0;
		while (it < max_iters && resid/resid0 > tol){

			// update some constants
			rholast = rho;
			rho = rhat.dot(r);
			beta = (rho.div(rholast)).times((alpha.div(omega)));

			// calculate new search direction
			p = r.plus((p.minus(v.times(omega))).times(beta));

			// matvec
			v = A.MatVec(p);

			// new step size
			alpha = rho.div(rhat.dot(v));

			// some other vars
			s = r.minus(v.times(alpha));
			t = A.MatVec(s);
			omega = (t.dot(s)).div(t.dot(t));
			x = x.plus((p.times(alpha)).plus(s.times(omega)));

			// new residual
			r = s.minus(t.times(omega));
			residm = r.dot(r);
			resid = residm.mod();

			System.out.println("iter: "+it+" resid/r0: "+(resid/resid0));

			it++;
		}

		return x;
	}


	public DCVector solve(Preconditioner pc, SDCMatrix A, DCVector b, DCVector x){

		x = new DCVector(b.size());
		x.assign(new Complex(0.0, 0.0));
		DCVector r = b.minus(A.MatVec(x));
		DCVector rhat = r.copy();
		DCVector v = x.copy();
		DCVector y;
		DCVector p = x.copy();
		DCVector s, t, z;
		DCVector k1inv_s, k1inv_t;

		Complex rho = new Complex(1.0, 0.0);
		Complex alpha = rho.copy();
		Complex omega = rho.copy();
		Complex beta;
		Complex rholast = rho.copy();

		// calculate initial residual magnitude
		Complex residm = r.dot(r);
		double resid0 = residm.mod();
		double resid = resid0;
		
		int it=0;
		while (it < max_iters && resid/resid0 > tol){

			// update some constants
			rholast = rho;
			rho = rhat.dot(r);
			beta = (rho.div(rholast)).times((alpha.div(omega)));

			// calculate new search direction
			p = r.plus((p.minus(v.times(omega))).times(beta));

			// use pc
			y = pc.solve(p);

			// matvec
			v = A.MatVec(y);

			// new step size
			alpha = rho.div(rhat.dot(v));

			// some other vars
			s = r.minus(v.times(alpha));
			k1inv_s = pc.solveLeft(s);
			z = pc.solveRight(k1inv_s);
			t = A.MatVec(z);
			k1inv_t = pc.solveLeft(t);
			omega = (k1inv_t.dot(k1inv_s)).div(k1inv_t.dot(k1inv_t));
			x = x.plus((y.times(alpha)).plus(z.times(omega)));

			// new residual
			r = s.minus(t.times(omega));
			residm = r.dot(r);
			resid = residm.mod();

			System.out.println("iter: "+it+" resid/r0: "+(resid/resid0));

			it++;
		}

		return x;
	}
}