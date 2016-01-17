package helmholtz;

import java.util.Vector;

public class FCycle{
	protected Vector<LexGrid> grids;
	protected int n1, n2;
	protected int nlevels;

	FCycle(LexGrid gfine, int nlev, int nu1, int nu2){
		nlevels = nlev;
		grids = new Vector<LexGrid>(nlevels+1);
		for (int i=0; i<=nlevels; i++){
			grids.addElement(gfine);
		}

		// create multigrids
		grids.set(nlevels, gfine);
		LexGrid g;
		for (int i=nlevels; i>0; i--){
			g = grids.get(i).coarsen();
			grids.set(i-1, g);
		}

		// number of smoothing steps parameters
		n1 = nu1;
		n2 = nu2;
	}

	public DCVector solve(DCVector vec){
		DCVector err, sln, resid;
		//DCVector v = vec.copy();	// rhs vector
		DCVector bl, rhs;

		// initial solution assumed to be zero at all levels
		Vector<DCVector> solns = new Vector<DCVector>(nlevels+1);
		for (int l=nlevels; l>=0; l--){
			bl = new DCVector(grids.get(l).numpoints());
			bl.assign(new Complex(0, 0));
			solns.add(bl);
		}

		// rhs
		Vector<DCVector> rhss = new Vector<DCVector>(nlevels+1);
		rhss.set(nlevels, vec.copy());

		// initial descent to the bottom level
		LexGrid g;
		for (int l=nlevels; l>0; l++){
			g = grids.get(l);
			sln = solns.get(l);
			rhs = rhss.get(l);

			// pre-smooth the solution
			sln = g.smooth(sln, rhs, n1);

			// compute residual
			resid = rhs.minus(apply_operator(l, sln));

			// restrict the residual
			err = g.restrict(resid);

			// the restricted residual becomes the new RHS
			rhss.set(l-1, err);
		}

		// solve at the bottom
		g = grids.get(0);
		solns.set(0, g.solve(rhss.get(0)));

		// cycle upwards with increasing reach
		for (int l=1; l<nlevels; l++){

			// ascend
			for (int u=1; u<=l; u++){
				g = grids.get(u);
				sln = solns.get(u);
				rhs = rhss.get(u);

				// interpolate the error
				err = g.interpolate(solns.get(u-1));

				// update the solution with error
				sln = sln.plus(err);

				// post-smooth
				sln = g.smooth(sln, rhs, n2);

			}

			// descend
			for (int d=l; d>0; d--){
				g = grids.get(d);
				sln = solns.get(d);
				rhs = rhss.get(d);

				// pre-smooth the solution
				sln = g.smooth(sln, rhs, n1);

				// compute residual
				resid = rhs.minus(apply_operator(l, sln));

				// restrict the residual
				err = g.restrict(resid);

				// the restricted residual becomes the new RHS
				rhss.set(d-1, err);

			}

			// solve at the bottom
			g = grids.get(0);
			solns.set(0, g.solve(rhss.get(0)));
		}

		// make final ascent to the top level
		for (int u=1; u<=nlevels; u++){
				g = grids.get(u);
				sln = solns.get(u);
				rhs = rhss.get(u);

				// interpolate the error
				err = g.interpolate(solns.get(u-1));

				// update the solution with error
				sln = sln.plus(err);

				// post-smooth
				sln = g.smooth(sln, rhs, n2);

		}


		return solns.get(nlevels);
	}


	protected DCVector apply_operator(int level, DCVector vec){

		LexGrid g;
		DCVector rvec, ivec, outr;
		int nlevels = grids.size()-1;

		ivec = vec;
		// interpolate up to top level
		for (int l=level; l<nlevels; l++){
			g = grids.get(l);
			ivec = g.interpolate(ivec);
		}

		// apply operator
		// THIS IS WRONG!!!
		outr = ivec;

		rvec = outr;
		// restrict back to original level
		for (int l=nlevels; l>level; l--){
			g = grids.get(l);
			rvec = g.restrict(rvec);
		}

		return rvec;
	}

}