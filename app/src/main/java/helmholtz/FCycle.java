package helmholtz;

public class FCycle{
	protected Vector<LexGrid> grids;
	protected int n1, n2;

	FCycle(LexGrid gfine, int nlevels, int nu1, int nu2){
		grids = new Vector<LexGrid>(nlevels+1);

		// create multigrids
		grids.add(nlevels, gfine);
		LexGrid g;
		for (int i=nlevels; i>0; i--){
			g = grids.get(i).coarsen();
			grids.add(i-1, g);
		}

		// number of smoothing steps parameters
		n1 = nu1;
		n2 = nu2;
	}

	public DCVector solve(DCVector vec){
		DCVector err, sln, resid;
		DCVector v = vec.copy();	// rhs vector
		DCVector bl;

		// initial solution assumed to be zero at all levels
		Vector<DCVector> solns = new Vector<DCVector>(nlevels+1);
		for (int l=nlevels; l>=0; l--){
			bl = new DCVector(grids.get(l).numpoints());
			bl.assign(new Complex(0, 0));
			solns.add(bl);
		}

		// initial descent to the bottom level
		LexGrid g;
		for (int l=nlevels; l>0; l++){
			g = grids.get(l);
			sln = solns.get(l);

			// pre-smooth the solution
			sln = g.smooth(sln);

			// compute residual
			resid = v.minus(apply_operator(l, sln));

			// restrict the residual
			err = g.restrict(resid);

			// the restricted residual becomes the new RHS
			v = resid;
		}

		// solve at the bottom
		g = grids.get(0);
		solns.get(0) = g.solve(v);

		// cycle upwards with increasing reach
		for (int l=1; l<nlevels; l++){

			// ascend
			for (int u=1; u<=l; u++){
				g = grids.get(u);
				sln = solns.get(u);

				// interpolate the error
				err = g.interpolate(solns.get(u-1));

				// update the solution with error
				sln = sln.plus(err);

				// post-smooth
				sln = g.smooth(sln);

			}

			// descend
			for (int d=l; d>0; d--){
				g = grids.get(d);
				sln = solns.get(d);

				// pre-smooth the solution
				sln = g.smooth(sln);

				// compute residual
				resid = v.minus(apply_operator(l, sln));

				// restrict the residual
				err = g.restrict(resid);

				// the restricted residual becomes the new RHS
				v = resid;

			}

			// solve at the bottom
			g = grids.get(0);
			solns.get(0) = g.solve(v);
		}

		// make final ascent to the top level
		for (int u=1; u<=nlevels; u++){
				g = grids.get(u);
				sln = solns.get(u);

				// interpolate the error
				err = g.interpolate(solns.get(u-1));

				// update the solution with error
				sln = sln.plus(err);

				// post-smooth
				sln = g.smooth(sln);

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
		//outr = 

		rvec = outr;
		// restrict back to original level
		for (int l=nlevels; l>level; l--){
			g = grids.get(l);
			rvec = g.restrict(rvec);
		}

	}

}