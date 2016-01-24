package helmholtz;

import java.util.Vector;

public class FCycle{
	protected Vector<LexGrid> grids;
	protected int n1, n2;
	protected int nlevels;
	protected FDOperator op0, opl;

	FCycle(FDOperator op, int nlev, int nu1, int nu2){
		nlevels = nlev;
		opl = op;
		LexGrid gfine = op.getGrid();
		grids = new Vector<LexGrid>(nlevels+1);
		for (int i=0; i<=nlevels; i++){
			grids.addElement(gfine);
		}

		// create multigrids and get coarse operator
		FDOperator ophold = opl;
		grids.set(nlevels, gfine);
		LexGrid g;
		for (int i=nlevels; i>0; i--){
			op0 = ophold.coarsen();
			g = grids.get(i).coarsen();
			grids.set(i-1, g);
			ophold = op0;
		}

		// number of smoothing steps parameters
		n1 = nu1;
		n2 = nu2;
	}

	public DCVector solve(DCVector vec){
		DCVector err, sln, resid;
		//DCVector v = vec.copy();	// rhs vector
		DCVector bl, rhs, tmp;
		GaussElim gel = new GaussElim();

		
		//grids.get(nlevels).print_vector(vec);

		// initial solution assumed to be zero at all levels
		Vector<DCVector> solns = new Vector<DCVector>(nlevels+1);
		for (int l=0; l<=nlevels; l++){
			bl = new DCVector(grids.get(l).numpoints());
			bl.assign(new Complex(0, 0));
			solns.addElement(bl);
		}
		// System.out.println("initial solutions: ");
		// for(int i=0; i<=nlevels; i++){
		// 	sln = solns.get(i);
		// 	System.out.print("level: "+i);
		// 	System.out.println(" sln: "+sln.dot(sln).toString());
		// }

		// rhs
		Vector<DCVector> rhss = new Vector<DCVector>(nlevels+1);
		for (int i=0; i<=nlevels; i++) rhss.addElement(vec);
		rhss.set(nlevels, vec.copy());

		// initial descent to the bottom level
		LexGrid g;
		for (int l=nlevels; l>0; l--){
			System.out.println(" ");
			System.out.println("level: "+l);
			g = grids.get(l);
			sln = solns.get(l);
			rhs = rhss.get(l);

			System.out.println("rhs: "+rhs.dot(rhs).toString());

			// pre-smooth the solution
			System.out.println("rough sln: "+sln.dot(sln).toString());
			sln = g.smooth(sln, rhs, n1);
			System.out.println("smoothed sln: "+sln.dot(sln).toString());
			solns.set(l,sln);

			// compute residual
			resid = rhs.minus(apply_operator(l, sln));
			System.out.println("resid: "+resid.dot(resid).toString());

			// restrict the residual
			err = g.restrict(resid);

			// update the containers
			// the restricted residual becomes the new RHS
			rhss.set(l-1, err);

		}

		// solve at the bottom by gaussian elimination
		solns.set(0, gel.solve(op0, rhss.get(0)));
		System.out.println(" ");
		System.out.println("..... BOTTOM LEVEL .....");
		resid = rhss.get(0).minus(apply_operator(0,solns.get(0)));
		System.out.println("resid: "+resid.dot(resid).toString());

		// cycle upwards with increasing reach
		for (int l=1; l<nlevels; l++){

			// ascend
			System.out.println(" ");
			System.out.println(".....ASCENDING.....");
			for (int u=1; u<=l; u++){
				System.out.println(" ");
				System.out.println("level: "+u);
				g = grids.get(u);
				sln = solns.get(u);
				rhs = rhss.get(u);

				// interpolate the error
				err = grids.get(u-1).interpolate(solns.get(u-1));

				// update the solution with error
				System.out.println("old sln: "+sln.dot(sln).toString());
				sln = sln.plus(err);
				
				// post-smooth
				sln = g.smooth(sln, rhs, n2);
				System.out.println("new sln: "+sln.dot(sln).toString());


				// update the containers
				solns.set(u, sln);

			}

			// descend
			System.out.println(" ");
			System.out.println(".....DESCENDING.....");
			for (int d=l; d>0; d--){
				System.out.println(" ");
				System.out.println("level: "+d);
				g = grids.get(d);
				sln = solns.get(d);
				rhs = rhss.get(d);

				System.out.println("rhs: "+rhs.dot(rhs).toString());


				// pre-smooth the solution
				System.out.println("rough sln: "+sln.dot(sln).toString());
				sln = g.smooth(sln, rhs, n1);
				System.out.println("smoothed sln: "+sln.dot(sln).toString());

				// compute residual
				resid = rhs.minus(apply_operator(d, sln));
				System.out.println("resid: "+resid.dot(resid).toString());


				// restrict the residual
				err = g.restrict(resid);

				// update the containers
				// the restricted residual becomes the new RHS
				rhss.set(d-1, err);
				solns.set(d,sln);

			}

			// solve at the bottom by gaussian elimination
			solns.set(0, gel.solve(op0, rhss.get(0)));
		}

		// make final ascent to the top level
		for (int u=1; u<=nlevels; u++){
				g = grids.get(u);
				sln = solns.get(u);
				rhs = rhss.get(u);

				// interpolate the error
				err = grids.get(u-1).interpolate(solns.get(u-1));

				// update the solution with error
				sln = sln.plus(err);

				// post-smooth
				sln = g.smooth(sln, rhs, n2);

				// update the containers
				solns.set(u, sln);

		}


		return solns.get(nlevels);
	}


	protected DCVector apply_operator(int level, DCVector vec){

		LexGrid g;
		DCVector rvec, ivec, outr;
		//int nlevels = grids.size();

		ivec = vec;
		// interpolate up to top level
		for (int l=level; l<nlevels; l++){
			g = grids.get(l);
			ivec = g.interpolate(ivec);
		}

		// apply operator
		outr = opl.multiply(ivec);

		rvec = outr;
		// restrict back to original level
		for (int l=nlevels; l>level; l--){
			//System.out.println("restrict level:"+l);
			g = grids.get(l);
			rvec = g.restrict(rvec);
		}

		//System.out.println("rvec: "+rvec.size());
		return rvec;
	}


}