package helmholtz;

public class FCycle{
	protected Vector<LexGrid> grids;
	protected int n1, n2;

	FCycle(LexGrid gfine, int nlevels, int nu1, int nu2){
		grids = new Vector<LexGrid>(nlevels+1);

		// create multigrids
		grids.add(0, gfine);
		LexGrid g;
		for (int i=0; i<nlevels; i++){
			g = grids.get(i).coarsen();
			grids.add(i=1, g);
		}

		// number of smoothing steps parameters
		n1 = nu1;
		n2 = nu2;
	}

	public DCVector solve(DCVector vec){
		DCVector err;
		DCVector v = vec.copy();

		// initial solution assumed to be zero

		// calculate the initial error

		// initial descent to the bottom level
		LexGrid g;
		for (int l=0; l<nlevels; l++){
			g = grids.get(l);

			// smooth
			err = g.smooth(err);

			// then restrict
			err = g.restrict(err);
		}

		// solve at the bottom
		g = grids.get(nlevels);
		g.solve(err);

		// cycle upwards with increasing reach
		for (int l=nlevels-1; l>0; l--){

			// ascend
			for (int u=nlevels-1; u>=l; u--){

			}

			// descend
			for (int d=l; d<nlevels; d++){

			}

			// solve at the bottom
		}

		// make final ascent to the top level
	}

}