package helmholtz;

import java.util.HashMap;

public class HelmholtzOperator2D implements FDOperator{
    protected double c0 = 2.99e+8;
	protected DCVector kv;
	protected int wid, hei;
	protected double dx;
	protected double omega;
	protected SDCMatrix smat;
	protected DDCMatrix dmat;
	protected LexGrid grid;

	public HelmholtzOperator2D(LexGrid g, double freq, DCVector kvec){
		kv = kvec;
		wid = g.rows();
		hei = g.cols();
		dx = g.res();
		grid = g;
		omega = freq*2*Math.PI;
	}

	public DCVector getKsq(){
		return kv;
	}

	public LexGrid getGrid(){return grid;};

	public HelmholtzOperator2D coarsen(){
		// if (kv == null) System.out.println("NULL KVECTOR");
		// System.out.println("grid size: "+(grid.rows()*grid.cols())+" kvector size: "+kv.size());
		DCVector knew = grid.restrict(kv);
		return new HelmholtzOperator2D(grid.coarsen(), omega/(2.0*Math.PI), knew);
	}

	public DCVector multiply(DCVector v){
		SDCMatrix s;
		if (smat == null){
			s = getSparseMatrix();
		}
		else{
			s = smat;
		}

		return s.MatVec(v);
	}

	public DDCMatrix getDenseMatrix(){
		if (dmat != null){
			return dmat;
		}

		int cind, lind, rind, uind, dind;

		int m = wid*hei;
		IndexPair idx = new IndexPair(0,0);

		DDCMatrix mat = new DDCMatrix(m, m, new Complex(0,0));
        Complex val = new Complex(-1.0 / Math.pow(dx,2), 0.0);

		for (int i=1; i<wid-1; i++){
			for (int j=1; j<hei-1; j++){

                cind = inds_to_global(i,j);
                lind = inds_to_global(i-1,j);
                rind = inds_to_global(i+1,j);
                uind = inds_to_global(i,j+1);
                dind = inds_to_global(i,j-1);

                // set the HashMap values
                // center
                idx = new IndexPair(cind, cind);
                mat.put(idx, (kv.at(cind).times(-1.0)).plus(4.0/(dx*dx)));

                // left
                idx = new IndexPair(cind, lind);
                mat.put(idx, val.copy());

                // right
                idx = new IndexPair(cind, rind);
                mat.put(idx, val.copy());

                // up
                idx = new IndexPair(cind, uind);
                mat.put(idx, val.copy());

                // down
                idx = new IndexPair(cind, dind);
                mat.put(idx, val.copy());

    //             //
				// mat.put(cind, cind, kv.at(cind));
				// mat.put(cind, lind, new Complex(1.0/(dx*dx), 0.0));
				// mat.put(cind, rind, new Complex(1.0/(dx*dx), 0.0));
				// mat.put(cind, uind, new Complex(1.0/(dx*dx), 0.0));
				// mat.put(cind, dind, new Complex(1.0/(dx*dx), 0.0));
			}
		}

		// FIRST ORDER absorbing boundaries
        // top boundary
        for (int i=0; i<wid-1; i++){

            cind = inds_to_global(i,hei-1);
            dind = inds_to_global(i,hei-2);
            
            idx = new IndexPair(cind, cind);
            mat.put(idx, new Complex(1.0/dx, -omega/c0));
        
            idx = new IndexPair(cind, dind);
            mat.put(idx, new Complex(-1.0/dx, 0.0));

        }

        // bottom boundary
        for (int i=1; i<wid; i++){
            cind = inds_to_global(i,0);
            uind = inds_to_global(i,1);
            
            idx = new IndexPair(cind, cind);
            mat.put(idx, new Complex(-1.0/dx, -omega/c0));
        
            idx = new IndexPair(cind, uind);
            mat.put(idx, new Complex(1.0/dx, 0.0));

        }

        // left boundary
        for (int j=0;j<hei-1; j++){
            cind = inds_to_global(0,j);
            rind = inds_to_global(1,j);
            
            idx = new IndexPair(cind, cind);
            mat.put(idx, new Complex(-1.0/dx, -omega/c0));
        
            idx = new IndexPair(cind, rind);
            mat.put(idx, new Complex(1.0/dx, 0.0));

        }

        // right boundary
        for (int j=1; j<hei; j++){
            cind = inds_to_global(wid-1,j);
            lind = inds_to_global(wid-2,j);
            
            idx = new IndexPair(cind, cind);
            mat.put(idx, new Complex(1.0/dx, -omega/c0));
        
            idx = new IndexPair(cind, lind);
            mat.put(idx, new Complex(-1.0/dx, 0.0));

        }

        dmat = mat;
        return dmat;
	}

	public SDCMatrix getSparseMatrix(){
		if (smat != null){
			return smat;
		}

		int m = wid*hei;

        // fill out a HashMap with nodes
        HashMap hmap = new HashMap<IndexPair, Complex>(5*m);

        IndexPair idx = new IndexPair(0,0);
        Complex val = new Complex(-1.0 / Math.pow(dx,2), 0.0);

        int cind, lind, rind, uind, dind, exind;
        
        for (int i=1; i<wid-1; i++){
            for (int j=1; j<hei-1; j++){

                cind = inds_to_global(i,j);
                lind = inds_to_global(i-1,j);
                rind = inds_to_global(i+1,j);
                uind = inds_to_global(i,j+1);
                dind = inds_to_global(i,j-1);

                // set the HashMap values
                // center
                idx = new IndexPair(cind, cind);
                hmap.put(idx, (kv.at(cind).times(-1.0)).plus(4.0/(dx*dx)));

                // left
                idx = new IndexPair(cind, lind);
                hmap.put(idx, val.copy());

                // right
                idx = new IndexPair(cind, rind);
                hmap.put(idx, val.copy());

                // up
                idx = new IndexPair(cind, uind);
                hmap.put(idx, val.copy());

                // down
                idx = new IndexPair(cind, dind);
                hmap.put(idx, val.copy());
            }
        }

        // SECOND ORDER absorbing boundaries
        double kwav = omega/c0;
        // top boundary
        for (int i=1; i<wid-1; i++){

            cind = inds_to_global(i,hei-1);
            dind = inds_to_global(i,hei-2);
            lind = inds_to_global(i-1, hei-1);
            rind = inds_to_global(i+1, hei-1);

            idx = new IndexPair(cind, cind);
            hmap.put(idx, new Complex(1.0/dx, -kwav - 1.0/(kwav*dx*dx)));
        
            idx = new IndexPair(cind, dind);
            hmap.put(idx, new Complex(-1.0/dx, 0.0));

            idx = new IndexPair(cind, lind);
            hmap.put(idx, new Complex(0.0, 1.0/(2.0*kwav*dx*dx)));

            idx = new IndexPair(cind, rind);
            hmap.put(idx, new Complex(0.0, 1.0/(2.0*kwav*dx*dx)));

        }

        // bottom boundary
        for (int i=1; i<wid-1; i++){
            cind = inds_to_global(i,0);
            uind = inds_to_global(i,1);
            lind = inds_to_global(i-1,0);
            rind = inds_to_global(i+1,0);

            idx = new IndexPair(cind, cind);
            hmap.put(idx, new Complex(-1.0/dx, -kwav - 1.0/(kwav*dx*dx)));
        
            idx = new IndexPair(cind, uind);
            hmap.put(idx, new Complex(1.0/dx, 0.0));

            idx = new IndexPair(cind, lind);
            hmap.put(idx, new Complex(0.0, 1.0/(2.0*kwav*dx*dx)));

            idx = new IndexPair(cind, rind);
            hmap.put(idx, new Complex(0.0, 1.0/(2.0*kwav*dx*dx)));

        }

        // left boundary
        for (int j=1;j<hei-1; j++){
            cind = inds_to_global(0,j);
            rind = inds_to_global(1,j);
            uind = inds_to_global(0,j+1);
            dind = inds_to_global(0,j-1);


            idx = new IndexPair(cind, cind);
            hmap.put(idx, new Complex(-1.0/dx, -kwav - 1.0/(kwav*dx*dx)));
        
            idx = new IndexPair(cind, rind);
            hmap.put(idx, new Complex(1.0/dx, 0.0));

            idx = new IndexPair(cind, uind);
            hmap.put(idx, new Complex(0.0, 1.0/(2.0*kwav*dx*dx)));

            idx = new IndexPair(cind, dind);
            hmap.put(idx, new Complex(0.0, 1.0/(2.0*kwav*dx*dx)));

        }

        // right boundary
        for (int j=1; j<hei-1; j++){
            cind = inds_to_global(wid-1,j);
            lind = inds_to_global(wid-2,j);
            uind = inds_to_global(wid-1,j+1);
            dind = inds_to_global(wid-1,j-1);

            idx = new IndexPair(cind, cind);
            hmap.put(idx, new Complex(1.0/dx, -kwav - 1.0/(kwav*dx*dx)));
        
            idx = new IndexPair(cind, lind);
            hmap.put(idx, new Complex(-1.0/dx, 0.0));

            idx = new IndexPair(cind, uind);
            hmap.put(idx, new Complex(0.0, 1.0/(2.0*kwav*dx*dx)));

            idx = new IndexPair(cind, dind);
            hmap.put(idx, new Complex(0.0, 1.0/(2.0*kwav*dx*dx)));

        }

        // corners
        // bottom left
        cind = inds_to_global(0,0);
        uind = inds_to_global(0,1);
        rind = inds_to_global(1,0);
        idx = new IndexPair(cind, cind);
        hmap.put(idx, new Complex(-2.0/dx, 0.0));
        idx = new IndexPair(cind, uind);
        hmap.put(idx, new Complex(1.0/dx, 0.0));
        idx = new IndexPair(cind, rind);
        hmap.put(idx, new Complex(1.0/dx, 0.0));

        // bottom right
        cind = inds_to_global(wid,0);
        uind = inds_to_global(wid,1);
        lind = inds_to_global(wid-1,0);
        idx = new IndexPair(cind, uind);
        hmap.put(idx, new Complex(1.0/dx, 0.0));
        idx = new IndexPair(cind, lind);
        hmap.put(idx, new Complex(-1.0/dx, 0.0));

        // top right
        cind = inds_to_global(wid,hei);
        dind = inds_to_global(wid,hei-1);
        lind = inds_to_global(wid-1,hei);
        idx = new IndexPair(cind, cind);
        hmap.put(idx, new Complex(2.0/dx, 0.0));
        idx = new IndexPair(cind, dind);
        hmap.put(idx, new Complex(-1.0/dx, 0.0));
        idx = new IndexPair(cind, lind);
        hmap.put(idx, new Complex(-1.0/dx, 0.0));

        // top left
        cind = inds_to_global(0,hei);
        dind = inds_to_global(0,hei-1);
        rind = inds_to_global(1,hei);
        idx = new IndexPair(cind, rind);
        hmap.put(idx, new Complex(1.0/dx, 0.0));
        idx = new IndexPair(cind, dind);
        hmap.put(idx, new Complex(-1.0/dx, 0.0));


        // make a matrix out of the hashmap
        smat = new SDCMatrix(m, m, hmap);

        return smat;
	}


	protected int inds_to_global(int i, int j){
		return j*hei + i;
	}
}