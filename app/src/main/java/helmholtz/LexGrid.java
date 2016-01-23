package helmholtz;

public class LexGrid extends Object{
	protected int mrows, ncols;
	protected double dx;			// grid spacing

	public LexGrid(int rows, int cols, double deltax){
		mrows = rows;
		ncols = cols;
		dx = deltax;
	}

	public LexGrid coarsen(){
		int newrows, newcols;

		newrows = (mrows-1)/2+1;
		newcols = (ncols-1)/2+1;

		LexGrid subgrid = new LexGrid(newrows, newcols, dx*2.0);
		return subgrid;
	}

	public int rows(){return mrows;};
	public int cols(){return ncols;};
	public double res(){return dx;};

	public int numpoints(){
		return mrows*ncols;
	}

	// restricts a grid of points onto a grid with
	// half as many points in each direction
	// Full weighting is used
	public DCVector restrict(DCVector longvec){

		int newrows, newcols;
		Complex c, u, d, l, r, ul, ur, dl, dr;
		Complex sum;

		newrows = (mrows-1)/2+1;
		newcols = (ncols-1)/2+1;
		DCVector newvec = longvec.copy();
		DCVector out = new DCVector(newrows*newcols);

		// apply full weighting for interior points
		// and injection for boundary points
		for (int i=1; i<mrows-1; i++){
			for (int j=1; j<ncols-1; j++){
				//System.out.println("i,j: "+i+","+j);

				c = longvec.at(j*ncols + i);
				u = longvec.at((j+1)*ncols + i);
				d = longvec.at((j-1)*ncols + i);
				l = longvec.at((j)*ncols + i-1);
				r = longvec.at((j)*ncols + i+1);
				ul = longvec.at((j+1)*ncols + i-1);
				ur = longvec.at((j+1)*ncols + i+1);
				dl = longvec.at((j-1)*ncols + i-1);
				dr = longvec.at((j-1)*ncols + i+1);


				sum = c.times(0.25);
				sum = sum.plus(u.times(0.125));
				sum = sum.plus(d.times(0.125));
				sum = sum.plus(l.times(0.125));
				sum = sum.plus(r.times(0.125));
				sum = sum.plus(ul.times(0.0625));
				sum = sum.plus(ur.times(0.0625));
				sum = sum.plus(dl.times(0.0625));
				sum = sum.plus(dr.times(0.0625));

				newvec.put(j*ncols+i, sum);
			}
		}

		int ir=0, jr=0;
		for (int i=0; i<mrows; i+=2){
			for (int j=0; j<ncols; j+=2){
				out.put(jr*newcols+ir, newvec.at(j*ncols+i));
				jr++;
			}
			jr=0;
			ir++;
		}

		return out;
	}

	// interpolates a grid of points from this grid onto 
	// a grid that has twice as many points
	public DCVector interpolate(DCVector shortvec){

		int newrows, newcols;
		Complex bl, br, tl, tr;
		Complex interp, im1, im2;

		// newrows = mrows*2-1;
		// newcols = ncols*2-1;
		newrows = 2*(mrows-1)+1;
		newcols = 2*(ncols-1)+1;
		DCVector out = new DCVector(newrows*newcols);
		out.assign(new Complex(0,0));

		// System.out.println("shortvec: "+shortvec.size());
		// System.out.println("rows: "+mrows+" cols: "+ncols);
		// System.out.println("newrows: "+newrows+" newcols: "+newcols);

		// exact interpolation at coarse grid points
		for (int i=0; i<mrows; i++){
			for (int j=0; j<ncols; j++){
				//System.out.println("loc:"+(j*ncols+i));
				// equivalent fine grid point
				out.put(2*j*newcols+2*i, shortvec.at(j*ncols+i).copy());

			}
		}

		// linear interpolation at horizontal coarse grid points
		for (int i=0; i<mrows-1; i++){
			for (int j=0; j<ncols; j++){

				// interpolate
				interp = (shortvec.at(j*ncols+i).plus(shortvec.at(j*ncols+i+1))).times(0.5);

				// equivalent fine grid point
				out.put(2*j*ncols+2*i+1, interp);

			}
		}

		// linear interpolation at vertical coarse grid points
		for (int i=0; i<mrows; i++){
			for (int j=0; j<ncols-1; j++){

				// interpolate
				interp = (shortvec.at(j*ncols+i).plus(shortvec.at((j+1)*ncols+i))).times(0.5);

				// equivalent fine grid point
				out.put((2*j+1)*ncols+2*i, interp);

			}
		}


		// bilinear interpolation for remaining points
		for (int i=1; i<newrows; i+=2){
			for (int j=1; j<newcols; j+=2){

				bl = shortvec.at((j-1)/2*ncols + (i-1)/2);
				br = shortvec.at((j-1)/2*ncols + (i+1)/2);
				tl = shortvec.at((j+1)/2*ncols + (i-1)/2);
				tr = shortvec.at((j+1)/2*ncols + (i+1)/2);


				im1 = (bl.plus(br)).times(0.5);
				im2 = (tl.plus(tr)).times(0.5);

				interp = (im1.plus(im2)).times(0.5);

				out.put(j*newcols+i, interp);
			}
		}

		return out;
	}


	// w-Jacobi relaxation as a smoother
	// u_(m+1) = 0.25*[h^2*f + uleft_m + uright_m + uup_m + udown_m]
	public DCVector smooth(DCVector rough, DCVector rhs, int ntimes){

		Complex u, d, l, r;
		Complex tmp;
		DCVector sm = rough.copy();//new DCVector(rough.size());
		DCVector prev = sm;
		//System.out.println("roughsize: "+rough.size());
		//System.out.println("rhssize: "+rhs.size());

		for (int iter=1; iter<=ntimes; iter++){

			for (int i=1; i<mrows-1; i++){
				for (int j=1; j<ncols-1; j++){
					//System.out.println("(i,j): "+i+", "+j);

					// add left, right, up, down
					u = prev.at((j+1)*ncols + i);
					d = prev.at((j-1)*ncols + i);
					l = prev.at((j)*ncols + i-1);
					r = prev.at((j)*ncols + i+1);

					tmp = u.plus(d).plus(l).plus(r);
					tmp = (tmp.plus(rhs.at(j*ncols+i).times(dx*dx))).times(0.25);

					sm.put(j*ncols+i, tmp.copy());
				}
			}

			prev = sm;
		}

		return sm;
	}

	public void print_vector(DCVector vec){

		for (int i=0; i<mrows; i++){
			for (int j=0; j<ncols; j++){
				System.out.print(vec.at(j*ncols+i));
				System.out.print(",");
			}
			System.out.println(" ");
		}
	}

}