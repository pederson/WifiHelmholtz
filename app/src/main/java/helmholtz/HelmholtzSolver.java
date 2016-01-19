package helmholtz;

// HashMap from java utils
import java.util.HashMap;

public class HelmholtzSolver {
    boolean layout_changed, source_changed;
    FloorLayout layout;
    double[] solution;

    DCVector rhs;
    SDCMatrix mat;
    HelmholtzOperator2D hOp;

    public HelmholtzSolver(FloorLayout flayout){
        layout = flayout;
        solution = new double[layout.fplan.num_cells_total];
        
        layout_changed = true;
        source_changed = true;

    }

    private void fillRHS(){
        if (!source_changed) return;

        // RHS
        rhs = new DCVector(layout.fplan.num_cells_total);
        rhs.assign(new Complex(0.0, 0.0));
        rhs.put(layout.fplan.reg_inds_to_global(layout.xloc_ind, layout.yloc_ind), new Complex(1.0, 0.0));
    
        source_changed = false;
    }

    private void fillMatrix(){
        if (!layout_changed) return;

        int m = layout.fplan.num_cells_total;

        // grid 
        LexGrid grid = new LexGrid(layout.fplan.num_width, layout.fplan.num_length, layout.fplan.res);

        // fill out the k-squared vector
        DCVector kv = new DCVector(m);

        double c0 = 2.99e+8;            // speed of light
        double eps0 = 8.85418782e-12;   // vacuum permittivity
        double mu0 = Math.PI*4.0e-7;    // vacuum permeability
        double ksq;
        double nsq;
        double omega;
        double sigma, wsq;

        IndexPair idx = new IndexPair(0,0);
        Complex val;

        int cind, lind, rind, uind, dind, exind;

        omega = 2.0*Math.PI*layout.wsource.freqHz;
        ksq = Math.pow(2.0*Math.PI*layout.wsource.freqHz/c0, 2);
        wsq = Math.pow(omega, 2);
        
        for (int i=0; i<layout.fplan.num_width; i++){
            for (int j=0; j<layout.fplan.num_length; j++){

                cind = layout.fplan.reg_inds_to_global(i,j);

                nsq = Math.pow(layout.fplan.get_permittivity(i,j),2);
                sigma = layout.fplan.get_conductivity(i,j);

                kv.put(cind, new Complex(ksq*nsq, -omega*mu0*sigma));

            }
        }

        // create helmholtz operator
        hOp = new HelmholtzOperator2D(grid, layout.wsource.freqHz, kv);

        // get sparse matrix
        mat = hOp.getSparseMatrix();

        // System.out.println("matrix dims: "+mat.rows()+" x "+mat.cols());
        // System.out.print("matrix nonzeros: " );
        // System.out.println(mat.num_nonzero());
        
        layout_changed = false;
    }


    public void solve(){
        if (!layout_changed && !source_changed) return;

        // first, fill in the matrix if not already filled
        fillMatrix();

        // fill RHS if not already done
        fillRHS();

        // some parameters
        double tol = 1.0e-4;  // tolerance on the residual
        int itermax = 1000;  // max iteration count

        // heavily damped operator PC
        // DCVector kvec = hOp.getKsq();
        // DCVector knew = kvec.times(new Complex(1.0, -0.5));
        // HelmholtzOperator2D pcOp = new HelmholtzOperator2D(hOp.getGrid(), layout.wsource.freqHz, knew);
        // Preconditioner pc = new MGPreconditioner(pcOp);
        Preconditioner pc = new JacobiPreconditioner(mat);
        DCVector soln = new DCVector(rhs.size());
        BiCGSTAB slvr = new BiCGSTAB();
        slvr.set_max_iters(itermax);
        slvr.set_tolerance(tol);
        soln = slvr.solve(pc, mat, rhs, soln);

        // calculate the resulting magnitude
        for (int i=0; i<layout.fplan.num_cells_total; i++){
            solution[i] = soln.at(i).mod();
        }



    }


    public void print_solution(){
        int cind;
        //System.out.println("Solution: ");

        double solnmax, solnmin;
        solnmax = solution[0]; solnmin = solution[0];
        for (int i=1; i<layout.fplan.get_num_length()*layout.fplan.get_num_width(); i++){
            if (solution[i] > solnmax) solnmax = solution[i];
            if (solution[i] < solnmin) solnmin = solution[i];
        }

        for (int j=layout.fplan.get_num_length()-1; j>=0; j--){
            for (int i=0; i<layout.fplan.get_num_width()-1; i++){
                cind = layout.fplan.reg_inds_to_global(i,j);
                System.out.print((int)(255*(solution[cind]-solnmin)/(solnmax-solnmin)));
                System.out.print(",");
            }
            cind = layout.fplan.reg_inds_to_global(layout.fplan.get_num_width()-1,j);
            System.out.println((int)(255*(solution[cind]-solnmin)/(solnmax-solnmin)));
        }
        
    }


}
