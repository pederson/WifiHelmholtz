package helmholtz;

// HashMap from java utils
import java.util.HashMap;

public class HelmholtzSolver {
    boolean layout_changed, source_changed;
    FloorLayout layout;
    double[] solution;

    DCVector rhs;
    SDCMatrix mat;

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

        // fill out a HashMap with nodes
        HashMap hmap = new HashMap<IndexPair, Complex>(5*m);

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
        
        for (int i=1; i<layout.fplan.num_width-1; i++){
            for (int j=1; j<layout.fplan.num_length-1; j++){

                cind = layout.fplan.reg_inds_to_global(i,j);
                lind = layout.fplan.reg_inds_to_global(i-1,j);
                rind = layout.fplan.reg_inds_to_global(i+1,j);
                uind = layout.fplan.reg_inds_to_global(i,j+1);
                dind = layout.fplan.reg_inds_to_global(i,j-1);

                nsq = Math.pow(layout.fplan.get_permittivity(i,j),2);
                sigma = layout.fplan.get_conductivity(i,j);

                // set the HashMap values
                // center
                idx = new IndexPair(cind, cind);
                val = new Complex(ksq*nsq - 4.0/Math.pow(layout.fplan.res, 2), -omega*mu0*sigma);
                hmap.put(idx, val);

                // left
                idx = new IndexPair(cind, lind);
                val = new Complex(1.0 / Math.pow(layout.fplan.res,2), 0.0);
                hmap.put(idx, val);

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

        // deal with boundary conditions
        // // top boundary
        // for (int i=0; i<layout.fplan.num_width; i++){
        //     cind = layout.fplan.reg_inds_to_global(2*i,2*(layout.fplan.num_length-1));
        //     matrix.setQuick(cind, cind,1.0);
        // }

        // // bottom boundary
        // for (int i=0; i<layout.fplan.num_width; i++){
        //     cind = layout.fplan.reg_inds_to_global(2*i,0);
        //     matrix.setQuick(cind, cind,1.0);
        // }

        // // left boundary
        // for (int j=0;j<layout.fplan.num_length; j++){
        //     cind = layout.fplan.reg_inds_to_global(0,2*j);
        //     matrix.setQuick(cind, cind, 1.0);
        // }

        // // right boundary
        // for (int j=0; j<layout.fplan.num_length; j++){
        //     cind = layout.fplan.reg_inds_to_global(2*(layout.fplan.num_width-1),2*j);
        //     matrix.setQuick(cind, cind,1.0);
        // }

        // make a matrix out of the hashmap
        mat = new SDCMatrix(m, m, hmap);

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
        double tol = 1.0e-5;  // tolerance on the residual
        int itermax = 300;  // max iteration count

        // System.out.println("About to MatVec");
        // DCVector soln = mat.MatVec(rhs);
        // //***** THIS IS NOT THE ACTUAL SOLUTION *** JUST A TEST
        // for (int i=0; i<20; i++){
        //     soln = mat.MatVec(rhs);
        //     rhs = soln;
        // }
        // //*********************************************
        // System.out.println("MatVec done");

        DCVector soln = new DCVector(rhs.size());
        BiCGSTAB slvr = new BiCGSTAB();
        slvr.set_max_iters(itermax);
        slvr.set_tolerance(tol);
        soln = slvr.solve(mat, rhs, soln);

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
            for (int i=0; i<layout.fplan.get_num_width(); i++){
                cind = layout.fplan.reg_inds_to_global(i,j);
                System.out.print((int)(9*(solution[cind]-solnmin)/(solnmax-solnmin)));
                System.out.print(",");
            }
            System.out.println(" ");
        }
        
    }


}
