package helmholtz;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix1D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;

public class HelmholtzSolver {

    FloorLayout layout;
    double[] solution;
    DoubleMatrix2D matrix;


    public HelmholtzSolver(FloorLayout flayout){
        layout = flayout;
        solution = new double[layout.fplan.num_cells_total];
        matrix = new SparseDoubleMatrix2D(layout.fplan.num_cells_total, layout.fplan.num_cells_total);
    }

    public void fillMatrix(){

        matrix.assign(0.0);     // assign zeros and trim to basically "clear" the matrix
        matrix.trimToSize();

        double c0 = 2.99e+8; // speed of light
        double ksq;
        double nsq;
        int cind, lind, rind, uind, dind;

        for (int i=1; i<layout.fplan.num_width-1; i++){
            for (int j=1; j<layout.fplan.num_length-1; j++){

                cind = layout.fplan.reg_inds_to_global(i,j);
                lind = layout.fplan.reg_inds_to_global(i-1,j);
                rind = layout.fplan.reg_inds_to_global(i+1,j);
                uind = layout.fplan.reg_inds_to_global(i,j+1);
                dind = layout.fplan.reg_inds_to_global(i,j-1);

                nsq = Math.pow(layout.fplan.get_real_refractive_index(i,j),2) - Math.pow(layout.fplan.get_imag_refractive_index(i,j),2);
                ksq = Math.pow(2.0*Math.PI*layout.wsource.freqHz/c0, 2);
                matrix.setQuick(cind, cind ,ksq/nsq - 4.0/Math.pow(layout.fplan.res, 2));
                matrix.setQuick(cind, lind, 1 / Math.pow(layout.fplan.res,2));
                matrix.setQuick(cind, rind, 1 / Math.pow(layout.fplan.res,2));
                matrix.setQuick(cind, uind, 1 / Math.pow(layout.fplan.res,2));
                matrix.setQuick(cind, dind, 1 / Math.pow(layout.fplan.res,2));
            }
        }

        // top boundary
        for (int i=0; i<layout.fplan.num_width; i++){
            cind = layout.fplan.reg_inds_to_global(i,layout.fplan.num_length-1);
            matrix.setQuick(cind, cind,1.0);
        }

        // bottom boundary
        for (int i=0; i<layout.fplan.num_width; i++){
            cind = layout.fplan.reg_inds_to_global(i,0);
            matrix.setQuick(cind, cind,1.0);
        }

        // left boundary
        for (int j=0;j<layout.fplan.num_length; j++){
            cind = layout.fplan.reg_inds_to_global(0,j);
            matrix.setQuick(cind, cind, 1.0);
        }

        // right boundary
        for (int j=0; j<layout.fplan.num_length; j++){
            cind = layout.fplan.reg_inds_to_global(layout.fplan.num_width-1,j);
            matrix.setQuick(cind, cind,1.0);
        }
    }


    public void solve(){
        //DoubleMatrix1D soln = new DenseDoubleMatrix1D(layout.fplan.num_cells_total);
        DoubleMatrix1D rhs = new SparseDoubleMatrix1D(layout.fplan.num_cells_total);

        rhs.setQuick(layout.fplan.reg_inds_to_global(layout.xloc_ind, layout.yloc_ind), 1);

        Algebra linalg = new Algebra();
        DoubleMatrix2D matinv = linalg.inverse(matrix);
        DoubleMatrix1D soln = linalg.mult(matinv, rhs);

        for (int i=0; i<layout.fplan.num_cells_total; i++){
            solution[i] = soln.get(i);
        }
    }


}
