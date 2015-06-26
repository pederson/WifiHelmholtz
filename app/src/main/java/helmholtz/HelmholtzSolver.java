package helmholtz;

// normal colt
// import cern.colt.matrix.DoubleMatrix1D;
// import cern.colt.matrix.DoubleMatrix2D;
// import cern.colt.matrix.impl.SparseDoubleMatrix1D;
// import cern.colt.matrix.linalg.Algebra;
// import cern.colt.matrix.impl.SparseDoubleMatrix2D;

// parallel colt
// API can be found here: http://incanter.org/docs/parallelcolt/api/
// import cern.colt.matrix.tdouble.DoubleMatrix1D;
// import cern.colt.matrix.tdouble.DoubleMatrix2D;
// import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;
// import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix1D;
// import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
// import cern.colt.matrix.tdouble.algo.solver.DoubleGMRES;
// import cern.colt.matrix.tdouble.algo.solver.DoubleBiCG;
// import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

import org.apache.commons.math3.complex.Complex;

public class HelmholtzSolver {

    FloorLayout layout;
    double[] solution;
    //DoubleMatrix2D matrix;


    public HelmholtzSolver(FloorLayout flayout){
        layout = flayout;
        solution = new double[layout.fplan.num_cells_total];
        //matrix = new SparseDoubleMatrix2D(layout.fplan.num_cells_total, layout.fplan.num_cells_total);
        //matrix = new SparseDoubleMatrix2D(layout.fplan.num_cells_total*2, layout.fplan.num_cells_total*2);
    }

    public void fillMatrix(){

        /*
        matrix.assign(0.0);     // assign zeros and trim to basically "clear" the matrix
        matrix.trimToSize();

        double c0 = 2.99e+8; // speed of light
        double eps0 = 8.85418782e-12; // vacuum permittivity
        double ksq;
        double nsq;
        double sigma, wsq;
        int cind, lind, rind, uind, dind, exind;

        ksq = Math.pow(2.0*Math.PI*layout.wsource.freqHz/c0, 2);


        for (int i=1; i<layout.fplan.num_width-1; i++){
            for (int j=1; j<layout.fplan.num_length-1; j++){

                nsq = Math.pow(layout.fplan.get_real_refractive_index(i,j),2);
                sigma = layout.fplan.get_imag_refractive_index(i,j);

                // real equation
                cind = layout.fplan.reg_inds_to_global(2*i,2*j);
                lind = layout.fplan.reg_inds_to_global(2*(i-1),2*j);
                rind = layout.fplan.reg_inds_to_global(2*(i+1),2*j);
                uind = layout.fplan.reg_inds_to_global(2*i,2*(j+1));
                dind = layout.fplan.reg_inds_to_global(2*i,2*(j-1));
                exind = layout.fplan.reg_inds_to_global(2*i, 2*j)+1;

                matrix.setQuick(cind, cind ,ksq*nsq - 4.0/Math.pow(layout.fplan.res, 2));
                matrix.setQuick(cind, lind, 1 / Math.pow(layout.fplan.res,2));
                matrix.setQuick(cind, rind, 1 / Math.pow(layout.fplan.res,2));
                matrix.setQuick(cind, uind, 1 / Math.pow(layout.fplan.res,2));
                matrix.setQuick(cind, dind, 1 / Math.pow(layout.fplan.res,2));
                //matrix.setQuick(cind, exind, 2.0*Math.PI*layout.wsource.freqHz*sigma/eps0);
                
                // // imaginary equation
                // cind = layout.fplan.reg_inds_to_global(2*i,2*j)+1;
                // lind = layout.fplan.reg_inds_to_global(2*(i-1),2*j)+1;
                // rind = layout.fplan.reg_inds_to_global(2*(i+1),2*j)+1;
                // uind = layout.fplan.reg_inds_to_global(2*i,2*(j+1))+1;
                // dind = layout.fplan.reg_inds_to_global(2*i,2*(j-1))+1;
                // exind = layout.fplan.reg_inds_to_global(2*i, 2*j);

                // matrix.setQuick(cind, cind ,ksq*nsq - 4.0/Math.pow(layout.fplan.res, 2));
                // matrix.setQuick(cind, lind, 1 / Math.pow(layout.fplan.res,2));
                // matrix.setQuick(cind, rind, 1 / Math.pow(layout.fplan.res,2));
                // matrix.setQuick(cind, uind, 1 / Math.pow(layout.fplan.res,2));
                // matrix.setQuick(cind, dind, 1 / Math.pow(layout.fplan.res,2));
                // //matrix.setQuick(cind, exind, -2.0*Math.PI*layout.wsource.freqHz*sigma/eps0);

                // //way that works for sure
                // cind = layout.fplan.reg_inds_to_global(i,j);
                // lind = layout.fplan.reg_inds_to_global(i-1,j);
                // rind = layout.fplan.reg_inds_to_global(i+1,j);
                // uind = layout.fplan.reg_inds_to_global(i,j+1);
                // dind = layout.fplan.reg_inds_to_global(i,j-1);

                // nsq = Math.pow(layout.fplan.get_real_refractive_index(i,j),2);
                // ksq = Math.pow(2.0*Math.PI*layout.wsource.freqHz/c0, 2);
                // wsq = Math.pow(2.0*Math.PI*layout.wsource.freqHz, 2);
                // matrix.setQuick(cind, cind ,ksq*nsq - 4.0/Math.pow(layout.fplan.res, 2));
                // matrix.setQuick(cind, lind, 1 / Math.pow(layout.fplan.res,2));
                // matrix.setQuick(cind, rind, 1 / Math.pow(layout.fplan.res,2));
                // matrix.setQuick(cind, uind, 1 / Math.pow(layout.fplan.res,2));
                // matrix.setQuick(cind, dind, 1 / Math.pow(layout.fplan.res,2));
            }
        }

        System.out.println("The matrix is filled!");

        // // top boundary
        // for (int i=0; i<layout.fplan.num_width; i++){
        //     cind = layout.fplan.reg_inds_to_global(i,layout.fplan.num_length-1);
        //     matrix.setQuick(cind, cind,1.0);
        // }

        // // bottom boundary
        // for (int i=0; i<layout.fplan.num_width; i++){
        //     cind = layout.fplan.reg_inds_to_global(i,0);
        //     matrix.setQuick(cind, cind,1.0);
        // }

        // // left boundary
        // for (int j=0;j<layout.fplan.num_length; j++){
        //     cind = layout.fplan.reg_inds_to_global(0,j);
        //     matrix.setQuick(cind, cind, 1.0);
        // }

        // // right boundary
        // for (int j=0; j<layout.fplan.num_length; j++){
        //     cind = layout.fplan.reg_inds_to_global(layout.fplan.num_width-1,j);
        //     matrix.setQuick(cind, cind,1.0);
        // }

        // top boundary
        for (int i=0; i<layout.fplan.num_width; i++){
            cind = layout.fplan.reg_inds_to_global(2*i,2*(layout.fplan.num_length-1));
            matrix.setQuick(cind, cind,1.0);
        }

        // bottom boundary
        for (int i=0; i<layout.fplan.num_width; i++){
            cind = layout.fplan.reg_inds_to_global(2*i,0);
            matrix.setQuick(cind, cind,1.0);
        }

        // left boundary
        for (int j=0;j<layout.fplan.num_length; j++){
            cind = layout.fplan.reg_inds_to_global(0,2*j);
            matrix.setQuick(cind, cind, 1.0);
        }

        // right boundary
        for (int j=0; j<layout.fplan.num_length; j++){
            cind = layout.fplan.reg_inds_to_global(2*(layout.fplan.num_width-1),2*j);
            matrix.setQuick(cind, cind,1.0);
        }
        */
    }


    public void solve(){
        //fillMatrix();

        // DoubleMatrix1D soln = new DenseDoubleMatrix1D(layout.fplan.num_cells_total*2);
        // DoubleMatrix1D rhs = new SparseDoubleMatrix1D(layout.fplan.num_cells_total*2);
        // rhs.setQuick(layout.fplan.reg_inds_to_global(layout.xloc_ind*2, layout.yloc_ind*2), 1);

        // DoubleMatrix1D soln = new DenseDoubleMatrix1D(layout.fplan.num_cells_total);
        // DoubleMatrix1D rhs = new SparseDoubleMatrix1D(layout.fplan.num_cells_total);
        // rhs.setQuick(layout.fplan.reg_inds_to_global(layout.xloc_ind, layout.yloc_ind), 1);

        // direct inverse calculation (long)
        // Algebra linalg = new Algebra();
        // DoubleMatrix2D matinv = linalg.inverse(matrix);
        // DoubleMatrix1D soln = linalg.mult(matinv, rhs);
        

        // iterative GMRES
        // DoubleGMRES gmrsolver = new DoubleGMRES(soln);
        // try{
        //     soln = gmrsolver.solve(matrix, rhs, soln);
        // }
        // catch (IterativeSolverDoubleNotConvergedException e){
        //     System.out.println("There was an issue in GMRES!");
        //     return;
        // }

        // // iterative BiCG
        // DoubleBiCG bicgsolver = new DoubleBiCG(soln);
        // try{
        //     soln = bicgsolver.solve(matrix, rhs, soln);
        // }
        // catch (IterativeSolverDoubleNotConvergedException e){
        //     System.out.println("There was an issue in BiCG!");
        //     System.out.println(e.getMessage());
        //     return;
        // }

        // set the right hand side
        double[] rhs = new double[layout.fplan.num_cells_total];
        for (int i=0; i<layout.fplan.num_cells_total; i++) rhs[i] = 0;
        rhs[layout.fplan.reg_inds_to_global(layout.xloc_ind, layout.yloc_ind)] = 1;

        double c0 = 2.99e+8; // speed of light
        double eps0 = 8.85418782e-12; // vacuum permittivity
        double ksq;
        double nsq;
        double hsq;
        double sigma, wsq, w;
        int cind, lind, rind, uind, dind;

        ksq = Math.pow(2.0*Math.PI*layout.wsource.freqHz/c0, 2);
        hsq = Math.pow(layout.fplan.res,2);
        w = 2.0*Math.PI*layout.wsource.freqHz;

        // Gauss-Seidel SOR
        double tol = 0.01; 
        int itermax = 100;
        Complex as, aw, an, ae, ap;
        as = new Complex(1.0, 0);
        aw = new Complex(1.0, 0);
        an = new Complex(1.0, 0);
        ae = new Complex(1.0, 0);
        Complex[] soln = new Complex[layout.fplan.num_cells_total];
        Complex czero = new Complex(0.0, 0.0);
        for (int i=0; i<layout.fplan.num_cells_total; i++) soln[i] = czero;
        Complex cinit = new Complex(0.5, 0.5);
        for (int i=1; i<layout.fplan.num_width-1; i++){
            for (int j=1; j<layout.fplan.num_length-1; j++){
                cind = layout.fplan.reg_inds_to_global(i,j);

                soln[cind] = cinit;
            }
        }
        int iter = 0;
        double err;
        double rcur, r0;
        Complex[] resid = new Complex[layout.fplan.num_cells_total];
        for (int i=0; i<layout.fplan.num_cells_total; i++) resid[i] = czero;
        // calculate the norm of the residual
        r0 = 0;
        for (int i=0; i<layout.fplan.num_cells_total; i++){
            r0 += soln[i].abs()*soln[i].abs();
        }
        // calculate the error err = norm(r)/r0
        err = r0;
        Complex cbuf = new Complex(0.0, 0.0);
        while (err > tol && iter < itermax){

            // iterate on the values
            for (int i=1; i<layout.fplan.num_width-1; i++){
                for (int j=1; j<layout.fplan.num_length-1; j++){

                    cind = layout.fplan.reg_inds_to_global(i,j);
                    lind = layout.fplan.reg_inds_to_global(i-1,j);
                    rind = layout.fplan.reg_inds_to_global(i+1,j);
                    uind = layout.fplan.reg_inds_to_global(i,j+1);
                    dind = layout.fplan.reg_inds_to_global(i,j-1);

                    nsq = Math.pow(layout.fplan.get_real_refractive_index(i,j),2);
                    sigma = layout.fplan.get_imag_refractive_index(i,j);
                    // nonzero matrix values
                    ap = new Complex(hsq*ksq*nsq - 4, hsq*w*sigma/eps0);
                    

                    // update
                    //cbuf.add(hsq*rhs[cind]);
                    soln[cind] = cbuf.add(hsq*rhs[cind]).subtract(soln[dind].multiply(as)).subtract(soln[lind].multiply(aw)).subtract(soln[uind].multiply(an)).subtract(soln[rind].multiply(ae)).divide(ap);
                    // soln[cind] = (hsq*rhs[cind] 
                    //             - soln[dind].multiply(as) 
                    //             - aw*soln[lind] 
                    //             - an*soln[uind] 
                    //             - ae*soln[rind])/ap;
                }
            }

            // calculate the residual r = b - A*x
            for (int i=1; i<layout.fplan.num_width-1; i++){
                for (int j=1; j<layout.fplan.num_length-1; j++){

                    cind = layout.fplan.reg_inds_to_global(i,j);
                    lind = layout.fplan.reg_inds_to_global(i-1,j);
                    rind = layout.fplan.reg_inds_to_global(i+1,j);
                    uind = layout.fplan.reg_inds_to_global(i,j+1);
                    dind = layout.fplan.reg_inds_to_global(i,j-1);

                    nsq = Math.pow(layout.fplan.get_real_refractive_index(i,j),2);
                    sigma = layout.fplan.get_imag_refractive_index(i,j);
                    // nonzero matrix values
                    ap = new Complex(hsq*ksq*nsq - 4, hsq*w*sigma/eps0);

                    resid[cind] = cbuf.add(hsq*rhs[cind]).subtract(soln[cind].multiply(ap)).subtract(soln[dind].multiply(as)).subtract(soln[lind].multiply(aw)).subtract(soln[uind].multiply(an)).subtract(soln[rind].multiply(ae));
                    // resid[cind] = rhs[cind] - (as*soln[dind]
                    //                          + an*soln[uind]
                    //                          + ae*soln[rind]
                    //                          + aw*soln[lind]
                    //                          + ap*soln[cind]);
                }
            }

            // calculate the norm of the residual
            rcur = 0;
            for (int i=0; i<layout.fplan.num_cells_total; i++){
                rcur += soln[i].abs()*soln[i].abs();
            }

            // calculate the error err = norm(r)/r0
            err = rcur/r0;

            System.out.print("Iteration: ");
            System.out.print(iter);
            System.out.print(" Error: ");
            System.out.println(err);

            // increment
            iter++;

        }

        // for (int i=0; i<layout.fplan.num_cells_total; i++){
        //     solution[i] = Math.sqrt(soln.get(2*i)*soln.get(2*i) + soln.get(2*i+1)*soln.get(2*i+1));
        // }
        // for (int i=0; i<layout.fplan.num_cells_total; i++){
        //     solution[i] = soln.get(i);
        // }
        for (int i=0; i<layout.fplan.num_cells_total; i++){
            solution[i] = soln[i].abs();
        }
    }

    public void print_solution(){
        int cind;
        System.out.println("Solution: ");

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
