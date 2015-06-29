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


        // iterative GMRES


        // iterative BiCG


        // set the right hand side
        double[] rhs = new double[layout.fplan.num_cells_total];
        for (int i=0; i<layout.fplan.num_cells_total; i++) rhs[i] = 0;
        rhs[layout.fplan.reg_inds_to_global(layout.xloc_ind, layout.yloc_ind)] = 1;

        double c0 = 2.99e+8; // speed of light
        double eps0 = 8.85418782e-12; // vacuum permittivity
        double mu0 = Math.PI*4.0e-7;    // vacuum permeability
        double ksq;
        double nsq;
        double hsq;
        double sigma, wsq, w;
        int cind, lind, rind, uind, dind;

        ksq = Math.pow(2.0*Math.PI*layout.wsource.freqHz/c0, 2);
        hsq = Math.pow(layout.fplan.res,2);
        w = 2.0*Math.PI*layout.wsource.freqHz;

        // Gauss-Seidel SOR
        double tol = 1.0e-5;  // tolerance
        int itermax = 50;  // max iteration count
        double relax = 1.0; // relaxation coefficient
        Complex as, aw, an, ae, ap;
        as = new Complex(1.0/hsq, 0);
        aw = new Complex(1.0/hsq, 0);
        an = new Complex(1.0/hsq, 0);
        ae = new Complex(1.0/hsq, 0);
        Complex cbuf = new Complex(0.0, 0.0);
        Complex csrc = new Complex(1.0, 0.0);
        Complex czero = new Complex(0.0, 0.0);
        Complex cinit = new Complex(0.5, 0.5);

        Complex[] soln = new Complex[layout.fplan.num_cells_total];
        for (int i=0; i<layout.fplan.num_cells_total; i++) soln[i] = czero;
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
        // calculate the initial residual r = b - A*x
        for (int i=1; i<layout.fplan.num_width-1; i++){
            for (int j=1; j<layout.fplan.num_length-1; j++){

                cind = layout.fplan.reg_inds_to_global(i,j);

                if (rhs[cind] > 0){
                    resid[cind] = czero;
                    continue;
                }

                lind = layout.fplan.reg_inds_to_global(i-1,j);
                rind = layout.fplan.reg_inds_to_global(i+1,j);
                uind = layout.fplan.reg_inds_to_global(i,j+1);
                dind = layout.fplan.reg_inds_to_global(i,j-1);

                nsq = Math.pow(layout.fplan.get_real_refractive_index(i,j),2);
                sigma = layout.fplan.get_imag_refractive_index(i,j);
                //sigma = 0.0;

                // nonzero matrix values
                ap = new Complex(ksq*nsq - 4.0/hsq, -w*mu0*sigma);

                
                resid[cind] = cbuf.add(rhs[cind]).subtract(soln[cind].multiply(ap)).subtract(soln[dind].multiply(as)).subtract(soln[lind].multiply(aw)).subtract(soln[uind].multiply(an)).subtract(soln[rind].multiply(ae));
                // resid[cind] = rhs[cind] - (as*soln[dind]
                //                          + an*soln[uind]
                //                          + ae*soln[rind]
                //                          + aw*soln[lind]
                //                          + ap*soln[cind]);
            }
        }
        // calculate the initial norm of the residual
        r0 = 0;
        for (int i=0; i<layout.fplan.num_cells_total; i++){
            r0 += resid[i].abs()*resid[i].abs();
        }
        // calculate the initial error err = norm(r)/r0
        err = r0;
        
        //****** THE MAIN LOOP ******************
        while (err > tol && iter < itermax){

            // increment
            iter++;

            // iterate on the values
            for (int j=1; j<layout.fplan.num_length-1; j++){
                for (int i=1; i<layout.fplan.num_width-1; i++){

                    cind = layout.fplan.reg_inds_to_global(i,j);

                    if (rhs[cind] > 0){
                        //System.out.println("SOURCE");
                        soln[cind] = csrc;
                        continue;
                    } 

                    lind = layout.fplan.reg_inds_to_global(i-1,j);
                    rind = layout.fplan.reg_inds_to_global(i+1,j);
                    uind = layout.fplan.reg_inds_to_global(i,j+1);
                    dind = layout.fplan.reg_inds_to_global(i,j-1);

                    nsq = Math.pow(layout.fplan.get_real_refractive_index(i,j),2);
                    sigma = layout.fplan.get_imag_refractive_index(i,j);
                    //sigma = 0.0;

                    // nonzero matrix values
                    ap = new Complex(ksq*nsq - 4.0/hsq, -w*mu0*sigma);
                    
                    if (iter==1 && cind == 20*layout.fplan.num_length+1){
                    //if (rhs[uind] > 0){
                    //if (cind == 481) {
                        System.out.println("******************** CIND(" + Double.toString(cind) + ") **********************");
                        System.out.print("freq (Hz): ");
                        System.out.println(layout.wsource.freqHz);
                        System.out.print("phi(P): ");
                        System.out.print(soln[cind].toString());
                        System.out.println(" A(P): " + ap.toString());
                        System.out.print("phi(S): " + soln[dind].toString());
                        System.out.println(" A(S): " + as.toString());
                        System.out.print("phi(N): " + soln[uind].toString());
                        System.out.println(" A(N): " + an.toString());
                        System.out.print("phi(E): " + soln[rind].toString());
                        System.out.println(" A(E): " + ae.toString());
                        System.out.print("phi(W): " + soln[lind].toString());
                        System.out.println(" A(W): " + aw.toString());
                        System.out.println("ksq: " + Double.toString(ksq));
                        System.out.println("nsq: " + Double.toString(nsq));
                        System.out.println("hsq: " + Double.toString(hsq));
                        System.out.println("sigma: " + Double.toString(sigma));
                        System.out.println("w: " + Double.toString(w));
                        System.out.print(" phi(n+1): ");
                        System.out.println(cbuf.add(rhs[cind]).subtract(soln[dind].multiply(as)).subtract(soln[lind].multiply(aw)).subtract(soln[uind].multiply(an)).subtract(soln[rind].multiply(ae)).divide(ap).multiply(relax).add(soln[cind].multiply(1-relax)).toString());

                        System.out.println("*********************************************************");
                    }

                    // update
                    soln[cind] = cbuf.add(rhs[cind]).subtract(soln[dind].multiply(as)).subtract(soln[lind].multiply(aw)).subtract(soln[uind].multiply(an)).subtract(soln[rind].multiply(ae)).divide(ap).multiply(relax).add(soln[cind].multiply(1-relax));

                }
            }

            // calculate the residual r = b - A*x
            for (int i=1; i<layout.fplan.num_width-1; i++){
                for (int j=1; j<layout.fplan.num_length-1; j++){

                    cind = layout.fplan.reg_inds_to_global(i,j);

                    if (rhs[cind] > 0){
                        resid[cind] = czero;
                        continue;
                    }

                    lind = layout.fplan.reg_inds_to_global(i-1,j);
                    rind = layout.fplan.reg_inds_to_global(i+1,j);
                    uind = layout.fplan.reg_inds_to_global(i,j+1);
                    dind = layout.fplan.reg_inds_to_global(i,j-1);

                    nsq = Math.pow(layout.fplan.get_real_refractive_index(i,j),2);
                    sigma = layout.fplan.get_imag_refractive_index(i,j);
                    //sigma = 0.0;

                    // nonzero matrix values
                    ap = new Complex(ksq*nsq - 4.0/hsq, -w*mu0*sigma);

                    resid[cind] = cbuf.add(rhs[cind]).subtract(soln[cind].multiply(ap)).subtract(soln[dind].multiply(as)).subtract(soln[lind].multiply(aw)).subtract(soln[uind].multiply(an)).subtract(soln[rind].multiply(ae));

                }
            }

            // calculate the norm of the residual
            rcur = 0;
            for (int i=0; i<layout.fplan.num_cells_total; i++){
                rcur += resid[i].abs()*resid[i].abs();
            }

            // calculate the error err = norm(r)/r0
            err = rcur/r0;

            System.out.print("Iteration: ");
            System.out.print(iter);
            System.out.print(" Error: ");
            System.out.println(err);

        }
        //************** END MAIN LOOP **********************

        // calculate the resulting magnitude
        for (int i=0; i<layout.fplan.num_cells_total; i++){
            solution[i] = soln[i].abs();
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
