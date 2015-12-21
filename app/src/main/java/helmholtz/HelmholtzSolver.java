package helmholtz;

// parallel colt
// API can be found here: http://incanter.org/docs/parallelcolt/api/
import cern.colt.matrix.tdcomplex.impl.SparseDComplexMatrix2D;
import cern.colt.matrix.tdcomplex.impl.DenseDComplexMatrix1D;
import cern.colt.matrix.tdcomplex.DComplexMatrix2D;
import cern.colt.matrix.tdcomplex.DComplexMatrix1D;
import cern.jet.math.tdcomplex.DComplexFunctions;
import cern.jet.math.tdcomplex.DComplex;
import cern.colt.list.tint.IntArrayList;
import java.util.ArrayList;

// Concurrent Hash Map from Java Utils
import java.util.concurrent.*;
//import java.util.concurrent.ConcurrrentHashMap;

public class HelmholtzSolver {

    FloorLayout layout;
    double[] solution;
    DComplexMatrix2D matrix;

    ConcurrentHashMap<Long, double[]> hmap;
    //DoubleMatrix2D matrix;

    public HelmholtzSolver(FloorLayout flayout){
        layout = flayout;
        solution = new double[layout.fplan.num_cells_total];
        
        // allocating the matrix here creates problems
        // The Colt interface hasn't yet supported constructing 
        // Sparse Complex matrices in this way.
        // Need to populate a ConcurrentHashMap instead, then construct
        // with the HashMap
        //matrix = new SparseDComplexMatrix2D(0, 0);
        //matrix = new SparseDComplexMatrix2D(layout.fplan.num_cells_total, layout.fplan.num_cells_total, 5*layout.fplan.num_cells_total, 0.3, 1.0);

   }

    public void fillMatrix(){

        int m = layout.fplan.num_cells_total;

        // fill out a ConcurrentHashMap with nodes
        System.out.println("about to instantiate HashMap");
        //hmap = new ConcurrentHashMap<Long, double[]>(5*m);
        //hmap = new ConcurrentHashMap<Long, double[]>();
        System.out.println("HashMap instantiated");
        // matrix.assign(0.0, 0.0);     // assign zeros and trim to basically "clear" the matrix
        // matrix.trimToSize();

        matrix = new SDCMatrix2D(m,m);

        double c0 = 2.99e+8;            // speed of light
        double eps0 = 8.85418782e-12;   // vacuum permittivity
        double mu0 = Math.PI*4.0e-7;    // vacuum permeability
        double ksq;
        double nsq;
        double omega;
        double sigma, wsq;

        Indices idx = new Indices(0,0);
        Long pos = new Long(0L);
        double[] val = {0.0, 0.0};

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

                // set for ParallelColt SparseDComplexMatrix2D
                matrix.setQuick(cind, cind ,ksq*nsq - 4.0/Math.pow(layout.fplan.res, 2), -omega*mu0*sigma);
                matrix.setQuick(cind, lind, 1 / Math.pow(layout.fplan.res,2), 0.0);
                matrix.setQuick(cind, rind, 1 / Math.pow(layout.fplan.res,2), 0.0);
                matrix.setQuick(cind, uind, 1 / Math.pow(layout.fplan.res,2), 0.0);
                matrix.setQuick(cind, dind, 1 / Math.pow(layout.fplan.res,2), 0.0);
            
                // set the HashMap values

                // // center
                // val[0] = ksq*nsq - 4.0/Math.pow(layout.fplan.res, 2);
                // val[1] = -omega*mu0*sigma;
                // pos = new Long((long)cind*m + cind);
                // hmap.put(pos, val);

                // // left
                // val[0] = 1.0 / Math.pow(layout.fplan.res,2);
                // val[1] = 0.0;
                // pos = new Long((long)cind*m + lind);
                // hmap.put(pos, val);

                // // right
                // pos = new Long((long)cind*m + rind);
                // hmap.put(pos, val);

                // // up
                // pos = new Long((long)cind*m + uind);
                // hmap.put(pos, val);

                // // down
                // pos = new Long((long)cind*m + dind);
                // hmap.put(pos, val);
            }
        }

        //System.out.println("HashMap size: "+hmap.size());

        // construct sparse matrix using hashmap
        System.out.println("about to allocate matrix");
        //matrix = new SDCMatrix2D(m, m, hmap, 1, 1, 1, 1);
        System.out.println("matrix is allocated");
        //matrix.setQuick(0, 0, 1.0, -1.0);

        // check nonzero entries
        IntArrayList rowList = new IntArrayList();
        IntArrayList colList = new IntArrayList();
        ArrayList<double[]> valList = new ArrayList<double[]>();
        matrix.getNonZeros(rowList, colList, valList);
        System.out.println("there are "+valList.size()+" nonzeros in the list");


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


        System.out.print("number cells total: " );
        System.out.println(layout.fplan.num_cells_total);
        System.out.print("matrix size: " );
        System.out.println(matrix.size());
        
    }


    public void solve(){

        // first, fill in the matrix if not already filled
        fillMatrix();

        // some parameters
        double tol = 1.0e-5;  // tolerance on the residual
        int itermax = 100;  // max iteration count


        // RHS
        System.out.println("Setting up b");
        DComplexMatrix1D b = new DenseDComplexMatrix1D(layout.fplan.num_cells_total);
        b.assign(0.0, 0.0); // zero out the vector
        b.setQuick(layout.fplan.reg_inds_to_global(layout.xloc_ind, layout.yloc_ind), 1.0, 0.0);

        // initial guess and residual
        System.out.println("Setting up x");
        DComplexMatrix1D x = new DenseDComplexMatrix1D(layout.fplan.num_cells_total);
        x.assign(0.0, 0.0);
        System.out.println("Setting up r");
        DComplexMatrix1D r = new DenseDComplexMatrix1D(layout.fplan.num_cells_total);
        r.assign(b);
        double[] minus_real = {-1.0, 0.0};
        double[] plus_real = {1.0, 0.0};
        System.out.print("matrix size: " );
        System.out.println(matrix.size());
        System.out.println("Calculating initial r");
        matrix.zMult(x, r, minus_real, plus_real, false);
        //r = b-matrix*x;
        double r0 = DComplex.abs(r.zDotProduct(r));
        double resid = r0;
        double[] alpha;
        double[] beta;

        // iterative GMRES


        // iterative BiCG


        // iterative CG
        System.out.println("Setting up d");
        DComplexMatrix1D d = new DenseDComplexMatrix1D(layout.fplan.num_cells_total);
        d.assign(r);
        System.out.println("matvec");
        DComplexMatrix1D matvec = new DenseDComplexMatrix1D(layout.fplan.num_cells_total);
        matrix.zMult(d, matvec);
        alpha = DComplex.plus(d.zDotProduct(r), d.zDotProduct(matvec));
        x.assign(b);
        x.assign(DComplexFunctions.mult(alpha));

        System.out.println("Starting iterations");
        int ctr=0;
        while(resid > tol && ctr < itermax)
        {
            // calculate residual
            r.assign(b);
            matrix.zMult(x, r, minus_real, plus_real, false);
            resid = DComplex.abs(r.zDotProduct(r));

            // pick new search direction
            beta = DComplex.div(r.zDotProduct(matvec), d.zDotProduct(matvec));
            beta = DComplex.mult(minus_real, beta);
            d.assign(DComplexFunctions.mult(beta));
            d.assign(r, DComplexFunctions.plus);
            //d = r - beta*d;

            // one matrix-vector product
            matrix.zMult(d, matvec);

            // calculate new step size
            alpha = DComplex.div(d.zDotProduct(r), d.zDotProduct(matvec));

            // update x
            d.assign(DComplexFunctions.mult(alpha));
            x.assign(d, DComplexFunctions.plus);
            //x += alpha*d;

            ctr++;
            System.out.println("iter: "+ctr+" resid: "+resid);
        }

        // extract the resulting magnitude
        for (int i=0; i<layout.fplan.num_cells_total; i++){
            solution[i] = DComplex.abs(x.getQuick(i));
        }


        // // set the right hand side
        // double[] rhs = new double[layout.fplan.num_cells_total];
        // for (int i=0; i<layout.fplan.num_cells_total; i++) rhs[i] = 0;
        // rhs[layout.fplan.reg_inds_to_global(layout.xloc_ind, layout.yloc_ind)] = 1;

        // double c0 = 2.99e+8; // speed of light
        // double eps0 = 8.85418782e-12; // vacuum permittivity
        // double mu0 = Math.PI*4.0e-7;    // vacuum permeability
        // double ksq;
        // double nsq;
        // double hsq;
        // double sigma, wsq, w;
        // int cind, lind, rind, uind, dind;

        // ksq = Math.pow(2.0*Math.PI*layout.wsource.freqHz/c0, 2);
        // hsq = Math.pow(layout.fplan.res,2);
        // w = 2.0*Math.PI*layout.wsource.freqHz;

        // // Gauss-Seidel SOR
        // double tol = 1.0e-5;  // tolerance
        // int itermax = 50;  // max iteration count
        // double relax = 1.0; // relaxation coefficient
        // Complex as, aw, an, ae, ap;
        // as = new Complex(1.0/hsq, 0);
        // aw = new Complex(1.0/hsq, 0);
        // an = new Complex(1.0/hsq, 0);
        // ae = new Complex(1.0/hsq, 0);
        // Complex cbuf = new Complex(0.0, 0.0);
        // Complex csrc = new Complex(1.0, 0.0);
        // Complex czero = new Complex(0.0, 0.0);
        // Complex cinit = new Complex(0.5, 0.5);

        // Complex[] soln = new Complex[layout.fplan.num_cells_total];
        // for (int i=0; i<layout.fplan.num_cells_total; i++) soln[i] = czero;
        // for (int i=1; i<layout.fplan.num_width-1; i++){
        //     for (int j=1; j<layout.fplan.num_length-1; j++){
        //         cind = layout.fplan.reg_inds_to_global(i,j);
        //         soln[cind] = cinit;
        //     }
        // }

        // int iter = 0;

        // double err;
        // double rcur, r0;
        // Complex[] resid = new Complex[layout.fplan.num_cells_total];
        // for (int i=0; i<layout.fplan.num_cells_total; i++) resid[i] = czero;
        // // calculate the initial residual r = b - A*x
        // for (int i=1; i<layout.fplan.num_width-1; i++){
        //     for (int j=1; j<layout.fplan.num_length-1; j++){

        //         cind = layout.fplan.reg_inds_to_global(i,j);

        //         if (rhs[cind] > 0){
        //             resid[cind] = czero;
        //             continue;
        //         }

        //         lind = layout.fplan.reg_inds_to_global(i-1,j);
        //         rind = layout.fplan.reg_inds_to_global(i+1,j);
        //         uind = layout.fplan.reg_inds_to_global(i,j+1);
        //         dind = layout.fplan.reg_inds_to_global(i,j-1);

        //         nsq = Math.pow(layout.fplan.get_real_refractive_index(i,j),2);
        //         sigma = layout.fplan.get_imag_refractive_index(i,j);
        //         //sigma = 0.0;

        //         // nonzero matrix values
        //         ap = new Complex(ksq*nsq - 4.0/hsq, -w*mu0*sigma);

                
        //         resid[cind] = cbuf.add(rhs[cind]).subtract(soln[cind].multiply(ap)).subtract(soln[dind].multiply(as)).subtract(soln[lind].multiply(aw)).subtract(soln[uind].multiply(an)).subtract(soln[rind].multiply(ae));
        //         // resid[cind] = rhs[cind] - (as*soln[dind]
        //         //                          + an*soln[uind]
        //         //                          + ae*soln[rind]
        //         //                          + aw*soln[lind]
        //         //                          + ap*soln[cind]);
        //     }
        // }
        // // calculate the initial norm of the residual
        // r0 = 0;
        // for (int i=0; i<layout.fplan.num_cells_total; i++){
        //     r0 += resid[i].abs()*resid[i].abs();
        // }
        // // calculate the initial error err = norm(r)/r0
        // err = r0;
        
        // //****** THE MAIN LOOP ******************
        // while (err > tol && iter < itermax){

        //     // increment
        //     iter++;

        //     // iterate on the values
        //     for (int j=1; j<layout.fplan.num_length-1; j++){
        //         for (int i=1; i<layout.fplan.num_width-1; i++){

        //             cind = layout.fplan.reg_inds_to_global(i,j);

        //             if (rhs[cind] > 0){
        //                 //System.out.println("SOURCE");
        //                 soln[cind] = csrc;
        //                 continue;
        //             } 

        //             lind = layout.fplan.reg_inds_to_global(i-1,j);
        //             rind = layout.fplan.reg_inds_to_global(i+1,j);
        //             uind = layout.fplan.reg_inds_to_global(i,j+1);
        //             dind = layout.fplan.reg_inds_to_global(i,j-1);

        //             nsq = Math.pow(layout.fplan.get_real_refractive_index(i,j),2);
        //             sigma = layout.fplan.get_imag_refractive_index(i,j);
        //             //sigma = 0.0;

        //             // nonzero matrix values
        //             ap = new Complex(ksq*nsq - 4.0/hsq, -w*mu0*sigma);
                    
        //             if (iter==1 && cind == 20*layout.fplan.num_length+1){
        //             //if (rhs[uind] > 0){
        //             //if (cind == 481) {
        //                 System.out.println("******************** CIND(" + Double.toString(cind) + ") **********************");
        //                 System.out.print("freq (Hz): ");
        //                 System.out.println(layout.wsource.freqHz);
        //                 System.out.print("phi(P): ");
        //                 System.out.print(soln[cind].toString());
        //                 System.out.println(" A(P): " + ap.toString());
        //                 System.out.print("phi(S): " + soln[dind].toString());
        //                 System.out.println(" A(S): " + as.toString());
        //                 System.out.print("phi(N): " + soln[uind].toString());
        //                 System.out.println(" A(N): " + an.toString());
        //                 System.out.print("phi(E): " + soln[rind].toString());
        //                 System.out.println(" A(E): " + ae.toString());
        //                 System.out.print("phi(W): " + soln[lind].toString());
        //                 System.out.println(" A(W): " + aw.toString());
        //                 System.out.println("ksq: " + Double.toString(ksq));
        //                 System.out.println("nsq: " + Double.toString(nsq));
        //                 System.out.println("hsq: " + Double.toString(hsq));
        //                 System.out.println("sigma: " + Double.toString(sigma));
        //                 System.out.println("w: " + Double.toString(w));
        //                 System.out.print(" phi(n+1): ");
        //                 System.out.println(cbuf.add(rhs[cind]).subtract(soln[dind].multiply(as)).subtract(soln[lind].multiply(aw)).subtract(soln[uind].multiply(an)).subtract(soln[rind].multiply(ae)).divide(ap).multiply(relax).add(soln[cind].multiply(1-relax)).toString());

        //                 System.out.println("*********************************************************");
        //             }

        //             // update
        //             soln[cind] = cbuf.add(rhs[cind]).subtract(soln[dind].multiply(as)).subtract(soln[lind].multiply(aw)).subtract(soln[uind].multiply(an)).subtract(soln[rind].multiply(ae)).divide(ap).multiply(relax).add(soln[cind].multiply(1-relax));

        //         }
        //     }

        //     // calculate the residual r = b - A*x
        //     for (int i=1; i<layout.fplan.num_width-1; i++){
        //         for (int j=1; j<layout.fplan.num_length-1; j++){

        //             cind = layout.fplan.reg_inds_to_global(i,j);

        //             if (rhs[cind] > 0){
        //                 resid[cind] = czero;
        //                 continue;
        //             }

        //             lind = layout.fplan.reg_inds_to_global(i-1,j);
        //             rind = layout.fplan.reg_inds_to_global(i+1,j);
        //             uind = layout.fplan.reg_inds_to_global(i,j+1);
        //             dind = layout.fplan.reg_inds_to_global(i,j-1);

        //             nsq = Math.pow(layout.fplan.get_real_refractive_index(i,j),2);
        //             sigma = layout.fplan.get_imag_refractive_index(i,j);
        //             //sigma = 0.0;

        //             // nonzero matrix values
        //             ap = new Complex(ksq*nsq - 4.0/hsq, -w*mu0*sigma);

        //             resid[cind] = cbuf.add(rhs[cind]).subtract(soln[cind].multiply(ap)).subtract(soln[dind].multiply(as)).subtract(soln[lind].multiply(aw)).subtract(soln[uind].multiply(an)).subtract(soln[rind].multiply(ae));

        //             System.out.println("residual: " + Double.toString(resid));
           
        //         }
        //     }

        //     // calculate the norm of the residual
        //     rcur = 0;
        //     for (int i=0; i<layout.fplan.num_cells_total; i++){
        //         rcur += resid[i].abs()*resid[i].abs();
        //     }

        //     // calculate the error err = norm(r)/r0
        //     err = rcur/r0;

        //     System.out.print("Iteration: ");
        //     System.out.print(iter);
        //     System.out.print(" Error: ");
        //     System.out.println(err);

        // }
        // //************** END MAIN LOOP **********************

        // // calculate the resulting magnitude
        // for (int i=0; i<layout.fplan.num_cells_total; i++){
        //     solution[i] = soln[i].abs();
        // }


    }

    // absolute value of complex number
    public double cabs(double[] z){
        return Math.sqrt(z[0]*z[0] + z[1]*z[1]);
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
