package helmholtz;

import java.util.HashMap;
import java.util.TreeMap;
import java.util.Map;

import java.util.Set;
import java.util.Iterator;

/*
 * Sparse Double Complex Matrix implementation
 * using CRS format
*/
public class SDCMatrix extends Object{
	private Complex[] values;
	private int[] col_ind;
	private int[] row_ptr;
	private int size;			// number of non-zero elements
	private int mrows, ncols;	// matrix dimensions

	public SDCMatrix(int m, int n, HashMap<IndexPair, Complex> map){
	
		// set metadata
		mrows = m;
		ncols = n;
		size = map.size();

		// initialize row_ptr size
		row_ptr = new int[n+1];

		// initialize values and col_ind
		values = new Complex[size];
		col_ind = new int[size];

		//System.out.println("unsorted hashmap size: "+map.size());
		
		// run through the HashMap to create CRS format
		// sort the hashmap on row index by making a TreeMap
		Map<IndexPair, Complex> sorted = new TreeMap<IndexPair, Complex>(map);
		Set set = sorted.entrySet();
		Iterator iterator = set.iterator();

		// iterate through the sorted map to set the row pointer
		Map.Entry mf = (Map.Entry)iterator.next();
		IndexPair indf = (IndexPair)mf.getKey();
		int rpos = 0;
		for (int i=0; i<mrows; i++){
			
			if (i < indf.i){
				row_ptr[i] = rpos;
			}
			else{
				if (iterator.hasNext()){
					// iterate through the set until
					// we get a new row 
					// (not necessarily consecutive)
					row_ptr[i] = rpos;
					while (indf.i == i && iterator.hasNext()){
						mf = (Map.Entry)iterator.next();
						indf = (IndexPair)mf.getKey();
						rpos++;
					}
				}
				else{
					row_ptr[i] = rpos+1;
				}
			}
		}

		// iterate through the sorted map and
		Complex val; 
		iterator = set.iterator();
		mf = (Map.Entry)iterator.next();
		indf = (IndexPair)mf.getKey();
		//row_ptr[0] = 0;
		col_ind[0] = indf.j;
		val = (Complex) mf.getValue();
		values[0] = val.copy();
		int row_pos = indf.i;
		int row_last = row_pos;
		int pos = 1;

		// System.out.println("first row: "+row_pos);
		// System.out.println("set has "+set.size()+" values");
		
		while(iterator.hasNext()) {
		  Map.Entry me = (Map.Entry)iterator.next();

		  // get key
		  IndexPair indp =(IndexPair)me.getKey();

		  // set row (if new)
		  //System.out.println("row: "+indp.i);
		  // row_pos = indp.i;
		  // if (row_pos != row_last){
		  // 	row_ptr[rpos] = pos;
		  // 	rpos++;
		  // }

		  // set col
		  col_ind[pos] = indp.j;

		  // set value
		  val =(Complex)me.getValue();
		  values[pos] = val.copy();

		  // advance positions
		  pos++;
		  row_last = row_pos;
		}
	}

	public int num_nonzero(){return size;};
	public int rows(){return mrows;};
	public int cols(){return ncols;};

	public DCVector MatVec(DCVector vec){
		DCVector result = new DCVector(vec.size());
		Complex czero = new Complex(0.0, 0.0);
		result.assign(czero);

		for (int i=0; i<ncols; i++){
			result.put(i, czero.copy());
			for (int j=row_ptr[i]; j<row_ptr[i+1]; j++){

				//if (i<100) System.out.print(" column: "+j);
				if (false){
					System.out.print(" endptr: "+row_ptr[i+1]);
					System.out.print(" currentresult: "+result.at(i).toString());
					System.out.print(" xval: "+vec.at(col_ind[j]).toString());
					System.out.println(" matrixval: "+values[j].toString());
				}

				result.put(i, result.at(i).plus(values[j].times(vec.at(col_ind[j]))).copy());
			}
		}

		return result;
	}

	public DCVector TransMatVec(DCVector vec){
		DCVector result = new DCVector(vec.size());
		Complex czero = new Complex(0.0, 0.0);
		result.assign(czero);

		for (int j=0; j<ncols; j++){
			for (int i=row_ptr[j]; i<row_ptr[j+1]; i++){
				result.put(col_ind[i], result.at(i).plus(values[i].times(vec.at(j))));
			}
		}

		return result;
	}


	public DCVector solveDiagonal(DCVector vec){
		DCVector out = new DCVector(vec.size());
		out.assign(new Complex(0.0, 0.0));

		// iterate through the matrix
		// if a value is on the diagonal, 
		// divide the corresponding vec value by it
		// and return as result. 
		for (int i=0; i<ncols; i++){
			for (int j=row_ptr[i]; j<row_ptr[i+1]; j++){

				if (col_ind[j] == i){
					out.put(i, vec.at(i).div(values[j]));
				}
			}
		}

		return out;
	}



	public void printPreview(){
		System.out.print("val: ");
		for (int i=0; i<10; i++){
			System.out.print(", "+values[i].toString());
		}
		System.out.println("");
		System.out.print("col: ");
		for (int i=0; i<10; i++){
			System.out.print(", "+col_ind[i]);
		}
		System.out.println("");
		System.out.print("row_ptr: ");
		for (int i=0; i<10; i++){
			System.out.print(", "+row_ptr[i]);
		}
		System.out.println("");
	}

}

