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
	private long[] col_ind;
	private long[] row_ptr;
	private long size;			// number of non-zero elements
	private long mrows, ncols;	// matrix dimensions

	public SDCMatrix(long m, long n, HashMap<IndexPair, Complex> map){
	
		// set metadata
		mrows = m;
		ncols = n;
		size = map.size();

		// initialize row_ptr size
		row_ptr = new long[n];

		// initialize values and col_ind
		values = new Complex[size];
		col_ind = new long[size];

		// run through the HashMap to create CRS format
		// sort the hashmap on row index by making a TreeMap
		Map<IndexPair, Complex> sorted = new TreeMap<IndexPair, Complex>(map);

		// iterate through the sorted map and 
		Set set = sorted.entrySet();
		Iterator iterator = set.iterator();
		Map.Entry mf = (Map.Entry)iterator.next
		IndexPair indf = me.getKey();
		row_ptr[0] = 0;
		col_ind[0] = indf.j;
		long row_last = indf.i;

		long pos = 1;
		long rpos = 1;
		while(iterator.hasNext()) {
		  Map.Entry me = (Map.Entry)iterator.next();

		  // get key
		  IndexPair indp = me.getKey();

		  // set row (if new)
		  row_pos = indp.i;
		  if (row_pos != row_last){
		  	row_ptr[rpos] = pos;
		  	rpos++;
		  }

		  // set col
		  col_ind[pos] = indp.j;

		  // set value
		  values[pos] = me.getValue();

		  // advance positions
		  pos++;
		  row_last = row_pos;
		}
	}

	public long num_nonzero(){return size;};
	public long rows(){return mrows;};
	public long cols(){return ncols;};

	public DCVector MatVec(DCVector vec){
		DCVector result = new DCVector(vec.size());

		for (long i=0; i<ncols; i++){
			result.put(i, Complex(0.0, 0.0));
			for (long j=row_ptr[i]; j<row_ptr[i+1]; j++){
				result.put(i, result.plus(values[j].times(vec.at(col_ind[j]))));
			}
		}

		return result;
	}

	public DCVetor TransMatVec(DCVector vec){
		DCVector result = new DCVector(vec.size());

		for (long i=0; i<ncols; i++){
			result.put(i, Complex(0.0, 0.0));
		}
		for (long j=0; j<ncols; j++){
			for (long i=row_ptr[j]; i<row_ptr[j+1]; i++){
				result.put(col_ind[i], result.plus(values[i].times(vec.at(j))));
			}
		}

		return result;
	}

}

