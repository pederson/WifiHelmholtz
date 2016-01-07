package helmholtz;



public class IndexPair extends Object implements Comparable<IndexPair>{
	public int i;
	public int j;

	public IndexPair(int i_ind, int j_ind){
		i = i_ind;
		j = j_ind;
	}

	@Override
	public int compareTo(IndexPair ind){
		if (this.i < ind.i) return -1;
		else return 1;
	}
}