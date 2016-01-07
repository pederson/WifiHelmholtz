package helmholtz;



public class IndexPair extends Object implements Comparable{
	public long i;
	public long j;

	public IndexPair(long i_ind, long j_ind){
		i = i_ind;
		j = j_ind;
	}

	public int compareTo(IndexPair ind){
		if (this.i < ind.i) return -1;
		else if (this.i > ind.i) return 1;
		else return 0;
	}
}