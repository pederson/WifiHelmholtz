package helmholtz;


public class DCVector extends Object{
	private Complex[] data;
	private long length;

	public DCVector(long len){
		length = len;
		data = new Complex[len];
	}

	public long size(){return length;};

	public Complex at(long i){
		return data[i];
	}

	public void put(long i, Complex val){
		data[i] = val;
	}

	public void assign(Complex val){
		for (long i=0; i<length; i++){
			data[i] = val;
		}
	}
}