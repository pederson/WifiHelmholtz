package helmholtz;


public class DCVector extends Object{
	private Complex[] data;
	private int length;

	public DCVector(int len){
		length = len;
		data = new Complex[len];
	}

	public int size(){return length;};

	public Complex at(int i){
		return data[i];
	}

	public void put(int i, Complex val){
		data[i] = val;
	}

	public void assign(Complex val){
		for (int i=0; i<length; i++){
			data[i] = val.copy();
		}
	}
}