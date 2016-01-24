package helmholtz;

import java.util.Random;

public class DCVector extends Object{
	private Complex[] data;
	private int length;

	public DCVector(int len){
		length = len;
		data = new Complex[len];
	}

	public DCVector copy(){
		DCVector out = new DCVector(length);
		for (int i=0; i<length; i++){
			out.put(i, data[i].copy());
		}
		return out;
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

	public DCVector plus(DCVector vec){
		DCVector out = new DCVector(length);
		for (int i=0; i<length; i++){
			out.put(i, data[i].plus(vec.at(i)));
		}
		return out;
	}

	public DCVector minus(DCVector vec){
		DCVector out = new DCVector(length);
		for (int i=0; i<length; i++){
			out.put(i, data[i].minus(vec.at(i)));
		}
		return out;
	}

	public DCVector times(Complex val){
		DCVector out = new DCVector(length);
		for (int i=0; i<length; i++){
			out.put(i, data[i].times(val));
		}
		return out;
	}

	public DCVector div(Complex val){
		DCVector out = new DCVector(length);
		for (int i=0; i<length; i++){
			out.put(i, data[i].div(val));
		}
		return out;
	}

	public Complex dot(DCVector vec){
		Complex out = new Complex(0.0, 0.0);
		for (int i=0; i<length; i++){
			out = out.plus(data[i].times(vec.at(i).conj()));
		}
		return out;
	}

	public void print(){
		for(int i=0; i<length; i++){
			System.out.println("i: "+i+" val: "+data[i].toString());
		}
	}

	public void print2D(int mrows, int ncols){
		for (int i=0; i<mrows; i++){
			for (int j=0; j<ncols; j++){
				System.out.print(data[j*ncols+i].toString());
				System.out.print(",");
			}
			System.out.println(" ");
		}
	}

	public static DCVector rand(int m){
		Random rgen = new Random();
		DCVector out = new DCVector(m);

		for (int i=0; i<m; i++){
			out.put(i, new Complex(rgen.nextDouble(), rgen.nextDouble()));
		}
		return out;
	}

}