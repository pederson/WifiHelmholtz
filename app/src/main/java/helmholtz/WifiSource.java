package helmholtz;

public class WifiSource {

    boolean double_antenna;
    int double_antenna_orientation;            // 0 = horizontal, 1 = 45 degrees CCW

    double freqHz;          // either 2.4 GHz or 5.0 GHz

    public WifiSource(boolean freq_high, boolean antenna_double){
        if (freq_high){
            freqHz = 5.0e+9;
        }
        else freqHz = 2.4e+9;

        if (antenna_double){
            double_antenna = true;
        }
        else{
            double_antenna = false;
        }

        double_antenna_orientation = 0;
    }

}
