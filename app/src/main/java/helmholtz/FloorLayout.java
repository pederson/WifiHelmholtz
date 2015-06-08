package helmholtz;

public class FloorLayout {

    FloorPlan fplan;
    WifiSource wsource;
    int xloc_ind, yloc_ind;

    public void set_floorplan(FloorPlan inplan){
        fplan = inplan;
    }

    public void set_source(WifiSource insource){
        wsource = insource;
    }

    public void set_source_location(int i, int j){
        xloc_ind = i;
        yloc_ind = j;
    }
}
