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

    public void print_layout(){
        int cind;
        System.out.println("FloorLayout: ");
        for (int j=fplan.get_num_length()-1; j>=0; j--){
            for (int i=0; i<fplan.get_num_width(); i++){
                cind = fplan.reg_inds_to_global(i,j);
                if (i==xloc_ind && j==yloc_ind){
                    System.out.print(-1);
                    continue;
                } 
                System.out.print(fplan.material[cind]);
                System.out.print(",");
            }
            System.out.println(" ");
        }
    }
}
