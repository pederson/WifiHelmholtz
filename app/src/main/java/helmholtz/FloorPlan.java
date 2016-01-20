package helmholtz;

public class FloorPlan {

    String floorplan_name;

    // material index for each pixel
    int[] material;     // 0 for empty space
                        // 1 for concrete
                        // 2 - brick
                        // 3 - wood
                        // 4 - glass
                        // 5 - drywall
                        // more materials can be added

    // see the report: "Building Materials and Propagation: Final Report"
    // by Ofcom for source of these values
    // everything assumes 5 GHz (OK b/c property variations b/w 5 and 2.5 GHz are small)
    // relative permittivity
    double[] permittivity = {1.0,   // free space
                             5.31,  // concrete
                             3.75,  // brick
                             1.99,  // wood
                             6.27,  // glass
                             2.94   // drywall
                             }; 

    // conductivity
    double[] conductivity = {0.0,       // free space
                             1.1996e-1, // concrete
                             3.8e-2,    // brick
                             2.638e-2,  // wood
                             2.93e-2,   // glass
                             3.62e-2    // drywall
                             };

    double length;      // total length (meters)
    double width;       // total width (meters)
    double res;         // size of one pixel (meters)
    int num_width;      // number of pixels in the x direction
    int num_length;     // number of pixels in the y direction
    int num_cells_total;    // total number of pixels

    public FloorPlan(double swidth, double slength, double sres, String fname){
        floorplan_name = fname;
        width = swidth;
        length = slength;
        res = sres;
        num_width = (int)(width/res);
        num_length = (int)(length/res);

        if (num_width%2 == 0) num_width++;
        if (num_length%2 == 0) num_length++;
        num_cells_total = num_width*num_length;

        material = new int[num_cells_total];
    }


    int reg_inds_to_global(int i, int j){
        return j*num_width + i;
    }

    public void set_pixel_material(int material_ind, int i, int j){
        material[reg_inds_to_global(i,j)] = material_ind;
    }

    public int get_pixel_material(int i, int j){
        return material[reg_inds_to_global(i,j)];
    }

    public double get_permittivity(int i, int j){
        return permittivity[material[reg_inds_to_global(i,j)]];
    }

    public double get_conductivity(int i, int j){
        return conductivity[material[reg_inds_to_global(i,j)]];
    }

    public int get_num_width(){
        return num_width;
    }

    public int get_num_length(){
        return num_length;
    }

    public void print_floorplan(){
        int cind;
        System.out.println("FloorPlan: ");
        for (int j=num_length-1; j>=0; j--){
            for (int i=0; i<num_width; i++){
                cind = reg_inds_to_global(i,j);
                System.out.print(material[cind]);
                System.out.print(",");
            }
            System.out.println(" ");
        }
        
    }

}
