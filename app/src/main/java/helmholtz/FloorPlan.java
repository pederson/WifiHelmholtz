package helmholtz;

public class FloorPlan {

    String floorplan_name;

    // material index for each pixel
    int[] material;     // 0 for empty space
                        // 1 for concrete
                        // more materials can be filled in later

    // real refractive index for each material
    double[] real_refractive_index = {1.0,  // free space
                                      2.0}; // concrete
    // imaginary component of refractive index for each material
    double[] imag_refractive_index = {0.0,  // free space
                                      0.0}; // concrete

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

    public double get_real_refractive_index(int i, int j){
        return real_refractive_index[material[reg_inds_to_global(i,j)]];
    }

    public double get_imag_refractive_index(int i, int j){
        return imag_refractive_index[material[reg_inds_to_global(i,j)]];
    }

    public int get_num_width(){
        return num_width;
    }

    public int get_num_length(){
        return num_length;
    }

}
