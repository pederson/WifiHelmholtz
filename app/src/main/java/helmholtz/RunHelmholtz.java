package helmholtz;

//import helmholtz.FloorPlan;

public class RunHelmholtz{


	public static void main(String args[]){
		//FloorPlan fplan = new FloorPlan(2.0, 2.0, 0.00625, "TestPlan");
		FloorPlan fplan = new FloorPlan(1.0, 1.0, 0.01, "TestPlan");
		FloorLayout flayout = new FloorLayout();
		int concretelayers = 5;
		// fill in concrete on the borders
		// top and bottom ----> 5 layers
		for (int i=0; i<fplan.get_num_width(); i++){
			for (int jlow=0; jlow < concretelayers; jlow++){
				fplan.set_pixel_material(1, i, jlow);
			}
			for (int jhi=fplan.get_num_length()-1; jhi>fplan.get_num_length()-1 - concretelayers; jhi--){
				fplan.set_pixel_material(1, i, jhi);
			}

		}
		// left and right ----> 5 layers
		for (int j=0; j<fplan.get_num_length(); j++){
			for (int ilow=0; ilow < concretelayers; ilow++){
				fplan.set_pixel_material(1, ilow, j);
			}
			for (int ihi=fplan.get_num_width()-1; ihi>fplan.get_num_width()-1 - concretelayers; ihi--){
				fplan.set_pixel_material(1, ihi, j);
			}

		}

		// add a "doorway" 3/4 of the way to the top
		int j34 = (int)(0.75*(double)fplan.get_num_length());
		for (int i=0; i<(int)(0.75*(double)fplan.get_num_width()); i++){
			fplan.set_pixel_material(1, i, j34);
			fplan.set_pixel_material(1, i, j34+1);
			fplan.set_pixel_material(1, i, j34-1);
			fplan.set_pixel_material(1, i, j34+2);
			fplan.set_pixel_material(1, i, j34-2);
			fplan.set_pixel_material(1, i, j34+3);
			fplan.set_pixel_material(1, i, j34-3);
		}

		// print it
		fplan.print_floorplan();

		// add a source in middle of unit
		WifiSource source = new WifiSource(false, false);
		flayout.set_source(source);
		flayout.set_source_location(fplan.get_num_width()/2, fplan.get_num_length()/2);

		flayout.set_floorplan(fplan);

		flayout.print_layout();

		HelmholtzSolver hsolve = new HelmholtzSolver(flayout);
		hsolve.solve();
		hsolve.print_solution();
	}
}