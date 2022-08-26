// Cylindrical wall example code
// together with cone zone...
// Also output some useful pack statistics.
//
// Author   : Josh
// Date     : February 3th 2022
//

#include "voro++.hh"
using namespace voro;

// Set up constants for the container geometry
const double x_min=-0.9,x_max=0.9;
const double y_min=-0.9,y_max=0.9;
const double z_min=0;
// const double z_max=0.4;
const double z_max=2.24;

const double pi=3.1415926535897932384626433832795;

// Set the computational grid size
const int n_x=23,n_y=23,n_z=25;
const double xc = 0.0, yc = 0.0, zc = 0.0;
const double xa = 0.0, ya = 0.0, za = 1.0;
const double rc = 0.9;

int main() {
	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block.
	container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);

	// Add a cylindrical wall to the container (wall.hh)
	// center     direction   r
	// (xc,yc,zc) (xa,ya,za) (radius)
	double cone_h = (0.9-0.25)/sqrt(3);
	wall_cylinder cyl(0.0, 0.0, cone_h, 0.0, 0.0, 1.0, 0.9);
	con.add_wall(cyl);

	double cone_z = -0.25/sqrt(3);
	double ang1 = 60.0/180.0*pi;
	wall_cone cone(0, 0, cone_z, 0, 0, 1, ang1);
	con.add_wall(cone);

	// Import the particles from a file
	con.import("pebble_pos.dat");
	// con.import("pebble_pos_cone.dat");

	//-Custom output routines.
	// Do a custom output routine to store packing statistics
	con.print_custom(
		// %i: ID
		// %q: The position vector of the particle, short for “%x %y %z”.
		// %v: The volume of the Voronoi cell.
		"ID=%i, pos=(%x,%y,%z), volume=%v", "packing_statistics.dat");

	con.print_custom("%i %q %v", "packing_statistics_MATLAB.dat");

	//-Output the particle positions and Voronoi cells in Gnuplot format
	con.draw_particles("cylinder_p.gnu");
	con.draw_cells_gnuplot("cylinder_v.gnu");

	// Output the particle positions and Voronoi cells in POV-Ray format
	con.draw_particles_pov("cylinder_p.pov");
	con.draw_cells_pov("cylinder_v.pov");

	// Compute the volume of the Voronoi cells and compare it to the
    // exact frustum volume
    // evol=pi*1*(0.5*0.5+0.5*1+1*1)/3;
    // vvol=con.sum_cell_volumes();
    // printf("Exact frustum volume : %g\n"
    //        "Voronoi cell volume  : %g\n"
    //        "Difference           : %g\n",evol,vvol,vvol-evol);

	int i;
	double x,y,z,r;
	voronoicell_neighbor c;
	std::vector<int> neigh;
	std::vector<double> vd;

	FILE *f1=safe_fopen("result.dat","w");

	c_loop_all cl(con);
	// if(cl.start()) do if(con.compute_cell(c,cl)) {
	// 	cl.pos(i,x,y,z,r);
	// 	c.neighbors(neigh);

	// 	// Get face areas
	// 	c.face_areas(vd);

	// 	// Print out the information about neighbors
	// 	if (neigh.size() != vd.size()) fprintf(f1, "err, ");
	// 	fprintf(f1, "Particle %d at (% 1.3f,% 1.3f,% 1.3f):",i,x,y,z);
	// 	for(unsigned int j=0;j<neigh.size();j++) fprintf(f1, " |%d %f|",neigh[j], vd[j]);
	// 	// puts("");
	// 	fprintf(f1, "\n");
	// } while (cl.inc());

	// Headers
	fprintf(f1, "id     x      y      z      id2    face_a\n");

	if(cl.start()) do if(con.compute_cell(c,cl)) {
		// Get particle position & pid
		cl.pos(i,x,y,z,r);
		// Get neighbor
		c.neighbors(neigh);
		// Get face areas
		c.face_areas(vd);

		// Print out the information about neighbors
		if (neigh.size() != vd.size()) printf("err, neighbor size != face areas size\n");
		for(unsigned int j=0;j<neigh.size();j++) fprintf(f1, "%d %6f %6f %6f %d %6f\n", i, x, y, z, neigh[j], vd[j]);
		fprintf(f1, "\n");
	} while (cl.inc());

	fclose(f1);
}
