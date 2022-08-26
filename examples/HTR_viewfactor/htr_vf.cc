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
// const double xc = 0.0, yc = 0.0, zc = 0.0;
// const double xa = 0.0, ya = 0.0, za = 1.0;
// const double rc = 0.9;

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

	//-Custom output routines.
	// Do a custom output routine to store packing statistics
	con.print_custom(
		// %i: ID
		// %q: The position vector of the particle, short for “%x %y %z”.
		// %v: The volume of the Voronoi cell.
		"ID=%i, pos=(%x,%y,%z), volume=%v", "packing_statistics.dat");

	// con.print_custom("%i %q %v", "packing_statistics_MATLAB.dat");

	// Output the particle positions and Voronoi cells in Gnuplot format
	// con.draw_particles("cylinder_p.gnu");
	// con.draw_cells_gnuplot("cylinder_v.gnu");

	// Output the particle positions and Voronoi cells in POV-Ray format
	// con.draw_particles_pov("cylinder_p.pov");
	// con.draw_cells_pov("cylinder_v.pov");

	// File Handler ==================================================
	FILE *f1=safe_fopen("viewfactor_result.dat","w");
	// Headers
	fprintf(f1, "id     x      y      z      id2    face_a\n");
	// write temp cell vertices and face vertices
	FILE *f2=safe_fopen("temp_voronoi_vertices.dat","w");
	FILE *f3=safe_fopen("temp_voronoi_facevtid.dat","w");

	int pid;
	double x,y,z;
	c_loop_all cl(con);
	voronoicell_neighbor c;
	std::vector<int> neigh, fvt, fo;
	std::vector<double> vt;
	std::vector<double> face_areas, face_normals;

	if(cl.start()) do if(con.compute_cell(c,cl)) {
		// Get particle position & pid
		cl.pos(x,y,z); pid = cl.pid();

		// Gather information about the computed Voronoi cell
		c.neighbors(neigh);
		c.face_areas(face_areas);
		c.normals(face_normals);
		c.vertices(x,y,z,vt);
		c.face_vertices(fvt);
		c.face_orders(fo);

		// Print out the information about neighbors
		printf("Particle %d: \n", pid);
		printf("Number of neighbours: %lu\n", neigh.size());
		printf("Number of faces: %d\n", c.number_of_faces());
		puts("");

		// Face Area
		// printf("face_areas[0]: %f\n", face_areas[0]);
		printf("Output Face areas:\n");
		c.output_face_areas();
		puts("");
		puts("");

		// Vertices
		printf("Number of vertices: %lu\n", vt.size());
		printf("Output Vertices:\n");
		c.output_vertices();
		puts("");
		puts("");
		// Output Vertices to F2
		for (unsigned int tmp_i = 0; tmp_i < vt.size()/3; tmp_i++) {
			fprintf(f2, "%d %.18f %.18f %.18f\n", tmp_i, vt[tmp_i*3], vt[tmp_i*3+1], vt[tmp_i*3+2]);
		}
		fprintf(f2, "%d %.18f %.18f %.18f\n", pid, x, y, z);

		// Face Vertices
		printf("Number of face vertices size: %lu\n", fvt.size());
		c.output_face_vertices();
		puts("");
		puts("");
		// Output Face Vertices to F3
		int tmp_count = 0;
		int num_vertices = 0;
		for (unsigned int tmp_i = 0; tmp_i < fvt.size(); tmp_i++) {
			if (tmp_count == 0) num_vertices = fvt[tmp_i];
			if (tmp_count < num_vertices) {
				fprintf(f3, "%d ", fvt[tmp_i]);
				tmp_count++;
			} else {
				fprintf(f3, "%d\n", fvt[tmp_i]);
				tmp_count = 0;
			}
		}

		// Face Orders
		printf("Number of face orders: %lu\n", fo.size());
		c.output_face_orders();
		puts("");
		puts("");

		printf("====================================================================\n");
		// 获取表面面积向量

		int tmp = 0; // index for Face Vertices
		// 循环每个面
		for(unsigned int j=0;j<neigh.size();j++) {
			// printf("%d\n", fo[j] - fvt[tmp]);
			// Face Index
			int ii, jj, kk;
			ii = j*3;
			jj = j*3+1;
			kk = j*3+2;

			// Calculate Y L H;
			double Y, H, V;
			printf("vertices of current face: %d\n", fo[j]);

			// Calculate plane equation (face) ================================================================
			// params: A B C D
			printf("Getting plane using Vertice: %d, %d, %d\n", fvt[tmp+2], fvt[tmp+3], fvt[tmp+4]);
			double x1 = vt[fvt[tmp+2]*3];
			double y1 = vt[fvt[tmp+2]*3+1];
			double z1 = vt[fvt[tmp+2]*3+2];
			double x2 = vt[fvt[tmp+3]*3];
			double y2 = vt[fvt[tmp+3]*3+1];
			double z2 = vt[fvt[tmp+3]*3+2];
			double x3 = vt[fvt[tmp+4]*3];
			double y3 = vt[fvt[tmp+4]*3+1];
			double z3 = vt[fvt[tmp+4]*3+2];
			// Plane params
			double A = ( (y2 - y1) * (z3 - z1) - (z2 - z1) * (y3 - y1) );
			double B = ( (z2 - z1) * (x3 - x1) - (x2 - x1) * (z3 - z1) );
			double C = ( (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1) );
			double D = ( 0 - (A * x1 + B * y1 + C * z1) );

			// check error with internal face normals
			double err_a = fabs(A / sqrt(A*A+B*B+C*C) + face_normals[ii]);
			double err_b = fabs(B / sqrt(A*A+B*B+C*C) + face_normals[jj]);
			double err_c = fabs(C / sqrt(A*A+B*B+C*C) + face_normals[kk]);
			if (fmax(fabs(err_a), fmax(fabs(err_b), fabs(err_c))) > 1e-10) {
				printf("Normal Vector (inter): (%f, %f, %f)\n", face_normals[ii], face_normals[jj], face_normals[kk]);
				printf("Normal Vector (pt123): (%f, %f, %f)\n", A / sqrt(A*A+B*B+C*C),
												  B / sqrt(A*A+B*B+C*C),
												  C / sqrt(A*A+B*B+C*C));
				printf("Normal Vector diff   : (%g, %g, %g)\n", err_a, err_b, err_c);
				printf("Normal Vector diff is too large!\n");
				exit(1);
			} else {
				printf("Normal Vector of face correct!\n");
			}

			H = fabs(A * x + B * y + C * z + D) / sqrt(A * A + B * B + C * C);
			// printf("H: %f\n", H);

			// Calculate projection point of sphere center ====================================================
			// params: p_x, p_y, p_z
			// double p_x = face_normals[ii]*H + x;
			// double p_y = face_normals[jj]*H + y;
			// double p_z = face_normals[kk]*H + z;
			double p_x = -A / sqrt(A*A+B*B+C*C)*H + x;
			double p_y = -B / sqrt(A*A+B*B+C*C)*H + y;
			double p_z = -C / sqrt(A*A+B*B+C*C)*H + z;

			// Check if on the plane
			double testp = A*p_x + B*p_y + C*p_z + D;
			if (fabs(testp) > 1e-15) {
				p_x = A / sqrt(A*A+B*B+C*C)*H + x;
				p_y = B / sqrt(A*A+B*B+C*C)*H + y;
				p_z = C / sqrt(A*A+B*B+C*C)*H + z;
				testp = A*p_x + B*p_y + C*p_z + D;
			}
			printf("Projection point of sphere center: (%f, %f, %f)\n", p_x, p_y, p_z);
			if (fabs(testp) > 1e-15) {
				printf("check if projection pt on plane: %g\n", testp);
				printf("Projection pt is not on the plane!\n");
				exit(1);
			} else {
				printf("Projection point correct!\n");
			}

			// double test1 = A*vt[fvt[tmp+4]*3] + B*vt[fvt[tmp+4]*3+1] + C*vt[fvt[tmp+4]*3+2] + D;
			double test1 = A*vt[fvt[tmp+1]*3] + B*vt[fvt[tmp+1]*3+1] + C*vt[fvt[tmp+1]*3+2] + D;
			if (fabs(test1) > 1e-15) {
				printf("check if First point on plane : %g\n", test1);
				printf("First point is not on the plane!\n");
				exit(1);
			} else {
				printf("First point correct!\n");
			}


			// 测试：计算三角形面积之和(face area) ================================================================
			// params: total_s_all
			double tmp_s_all = 0.0;
			for (unsigned int l=1; l<fo[j]-1; l++) {
				int pt_1 = 1;
				int pt_2 = l+1;
				int pt_3 = l+2;
				printf("calculating triangle: %d, %d, %d\n", fvt[tmp+pt_1], fvt[tmp+pt_2], fvt[tmp+pt_3]);
				double tmp_a = sqrt(pow(vt[fvt[tmp+pt_1]*3+0] - vt[fvt[tmp+pt_2]*3+0], 2) +
								    pow(vt[fvt[tmp+pt_1]*3+1] - vt[fvt[tmp+pt_2]*3+1], 2) +
								    pow(vt[fvt[tmp+pt_1]*3+2] - vt[fvt[tmp+pt_2]*3+2], 2));
				double tmp_b = sqrt(pow(vt[fvt[tmp+pt_2]*3+0] - vt[fvt[tmp+pt_3]*3+0], 2) +
								    pow(vt[fvt[tmp+pt_2]*3+1] - vt[fvt[tmp+pt_3]*3+1], 2) +
								    pow(vt[fvt[tmp+pt_2]*3+2] - vt[fvt[tmp+pt_3]*3+2], 2));
				double tmp_c = sqrt(pow(vt[fvt[tmp+pt_3]*3+0] - vt[fvt[tmp+pt_1]*3+0], 2) +
								    pow(vt[fvt[tmp+pt_3]*3+1] - vt[fvt[tmp+pt_1]*3+1], 2) +
								    pow(vt[fvt[tmp+pt_3]*3+2] - vt[fvt[tmp+pt_1]*3+2], 2));
				double tmp_p = 0.5*(tmp_a+tmp_b+tmp_c);
				double tmp_s = sqrt(tmp_p*(tmp_p-tmp_a)*(tmp_p-tmp_b)*(tmp_p-tmp_c));
				tmp_s_all += tmp_s;
			}
			// printf("total face area by triangles: %f\n", tmp_s_all);
			// printf("face area (internal)        : %f\n", face_areas[j]);
			if (fabs(tmp_s_all - face_areas[j]) > 1e-15) {
				printf("face area error: %f\n", tmp_s_all - face_areas[j]);
			} else {
				printf("Face area correct!\n");
			}


			puts("");


			// 循环每个顶点 =========================================================================================
			// params: a, b, c, p, S, total_s
			printf("iterating over vertices:\n");
			int k;
			double total_s = 0;
			for(k=1;k<=fo[j];k++) {
				int iii, jjj, kkk, iiin, jjjn, kkkn;
				iii = fvt[k+tmp]*3;
				jjj = fvt[k+tmp]*3+1;
				kkk = fvt[k+tmp]*3+2;
				int k_next;
				if (k == fo[j]) k_next = 1; else k_next = k+1;
				iiin = fvt[k_next+tmp]*3;
				jjjn = fvt[k_next+tmp]*3+1;
				kkkn = fvt[k_next+tmp]*3+2;
				// printf("calculate area of triangle: p, %d, %d\n", fvt[k+tmp], fvt[k_next+tmp]);

				// length of three edges
				double tmp_a, tmp_b, tmp_c;
				tmp_a = sqrt(pow(p_x - vt[iii],  2)+pow(p_y - vt[jjj],  2)+pow(p_z - vt[kkk],  2));
				tmp_b = sqrt(pow(p_x - vt[iiin], 2)+pow(p_y - vt[jjjn], 2)+pow(p_z - vt[kkkn], 2));
				tmp_c = sqrt(pow(vt[iii] - vt[iiin], 2) +pow(vt[jjj] - vt[jjjn], 2) + pow(vt[kkk] - vt[kkkn], 2));
				double tmp_p = 0.5*(tmp_a+tmp_b+tmp_c);
				double tmp_s = sqrt(tmp_p*(tmp_p-tmp_a)*(tmp_p-tmp_b)*(tmp_p-tmp_c));
				printf("area, %f\n", tmp_s);
				total_s += tmp_s;

				Y = tmp_s*2.0/tmp_c;

				double B1, B3;
				B1 = Y/H;
				B3 = tmp_a/Y;
				V = 0.25/pi*acos(1/B3)-0.125/pi*asin( ((1-pow(B1,2))*pow(B3,2)-2) / ((1+pow(B1,2))*pow(B3,2)) );
			}
			tmp += k;
			// printf("tmp var: %d\n", tmp);
			printf("face area (internal): %f\n", face_areas[j]);
			printf("face area (me)      : %f\n", total_s);
			printf("-------------------------------------");
			puts("");

			// if (j == 2) break;
		}

		break;
	} while (cl.inc());

	fclose(f1);
	fclose(f2);
	fclose(f3);
}
