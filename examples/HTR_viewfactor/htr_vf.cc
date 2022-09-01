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

	// Output the particle positions and Voronoi cells in Gnuplot format
	// con.draw_particles("cylinder_p.gnu");
	// con.draw_cells_gnuplot("cylinder_v.gnu");

	// Output the particle positions and Voronoi cells in POV-Ray format
	// con.draw_particles_pov("cylinder_p.pov");
	// con.draw_cells_pov("cylinder_v.pov");

	// File Handler ==================================================
	FILE *f1=safe_fopen("viewfactor_result.dat","w");
	// Headers
	fprintf(f1, "id1     id2     d                     VF1                   Ai                    VF2\n");
	FILE *f4=safe_fopen("view_factor_peb.dat", "w");
	// write temp cell vertices and face vertices
	FILE *f2=safe_fopen("temp_voronoi_vertices.dat","w");
	FILE *f3=safe_fopen("temp_voronoi_facevtid.dat","w");

	int pid;
	double x,y,z;
	c_loop_all cl(con);
	c_loop_all cl2(con);
	voronoicell_neighbor c;
	voronoicell c2;
	std::vector<int> neigh, fvt, fo;
	std::vector<double> vt;
	std::vector<double> face_areas, face_normals;

	printf("%d\n", con.total_particles());

	int index = 0;

	if(cl.start()) do if(con.compute_cell(c,cl)) {
		index++;
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
		if (index % 100 == 0) {
			printf("Particle %d - %d: ", index, pid);
			printf("Number of faces: %d\n", c.number_of_faces());
		}

		// Face Area
		// printf("Output Face areas:\n");
		// c.output_face_areas();
		// puts("");
		// puts("");

		// Vertices
		// printf("Number of vertices: %lu\n", vt.size());
		// printf("Output Vertices:\n");
		// c.output_vertices();
		// puts("");
		// puts("");
		// Output Vertices to F2 FILE
		// for (unsigned int tmp_i = 0; tmp_i < vt.size()/3; tmp_i++) {
		// 	fprintf(f2, "%d %.18f %.18f %.18f\n", tmp_i, vt[tmp_i*3], vt[tmp_i*3+1], vt[tmp_i*3+2]);
		// }
		// fprintf(f2, "%d %.18f %.18f %.18f\n", pid, x, y, z);

		// Face Vertices
		// printf("Number of face vertices size: %lu\n", fvt.size());
		// c.output_face_vertices();
		// puts("");
		// puts("");
		// Output Face Vertices to F3 FILE
		// int tmp_count = 0;
		// int num_vertices = 0;
		// for (unsigned int tmp_i = 0; tmp_i < fvt.size(); tmp_i++) {
		// 	if (tmp_count == 0) num_vertices = fvt[tmp_i];
		// 	if (tmp_count < num_vertices) {
		// 		fprintf(f3, "%d ", fvt[tmp_i]);
		// 		tmp_count++;
		// 	} else {
		// 		fprintf(f3, "%d\n", fvt[tmp_i]);
		// 		tmp_count = 0;
		// 	}
		// }

		// Face Orders
		// printf("Number of face orders: %lu\n", fo.size());
		// c.output_face_orders();
		// puts("");
		// puts("");

		// 获取表面面积向量

		double V1_all = 0, V2_all = 0;
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
			double Y, H, V1 = 0, V2 = 0;
			// printf("vertices of current face: %d\n", fo[j]);

			// Calculate plane equation (face) ================================================================
			// params: A B C D
			// printf("Getting plane using Vertice: %d, %d, %d\n", fvt[tmp+2], fvt[tmp+3], fvt[tmp+4]);
			int id_f1, id_f2, id_f3;
			if (fo[j] > 3)
			{
				id_f1 = tmp + 2;
				id_f2 = tmp + 3;
				id_f3 = tmp + 4;
			}
			else if (fo[j] == 3)
			{
				id_f1 = tmp + 1;
				id_f2 = tmp + 2;
				id_f3 = tmp + 3;
			}
			else
			{
				printf("vertices number of face is error!\n");
				exit(1);
			}

			// if no normal vec
			if (fmax(fmax(fabs(face_normals[ii]), fabs(face_normals[jj])), fabs(face_normals[kk])) < 1e-5) {
				tmp += fo[j]+1;
				continue;
			}
			double x1 = vt[fvt[id_f1]*3];
			double y1 = vt[fvt[id_f1]*3+1];
			double z1 = vt[fvt[id_f1]*3+2];
			double x2 = vt[fvt[id_f2]*3];
			double y2 = vt[fvt[id_f2]*3+1];
			double z2 = vt[fvt[id_f2]*3+2];
			double x3 = vt[fvt[id_f3]*3];
			double y3 = vt[fvt[id_f3]*3+1];
			double z3 = vt[fvt[id_f3]*3+2];
			// Plane params
			double A = ( (y2 - y1) * (z3 - z1) - (z2 - z1) * (y3 - y1) );
			double B = ( (z2 - z1) * (x3 - x1) - (x2 - x1) * (z3 - z1) );
			double C = ( (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1) );
			double D = ( 0 - (A * x1 + B * y1 + C * z1) );

			// check error with internal face normals
			double err_a0 = fabs(A / sqrt(A*A+B*B+C*C) + face_normals[ii]);
			double err_a1 = fabs(A / sqrt(A*A+B*B+C*C) - face_normals[ii]);
			double err_b0 = fabs(B / sqrt(A*A+B*B+C*C) + face_normals[jj]);
			double err_b1 = fabs(B / sqrt(A*A+B*B+C*C) - face_normals[jj]);
			double err_c0 = fabs(C / sqrt(A*A+B*B+C*C) + face_normals[kk]);
			double err_c1 = fabs(C / sqrt(A*A+B*B+C*C) - face_normals[kk]);
			if (fmax(fmin(err_a0, err_a1), fmax(fmin(err_b0, err_b1), fmin(err_c0, err_c1))) > 1e-10) {
				printf("ERROR on face %d\n", j);
				printf("vertices id: %d, %d, %d, ...\n", fvt[id_f1], fvt[id_f2], fvt[id_f3]);
				printf("Normal Vector (inter): (%f, %f, %f)\n", face_normals[ii], face_normals[jj], face_normals[kk]);
				printf("Normal Vector (pt123): (%f, %f, %f)\n", A / sqrt(A*A+B*B+C*C),
												  B / sqrt(A*A+B*B+C*C),
												  C / sqrt(A*A+B*B+C*C));
				printf("Normal Vector diff is too large!\n");
				exit(1);
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
			if (fabs(testp) > 1e-10) {
				printf("Projection point of sphere center: (%f, %f, %f)\n", p_x, p_y, p_z);
				printf("check if projection pt on plane: %g\n", testp);
				printf("Projection pt is not on the plane!\n");
				exit(1);
			}

			// double test1 = A*vt[fvt[tmp+4]*3] + B*vt[fvt[tmp+4]*3+1] + C*vt[fvt[tmp+4]*3+2] + D;
			double test1 = A*vt[fvt[tmp+1]*3] + B*vt[fvt[tmp+1]*3+1] + C*vt[fvt[tmp+1]*3+2] + D;
			if (fabs(test1) > 1e-10) {
				printf("check if First point on plane : %g\n", test1);
				printf("First point is not on the plane!\n");
				exit(1);
			}


			// 测试：计算三角形面积之和(face area) ================================================================
			// params: total_s_all
			double tmp_s_all = 0.0;
			for (unsigned int l=1; l<fo[j]-1; l++) {
				int pt_1 = 1;
				int pt_2 = l+1;
				int pt_3 = l+2;
				// printf("calculating triangle: %d, %d, %d\n", fvt[tmp+pt_1], fvt[tmp+pt_2], fvt[tmp+pt_3]);
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
			}


			// 循环每个顶点 =========================================================================================
			// params: a, b, c, p, S, total_s
			// printf("iterating over vertices:\n");
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
				// printf("area, %f\n", tmp_s);
				total_s += tmp_s;

				Y = tmp_s*2.0/tmp_c;

				double B1, B3;
				B1 = Y/H;
				B3 = tmp_a/Y;
				V1 += 0.25/pi*acos(1/B3)-1.0/8.0/pi*asin( ((1.-pow(B1,2))*pow(B3,2)-2.) / ((1.+pow(B1,2))*pow(B3,2)) ) - 1./16.;

				B3 = tmp_b/Y;
				V1 += 0.25/pi*acos(1/B3)-1.0/8.0/pi*asin( ((1.-pow(B1,2))*pow(B3,2)-2.) / ((1.+pow(B1,2))*pow(B3,2)) ) - 1./16.;
			}

			if (neigh[j] > 0)	V1_all += V1;

			// 用方法2计算角系数
			double xj = 0, yj = 0, zj = 0;
			bool isNeighFound = false;
			// printf("%d, %d\n", pid, neigh[j]);
			// if(cl2.start()) do if (con.compute_cell(c2,cl2)) {
			if(cl2.start()) do {
				if (neigh[j]<0) { isNeighFound = true; break; }
				if (cl2.pid() == neigh[j]) {
					cl2.pos(xj, yj, zj);
					isNeighFound = true;
					break;
				}
			} while (cl2.inc());
			if (!isNeighFound) {
				printf("NEIGH NOT FOUND!\n");
				exit(1);
			}
			// printf("pos: (%g, %g, %g)\n", xj, yj, zj);
			double L = 0;
			double Ai = 0;
			if (neigh[j] > 0) {
				double R = 0.03;
				// printf("neigh pos: (%g, %g, %g)\n", xj, yj, zj);
				L = sqrt(pow(x - xj, 2) + pow(y - yj, 2) + pow(z - zj, 2));
				if (L < 0.059) {
					printf("ERR - short distance found: %g\n", L);
				}
				// double d_max = 0.5 * sqrt(1 + (4.0 * face_areas[j]) / (pi*L*L));
				// double cos_cita = 0.5 / (d_max);
				double cos_cita = sqrt(pi*L*L / (4*face_areas[j] + pi*L*L));
				double sin_cita = sqrt(1 - cos_cita*cos_cita);
				double L_prime = L - R*cos_cita - R*cos_cita ;
				if (L_prime < 0) printf("L_prime < 0, %f\n", L_prime);
				double R_prime = R * sin_cita;
				double Z = 2 + L_prime*L_prime / R_prime / R_prime;
				double F_ij_prime = 0.5 * (Z - sqrt(Z*Z - 4));
				V2 = (1+cos_cita) / 2.0 * F_ij_prime;
				// V2 = (pi*R_prime*R_prime) / face_areas[j] * F_ij_prime;
				// printf("Ai_p/Ai, method1: %g, method2: %g\n", (1+cos_cita) / 2.0, (pi*R_prime*R_prime) / face_areas[j]);

				// double RR = R_prime / L_prime;
				// F_ij_prime = 1.0+(1.-sqrt(4*RR*RR+1))/2/RR/RR;
				// V2 = (1+cos_cita) / 2.0 * F_ij_prime;

				Ai = 2*pi*R*R*(1-cos_cita);
			}

			if (neigh[j] > 0)	V2_all += V2;

			fprintf(f1, "%6d, %6d, %.18f, %.18f, %.18f, %.18f\n", pid, neigh[j], L, V1, Ai, V2);
			tmp += k;
			// printf("face area (internal): %g\n", face_areas[j]);
			// printf("face area (me)      : %g\n", total_s);
			// if (fabs(total_s - face_areas[j]) > 1e-5) {
			// 	printf("face area error: %g\n", total_s - face_areas[j]);
			// 	printf("face area (internal): %g\n", face_areas[j]);
			// 	printf("face area (me)      : %g\n", total_s);
			// 	printf("-------------------------------------\n");
			// }

			// if (j == 2) break;
		}

		// printf("TOTAL VIEW FACTOR FOR FACE: %g\n", V1_all);
		// fprintf(f4, "%.19f\n", V1_all);
		fprintf(f4, "%g, %g\n", V1_all, V2_all);

		// if (index == 10) exit(1);
	} while (cl.inc());

	fclose(f1);
	fclose(f2);
	fclose(f3);
	fclose(f4);
}
