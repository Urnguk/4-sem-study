#include <iostream>
#include <fstream>
using namespace std;
#include <math.h>


class Molecule {
public:
	/* переменные */

	double x, y, z;
	double x_0, y_0, z_0;
	double Vx, Vy, Vz;
	double Ax = 0.;
	double Ay = 0.;
	double Az = 0.;
	double tmp_Ax = 0.;
	double tmp_Ay = 0.;
	double tmp_Az = 0.;
	double dt;

	/* функции */

	void input_Mol(double input_x, double input_y, double input_z, double input_v, double input_dt) {

		x = input_x;
		x_0 = x;
		y = input_y;
		y_0 = y;
		z = input_z;
		z_0 = z;

		if (rand() % 2) {
			Vx = input_v;
		}
		else {
			Vx = -1 * input_v;
		}
		if (rand() % 2) {
			Vy = input_v;
		}
		else {
			Vy = -1 * input_v;
		}
		if (rand() % 2) {
			Vz = input_v;
		}
		else {
			Vz = -1 * input_v;
		}

		dt = input_dt;

		return;
	}

	void move() {
		x += (Vx * dt) + (Ax * pow(dt, 2.) / 2);
		y += (Vy * dt) + (Ay * pow(dt, 2.) / 2);
		z += (Vz * dt) + (Az * pow(dt, 2.) / 2);

		tmp_Ax = Ax;
		tmp_Ay = Ay;
		tmp_Az = Az;

		Ax = 0;
		Ay = 0;
		Az = 0;

		return;
	}

	void respeed() {
		Vx += (Ax + tmp_Ax) * dt / 2;
		Vy += (Ay + tmp_Ay) * dt / 2;
		Vz += (Az + tmp_Az) * dt / 2;

		return;
	}
};
