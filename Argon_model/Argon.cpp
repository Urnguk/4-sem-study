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

double find_arg(double argue, double length) {
	double arg_1 = fmod(argue, length);
	double arg_2;
	if (arg_1 >= 0.) {
		arg_2 = arg_1 - length;
	}
	else {
		arg_2 = arg_1 + length;
	}
	if (abs(arg_1) >= abs(arg_2)) {
		return arg_2;
	}
	return arg_1;
}

double count_potential(Molecule *system, int N, double length, double temperature, int time) {
	double potential = 0;
	double d_x, d_y, d_z;
	double r;
	double F;

	double T0 = 300;

	for (int i = 0; i < N - 1; i++) {
		for (int j = i + 1; j < N; j++) {

			d_x = find_arg((system[j].x - system[i].x), length);
			d_y = find_arg((system[j].y - system[i].y), length);
			d_z = find_arg((system[j].z - system[i].z), length);

			r = pow((pow(d_x, 2) + pow(d_y, 2) + pow(d_z, 2)), 0.5);
			if (r == 0.) {
				r = 0.000001;
			}

			potential += (4 / pow(r, 12)) - (4 / pow(r, 6)); // сделали лишнее умножение, проверка Максвелла

			F = -1 * ((48 / pow(r, 14)) - (24 / pow(r, 8)));

			system[i].Ax += F * d_x;
			system[i].Ay += F * d_y;
			system[i].Az += F * d_z;
			system[j].Ax -= F * d_x;
			system[j].Ay -= F * d_y;
			system[j].Az -= F * d_z;

		}
	}
	return potential;
}
