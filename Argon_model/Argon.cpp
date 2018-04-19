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

int main() {

	/* вспомогательные константы, переменные и создание массивов */

	const int n = 10;
	int N = pow(n, 3);
	const int time_period = 100000;
	double const_Bolcman = 1.4 * pow(10, -23);
	double epsilon = 171.12 * pow(10, -29);

	int length = 8;
	double Temperature = 250;

	double axis_v = pow((const_Bolcman * Temperature / epsilon), 0.5);
	double time_step = pow(10, -6);
	double full_v;

	double d_2_r = 0;

	Molecule system[n*n*n];


	double Kinetic_energy;
	double Potential_energy;
	double Full_energy;

	ofstream fout_1;
	fout_1.open("Energy_data.txt");
	
	
	/* ввод данных молекул*/

	double length_step = length / n;
	int mol_number = 0;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				system[mol_number].input_Mol(length_step * i, length_step * j, length_step * k, axis_v, time_step);
				mol_number++;
			}
		}
	}
	
	/* запуск основного цикла */

	for (int current_time = 0; current_time < time_period; current_time++) {

		if (!(current_time % 1000)) {
			cout << current_time << endl;
		}

		for (mol_number = 0; mol_number < N; mol_number++) {
			system[mol_number].move();
		}

		Potential_energy = count_potential(system, N, length, Temperature, current_time);
		Kinetic_energy = 0;

		for (mol_number = 0; mol_number < N; mol_number++) {
			system[mol_number].respeed();
			full_v = pow(pow(system[mol_number].Vx, 2) + pow(system[mol_number].Vy, 2) + pow(system[mol_number].Vz, 2), 0.5);
			Kinetic_energy += pow(full_v, 2) / 2;
		}

		Temperature = (Kinetic_energy * 2 * epsilon) / (3 * const_Bolcman * N);

		Full_energy = Potential_energy + Kinetic_energy;
		fout_1 << current_time << " " << Kinetic_energy << " " << Potential_energy << " " << Full_energy << endl;
	}
	
	fout_1.close();
	
	return 0;
}
