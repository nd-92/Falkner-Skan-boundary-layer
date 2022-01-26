// Requisite includes
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>

// Requisite namespace
using std::cout;
using std::endl;
using std::ofstream;
using std::string;

//Input solver parameters
const unsigned Np = 10000001;
double beta = 0, etaMax = 10;
double initialGuessMin = 0, initialGuessMax = 1;
double solutionTolerance = 1e-15;

// Preallocate grids and create some constants
double solutionWallShear;
double eta[Np], f[Np], f1[Np], f2[Np];
double A[11], difference[11], finalValue[11];
bool exportData = 0;
double m = -beta / (beta - 2);
double mCoeff = 0.5 * (m + 1);
unsigned NpMinOne = Np - 1;

// Declaration of functions
double createSolutionGrid() {
	for (unsigned i = 0; i < Np; i++) {
		eta[i] = etaMax * double(i) / double(Np - 1);
		f[i] = 0;
		f1[i] = 0;
		f2[i] = 0;
	}
	return 0;
}
double createGuessGrid(double initialGuessMin, double initialGuessMax) {
	for (unsigned i = 0; i < 11; i++) {
		A[i] = initialGuessMin + ((initialGuessMax - initialGuessMin) * double(i) / 10);
		difference[i] = 0;
		finalValue[i] = 0;
	}
	return 0;
}
double rungeKuttaFalknerSkan(double wallShear) {
	// Assign f''(0) and define constants
	f2[0] = wallShear;
	// Runge-Kutta solution
	for (unsigned k = 0; k < NpMinOne; k++) {
		unsigned kPlusOne = k + 1;
		double dy = eta[kPlusOne] - eta[k];
		double halfdy = 0.5 * dy;
		double K_1 = f1[k];
		double L_1 = f2[k];
		double M_1 = -m;
		M_1 = M_1 - (mCoeff * f[k] * f2[k]);
		M_1 = M_1 + (m * f1[k] * f1[k]);
		double K_2 = f1[k] + (K_1 * halfdy);
		double L_2 = f2[k] + (L_1 * halfdy);
		double M_2 = -m;
		M_2 = M_2 - (mCoeff * (f[k] + (halfdy * K_1)) * (f2[k] + (halfdy * M_1)));
		M_2 = M_2 + (m * (f1[k] + (halfdy * L_1)) * (f1[k] + (halfdy * L_1)));
		double K_3 = f1[k] + (K_2 * halfdy);
		double L_3 = f2[k] + (L_2 * halfdy);
		double M_3 = -m;
		M_3 = M_3 - (mCoeff * (f[k] + (halfdy * K_2)) * (f2[k] + (halfdy * M_2)));
		M_3 = M_3 + (m * (f1[k] + (halfdy * L_2)) * (f1[k] + (halfdy * L_2)));
		double K_4 = f1[k] + (K_3 * dy);
		double L_4 = f2[k] + (L_3 * dy);
		double M_4 = -m;
		M_4 = M_4 - (mCoeff * (f[k] + (dy * K_3)) * (f2[k] + (dy * M_3)));
		M_4 = M_4 + (m * (f1[k] + (dy * L_3)) * (f1[k] + (dy * L_3)));
		f[kPlusOne] = f[k] + ((dy / 6) * (K_1 + K_4 + (2 * (K_2 + K_3))));
		f1[kPlusOne] = f1[k] + ((dy / 6) * (L_1 + L_4 + (2 * (L_2 + L_3))));
		f2[kPlusOne] = f2[k] + ((dy / 6) * (M_1 + M_4 + (2 * (M_2 + M_3))));
	}
	return 0;
}
double loopFalknerSkan(double &initialGuessMin, double &initialGuessMax) {
	// Create guess grids
	createGuessGrid(initialGuessMin, initialGuessMax);
	// First loop: solve over a range of initial guesses
	for (unsigned q = 0; q < 11; q++) {
		rungeKuttaFalknerSkan(double(A[q]));
		finalValue[q] = f1[NpMinOne];
		difference[q] = abs(1 - finalValue[q]);
	}
	// Get value of wall shear for u_inf = 1
	unsigned minIndex = 0;
	unsigned iMinOne, iPlusOne;
	double minVal = difference[minIndex];
	for (unsigned i = 1; i < 11; i++) {
		if (difference[i] < minVal) {
			iMinOne = i - 1;
			iPlusOne = i + 1;
			minVal = difference[i];
			minIndex = i;
			initialGuessMin = A[iMinOne];
			initialGuessMax = A[iPlusOne];
		}
	}
	// Return new guess
	return initialGuessMin, initialGuessMax;
}

int main()
{
	// Set initial clock time
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	// Set precision of cout
	cout.precision(15);
	// Initialise solution grids
	createSolutionGrid();
	// Converge guess for wall shear
	while (initialGuessMax - initialGuessMin > solutionTolerance) {
		loopFalknerSkan(initialGuessMin, initialGuessMax);
		cout << (initialGuessMax + initialGuessMin) / 2 << endl;
	}
	solutionWallShear = (initialGuessMax + initialGuessMin) / 2;
	// Solve for u_inf = 1
	rungeKuttaFalknerSkan(solutionWallShear);
	// Output solution wall shear to console
	cout << "Wall shear = " << f2[0] << endl;
	// Print elapsed time to console
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	cout << "Elapsed time is " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() * 0.000001 << " seconds." << endl;
	// Write solution to file
	ofstream writeFile("Solution.dat");
	if (writeFile.is_open()) {
		writeFile.precision(15);
		writeFile << "eta \t f \t f1 \t f2 \n";
		for (int i = 0; i < Np; i++) {
			writeFile << eta[i] << "\t" << f[i] << "\t" << f1[i] << "\t" << f2[i] << endl;
		}
		writeFile.close();
	}
	else {
		cout << "Unable to open file";
	}
}