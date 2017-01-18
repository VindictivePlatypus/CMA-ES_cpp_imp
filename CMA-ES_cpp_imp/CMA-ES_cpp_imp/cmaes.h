/*includes go here*/
#include "Eigen\Dense"
#include <math.h>
#include <random>
#include <chrono>
#include "VectorUtils.h"


/*declarations go here*/

using namespace Eigen;
using namespace std;

class CMAES {
public:
	CMAES(int s, double (*obj) (const double* x), const double* start = NULL, double sigma = 0.3, double stopF = 1e-10, unsigned int stopE = 1e100){
		// User defined parameters
		N = s;
		objective = obj;
		if (start) {
			curX = VectorXd(N);
			for (int i = 0; i < s; i++) {
				curX[i] = start[i];
			}
		}
		else curX = (VectorXd::Random(N) + VectorXd::Ones(N)) / 2;
		sig = sigma;
		stopFit = stopF;
		stopEval = stopE;
		
		// Strategy parameter setting: Selection
		lambda = 4 + floor(3 * log(N));
		double muTemp = (double)lambda / 2.;
		mu = floor(muTemp);
		weights = VectorXd(mu);
		for (int i = 1; i <= muTemp; i++) {
			weights[i - 1] = log(muTemp + .5) - log(i);
		}
		weights /= weights.sum();
		ArrayXd tmpV(weights);
		mueff = pow(weights.sum(),2) / tmpV.square().sum();

		// Strategy parameter setting: Adaptation
		cc = (4 + (mueff / N)) / (N + 4 + 2 * (mueff / N));
		cs = (mueff + 2) / (N + mueff + 5);
		c1 = 2 / (pow(N + 1.3, 2) + mueff);
		cmu = min(1 - c1, 2 * (mueff - 2 + (1 / mueff)) / (pow(N + 2, 2) + mueff));
		damps = 1 + 2 * max(0., sqrt((mueff - 1) / (N + 1)) - 1) + cs;

		// Initialize dynamic(internal) strategy parameters and constants
		pc = VectorXd::Zero(N);
		ps = VectorXd::Zero(N);
		B = MatrixXd::Identity(N,N);
		D = VectorXd::Ones(N);
		ArrayXd tmpV2(D);
		tmpV2 = tmpV2.square();
		VectorXd tmpV3(tmpV2);
		C = B * tmpV3.asDiagonal() * B.transpose();
		invSqrtC = B * D.asDiagonal().inverse() * B.transpose();
		eigenEval = 0;
		chiN = sqrt(N)*(1 - (1 / (4 * N)) + (1 / (21 * N * N)));
		
		// Other initializations
		arx = MatrixXd::Zero(N, lambda);
		arFitness = VectorXd::Zero(lambda);
		arIndexes = VectorXi::Zero(lambda);

		/**/
	}

	void step() {
		// inside "while counteval < stopeval" in matlab file
		//Generation
		for (int k = 0; k < lambda; k++) {
			unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
			default_random_engine generator(seed);
			normal_distribution<double> distr(0.0, 1.0);
			VectorXd temp(N);
			for (int i = 0; i < N; i++) {
				temp[i] = distr(generator);
			}
			arx.col(k) = curX + sig * B * (D.cwiseProduct(temp));
			arFitness[k] = objective(arx.col(k).data());
			++countEval;
		}
		//Sort
		for (int i = 0; lambda; i++) arIndexes[i] = i;
		quickSortwIndex(arFitness, arIndexes,0,lambda-1);
		xOld = VectorXd(curX);
		MatrixXd tempMat(N, mu);
		for (int i = 0; mu; i++) tempMat.col(i) = arx.col(arIndexes[i]);
		curX = tempMat*weights;

		//Cumulation
		ps = (1 - cs)*ps + sqrt(cs*(2 - cs)*mueff) * invSqrtC * (curX - xOld) / sig;
		hsig = ps.norm() / (sqrt(1-pow(1-cs,2*countEval/lambda)))/chiN < 1.4 + 2 / (N+1);
		pc = (1 - cc)*pc + hsig*sqrt(cc*(2-cc)*mueff)*(curX - xOld)/sig;
	}

	void cma() {

	}

	void setStartingPoint(const double* x) {
		if (countEval == 0) {
			//curX = VectorXd(x);
		}
	}

	void setSigma(double s) {
		sig = s;
	}

	void setStopFit(double sf) {
		stopFit = sf;
	}

	void setStopEval(unsigned int se) {
		stopEval = se;
	}

	void resetCMA() {
		countEval = 0;
	}

	


private:
	int N;
	double (*objective) (const double* x);
	double sig, stopFit, mueff, cc, cs, c1, cmu, damps, chiN;
	unsigned int stopEval, countEval, eigenEval;
	int lambda, mu;
	bool hsig;

	VectorXd weights;
	VectorXd curX, xOld;
	VectorXd pc, ps, D; // For D, think to use a temp matrix when calculating eigenvalues
	MatrixXd B, C, invSqrtC;

	MatrixXd arx;
	VectorXd arFitness;
	VectorXi arIndexes;
};