#include <iostream>
#include <random>
#include <chrono>
#include "cmaes.h"

using namespace std;


int main() {
	VectorXd v(5);
	for (int i = 0; i < 5; i++)
	{
		v[i] = 15-i;
	}
	/*
	int* t = v.data();
	cout << t[1];
	//v /= v.sum();
	//double k = v.sum();
	//ArrayXd v1(v);
	//cout << v1.square();
	VectorXd pc, ps, D; // For D, think to use a temp matrix when calculating eigenvalues
	MatrixXd B, C, invSqrtC;
	pc = VectorXd::Zero(N);
	ps = VectorXd::Zero(N);
	B = MatrixXd::Identity(N,N);
	VectorXd D = VectorXd::Ones(N);
	ArrayXd tmpV2(D);
	tmpV2 = tmpV2.square();
	VectorXd tmpV3(tmpV2);
	C = B * tmpV3.asDiagonal() * B.transpose();
	invSqrtC = B * D.asDiagonal().inverse() * B.transpose();
	//cout << v.pow(2);
	cout << pc << endl;
	cin.ignore();
	cout << ps << endl;
	cin.ignore();

	cout << D << endl;
	cin.ignore();
	cout << B << endl;
	cin.ignore();
	cout << C << endl;
	cin.ignore();
	cout << invSqrtC;
	cin.ignore();

	cout << test
	MatrixXd test(5, 3);
	test.col(0) = v;
	cout << test.col(0);
	//cout << v;;
	*/
	int N = 5;
	VectorXd D = VectorXd::Ones(N)*2;
	cout << v << endl;
	cout << sortVector(v);
	cout << v;
	cin.ignore();
	return 0;
}