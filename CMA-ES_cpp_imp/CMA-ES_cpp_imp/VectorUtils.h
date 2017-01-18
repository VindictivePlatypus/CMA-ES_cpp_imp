#include "Eigen\Dense"


using namespace Eigen;
using namespace std;


VectorXd subVec(VectorXd v, int start, int size) {
	VectorXd sub(size);
	for (int i = 0; i < size; i++) {
		sub[i] = v[start + i];
	}
	return sub;
}

void quickSortwIndex(VectorXd &v, VectorXi &index, int left, int right) {
	int i = left, j = right;
	int tmp;
	double tmpd;
	double pivot = v[(left + right) / 2];

	// partition
	while (i <= j) {
		while (v[i] < pivot)
			i++;
		while (v[j] > pivot)
			j--;
		if (i <= j) {
			tmpd = v[i];
			v[i] = v[j];
			v[j] = tmpd;	

			tmp = index[i];
			index[i] = index[j];
			index[j] = tmp;

			i++;
			j--;
		}
	};

	// recursion 
	if (left < j)
		quickSortwIndex(v, index, left, j);
	if (i < right)
		quickSortwIndex(v, index, i, right);
}

VectorXi sortVector(VectorXd& v) {
	VectorXi indexes = VectorXi::Zero(v.size());
	for (int i = 0; i < v.size(); i++) indexes[i] = i;

	quickSortwIndex(v, indexes, 0, v.size() - 1);

	return indexes;
}