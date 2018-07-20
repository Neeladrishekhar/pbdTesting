/*
cat tricep.1.node tricep.1.face > tricepTet
g++ -Wall -I/usr/include/eigen3 convObj.cpp -o convObj
./convObj < tricepTet > allObjs
*/
#include <cstdio>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;

int main(int argc, char const *argv[])
{
	std::vector<Vector3d> nodes;
	std::vector<Vector3i> faces;
	int numNodes, numFaces, t1, t2, t3, extra;
	double x, y, z;
	std::cin >> numNodes >> extra >> extra >> extra;
	while(numNodes--) {
		std::cin >> extra >> x >> y >> z;
		nodes.push_back(Vector3d(x,y,z));
	}
	std::cin >> numFaces >> extra;
	while(numFaces--) {
		std::cin >> extra >> t1 >> t2 >> t3 >> extra;
		faces.push_back(Vector3i(t1,t2,t3));
	}
	std::cout << "16\n" << (nodes.size()+faces.size()) << std::endl;
	for (int i = 0; i < nodes.size(); ++i) {
		std::cout << "v " << nodes[i].x() << " " << nodes[i].y() << " " << nodes[i].z() << std::endl;
	}
	for (int i = 0; i < faces.size(); ++i) {
		std::cout << "f " << faces[i].x()+1 << " " << faces[i].z()+1 << " " << faces[i].y()+1 << std::endl;
	}
	for (int i = 0; i < 15; ++i) {
		std::cout << "0 ";
	} std::cout << std::endl;
	return 0;
}