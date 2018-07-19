/**
* Y forward Z up for exporting from blender
*
* usage : g++ objInter.cpp -I/usr/include/eigen3 -o objInter
* ./objInter < particleIntermediate
* OR 
* cat particles lowObjs/allObjs > particleIntermediate && ./objInter < particleIntermediate > particleInterpolInfo && cat particles particleInterpolInfo > interpolatedParticles
*/

#include <cstdio>
#include <iostream>
#include <set>

#include "Particle.h"

int main(int argc, char **argv) {
	int N, n, numMuscles, numNeighbours;
	std::cin >> numMuscles;

	std::vector<int> muscles;
	n = numMuscles;
	while(n--) {
		std::cin >> N; muscles.push_back(N);
	}

	std::vector<Particle> particles(N);
	
	n = N;
	while(n--) { std::vector<double> info(7); double idx;
		std::cin >> idx >> info[0] >> info[1] >> info[2] >> info[3] >> info[4] >> info[5] >> info[6];
		particles[int(idx)] = Particle(info);
	}

	std::cin >> numNeighbours;
	n = N;	// not doing anything with the group information
	while(n--) { int index, num;
		std::cin >> index;
		for(int i = 0; i < numNeighbours; ++i) {
			std::cin >> num;
		}
	}

	std::cin >> numMuscles;
	if (numMuscles != muscles.size()) {
		printf("number of muscles do not match between particles and objects");
		return 1;
	}

	int prevMuscleNum = 0, numVerticies = 0, numFaces = 0;
	for (int muscleCount = 0; muscleCount < numMuscles; ++muscleCount) {
		std::vector<Vector3d> nodes; double nx,ny,nz;
		std::vector<Vector3i> triangles; int t1,t2,t3;
		char newLine;
		std::cin >> n;
		while(n--) {
			std::cin >> newLine; 
			switch (newLine) {
				case 'v': std::cin >> nx >> ny >> nz; nodes.push_back(Vector3d(nx,ny,nz)); break;
				case 'f': std::cin >> t1 >> t2 >> t3; triangles.push_back(Vector3i(t1-1,t2-1,t3-1)); break;
				default : break;
			}
		}

		char onlyTranslation = 'r';// 't' when we do not want rotation
		int prox = 3; // interpolation over 'prox' closest particles
		std::vector<int> temp(prox);
		std::vector<double> tempd2(prox);
		std::vector<Vector3d> tempdVec(prox);
		std::vector<std::vector<int> > closePartices;
		std::vector<std::vector<double> > closeParticesd2;
		std::vector<std::vector<Vector3d> > closeParticesdVec;
		for (int i = 0; i < nodes.size(); ++i) {
			closePartices.push_back(temp); closeParticesd2.push_back(tempd2); closeParticesdVec.push_back(tempdVec);
			for (int p = 0; p < prox; ++p) {
				bool first = true;
				for (int n = prevMuscleNum; n < muscles[muscleCount]; ++n) {
					bool seen = false;
					for (int j = 0; j < p; ++j) {
						if (n == closePartices[i][j]) {seen = true; break;}
					}
					if (seen) continue;
					Vector3d dist;
					if (onlyTranslation == 't') {
						dist = nodes[i] - particles[n].pos;
					} else {
						dist = particles[n].mat3_rot.inverse()*(nodes[i] - particles[n].pos);
					}
					// Vector3d dist = nodes[i] - particles[n].pos;
					if (first) {
						first = false; 
						closePartices[i][p] = n; closeParticesd2[i][p] = dist.squaredNorm(); closeParticesdVec[i][p] = dist;
						continue;
					}
					if (dist.squaredNorm() < closeParticesd2[i][p]) {
						closePartices[i][p] = n; closeParticesd2[i][p] = dist.squaredNorm(); closeParticesdVec[i][p] = dist;
					}
				}
			}
		}
		std::cout << nodes.size() << " " << prox << " " << onlyTranslation << std::endl;
		for (int i = 0; i < nodes.size(); ++i) {
			VectorXd w = VectorXd::Unit(prox, 0);
			for (int p = 0; p < prox; ++p) w[p] = 1/closeParticesd2[i][p];
			// for (int p = 0; p < prox; ++p) w[p] = 100 - closeParticesd2[i][p];
			for (int p = 0; p < prox; ++p) {
				std::cout << closePartices[i][p] << " " << (w[p]/w.sum()) << " " << closeParticesdVec[i][p].x() << " " << closeParticesdVec[i][p].y() << " " << closeParticesdVec[i][p].z() << " ";
			}
			std::cout << std::endl;
		}
		std::cout << triangles.size() << std::endl;
		for (int i = 0; i < triangles.size(); ++i) {
			std::cout << triangles[i].x() << " " << triangles[i].y() << " " << triangles[i].z() << " " << 
			(nodes[triangles[i].y()] - nodes[triangles[i].x()]).norm() << " " << 
			(nodes[triangles[i].z()] - nodes[triangles[i].y()]).norm() << " " << 
			(nodes[triangles[i].x()] - nodes[triangles[i].z()]).norm() << std::endl;
		}

		prevMuscleNum = muscles[muscleCount];
		numVerticies += nodes.size();
		numFaces += triangles.size();
	}
//	std::cout << "number of vertices = 14801\nnumber of faces = 29545" << std::endl;
	return 0;
}