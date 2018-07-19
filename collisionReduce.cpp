#include <cstdio>
#include <iostream>
#include <set>
#include "Particle.h"
#include "omp.h"

std::vector< std::vector<std::vector<int> > > Muscle_interP;
std::vector< std::vector<std::vector<double> > > Muscle_interPw;
std::vector< std::vector<std::vector<Vector3d> > > Muscle_interVec;
std::vector< std::vector<Vector3i> > Muscle_interTris;
std::vector< std::vector<Vector3d> > Muscle_interTrisLen;
bool onlyTranslation;
std::vector<int> muscles;

int main(int argc, char **argv) {
	int N, n, numMuscles, numNeighbours;
	std::cin >> numMuscles;

	n = numMuscles;
	while(n--) {
		std::cin >> N; muscles.push_back(N);
	}

	std::vector<Particle> particles(N);
	std::vector<std::set<int> > growUps(N);
	std::vector<std::vector<int> > groups(N);
	
	for (int i = 0; i < N; ++i) {
		std::vector<double> info(7); double idx;
		std::cin >> idx >> info[0] >> info[1] >> info[2] >> info[3] >> info[4] >> info[5] >> info[6];
		particles[int(idx)] = Particle(info);
	}
	// printf("%1.20f\n", particles[0].pos.z());
	std::cin >> numNeighbours;
	n = N;
	while(n--) { int index, num;
		std::cin >> index;
		for(int i = 0; i < numNeighbours; ++i) {
			std::cin >> num;
			growUps[index].insert(num);
			growUps[num].insert(index);
		}
	}
	for (int i = 0; i < N; ++i) {
		groups[i] = std::vector<int> (growUps[i].size());
		std::copy(growUps[i].begin(), growUps[i].end(), groups[i].begin());
	}

	for (int muscleCount = 0; muscleCount < numMuscles; ++muscleCount) {

		std::vector<std::vector<int> > interP;
		std::vector<std::vector<double> > interPw;
		std::vector<std::vector<Vector3d> > interVec;
		std::vector<Vector3i> interTris;
		std::vector<Vector3d> interTrisLen;

		int numNodes, prox, numTris; double x,y,z; char translation;
		std::cin >> numNodes >> prox >> translation;
		onlyTranslation = (translation == 't') ? true : false;
		std::vector<int> temp(prox);
		std::vector<double> tempw(prox);
		std::vector<Vector3d> tempdVec(prox);
		for (int i = 0; i < numNodes; ++i) {
			for (int p = 0; p < prox; ++p) {
				std::cin >> temp[p] >> tempw[p] >> x >> y >> z;
				tempdVec[p] = Vector3d(x,y,z);
			}
			interP.push_back(temp); interPw.push_back(tempw); interVec.push_back(tempdVec);
		}
		std::cin >> numTris; int t1,t2,t3;
		for (int i = 0; i < numTris; ++i) {
			std::cin >> t1 >> t2 >> t3 >> x >> y >> z;
			interTris.push_back(Vector3i(t1,t2,t3));
			interTrisLen.push_back(Vector3d(x,y,z));
		}

		Muscle_interP.push_back(interP);
		Muscle_interPw.push_back(interPw);
		Muscle_interVec.push_back(interVec);
		Muscle_interTris.push_back(interTris);
		Muscle_interTrisLen.push_back(interTrisLen);
	}

	double factor = 5;
	std::vector<std::vector<int> > probableMuscleCollIdx;
	std::vector<std::vector<int> > probableTriangleCollIdx;
	// #pragma omp parallel for
	for (int s = 0; s < particles.size(); ++s) {
		int muscleIndex;
		for (int i = 0; i < muscles.size(); ++i) {
			if (s < muscles[i]) {muscleIndex = i; break;}
		}
		std::vector<int> muscleIdxs;
		std::vector<int> triangleIdxs;
		
		for (int muscleCount = 0; muscleCount < Muscle_interTris.size(); ++muscleCount) {
			if (muscleCount == muscleIndex) {continue;} // do not bother with its own muscle
			
			std::vector<std::vector<int> > interP = Muscle_interP[muscleCount];
			std::vector<std::vector<double> > interPw = Muscle_interPw[muscleCount];
			std::vector<std::vector<Vector3d> > interVec = Muscle_interVec[muscleCount];
			std::vector<Vector3i> interTris = Muscle_interTris[muscleCount];
			// std::vector<Vector3d> interTrisLen = Muscle_interTrisLen[muscleCount];

			for (int i = 0; i < interTris.size(); ++i) {
				std::vector<Vector3d> tria(3);
				for (int t = 0; t < 3; ++t) { tria[t] = Vector3d::Zero();
					for (int p = 0; p < interP[interTris[i][t]].size(); ++p) {
						if (onlyTranslation) {
							tria[t] += interPw[interTris[i][t]][p] * (particles[interP[interTris[i][t]][p]].pos + interVec[interTris[i][t]][p]);
						} else {
							tria[t] += interPw[interTris[i][t]][p] * (particles[interP[interTris[i][t]][p]].pos + (particles[interP[interTris[i][t]][p]].mat3_rot*interVec[interTris[i][t]][p]));
						}
					}
				}
				if ((tria[0]-particles[s].pos).norm() < factor*particles[s].size.x() ||
					(tria[1]-particles[s].pos).norm() < factor*particles[s].size.x() ||
					(tria[2]-particles[s].pos).norm() < factor*particles[s].size.x() ) {
					muscleIdxs.push_back(muscleCount); triangleIdxs.push_back(i);
				}
			}
		}
		// std::cout << " " << s << " " << muscleIdxs.size() << std::endl;
		probableMuscleCollIdx.push_back(muscleIdxs); probableTriangleCollIdx.push_back(triangleIdxs);
	}
	for (int pmci = 0; pmci < probableMuscleCollIdx.size(); ++pmci) {
		std::cout << pmci << " " << probableMuscleCollIdx[pmci].size() << std::endl;
		for (int indx = 0; indx < probableMuscleCollIdx[pmci].size(); ++indx) {
			std::cout << probableMuscleCollIdx[pmci][indx] << " " << probableTriangleCollIdx[pmci][indx] << " ";
		} if (probableMuscleCollIdx[pmci].size() != 0) std::cout << std::endl;
	}

	return 0;
}