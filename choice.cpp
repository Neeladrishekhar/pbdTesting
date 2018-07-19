#include <cstdio>
#include <iostream>
#include <set>
#include <algorithm>
#include "Particle.h"

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
	
	std::vector<std::vector<double> > manyInfo(N);
	std::vector<Particle> particles(N);
	std::vector<std::set<int> > growUps(N);
	std::vector<std::vector<int> > groups(N);
	
	for (int i = 0; i < N; ++i) {
		std::vector<double> info(7); double idx;
		std::cin >> idx >> info[0] >> info[1] >> info[2] >> info[3] >> info[4] >> info[5] >> info[6];
		particles[int(idx)] = Particle(info);
		manyInfo[int(idx)] = info;
	}
	// printf("%1.20f\n", particles[0].pos.z());
	std::cin >> numNeighbours;
	n = N;
	while(n--) { int index, num;
		std::cin >> index;
		for(int i = 0; i < numNeighbours; ++i) {
			std::cin >> num;
			// growUps[index].insert(num);
			// growUps[num].insert(index);
			groups[index].push_back(num);
		}
	}
	// for (int i = 0; i < N; ++i) {
	// 	groups[i] = std::vector<int> (growUps[i].size());
	// 	std::copy(growUps[i].begin(), growUps[i].end(), groups[i].begin());
	// }

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

	std::vector<int> choiceMuscles;
	choiceMuscles.push_back(0);
	choiceMuscles.push_back(6);
	// choiceMuscles.push_back(10);
	// choiceMuscles.push_back(13);

	std::sort(choiceMuscles.begin(), choiceMuscles.end());
	int numMusclesNew = choiceMuscles.size();
	std::vector<int> musclesNew(numMusclesNew);
	musclesNew[0] = (choiceMuscles[0] == 0) ? muscles[0] : (muscles[choiceMuscles[0]]-muscles[choiceMuscles[0]-1]);
	for (int i = 1; i < numMusclesNew; ++i) {
		musclesNew[i] = musclesNew[i-1] + (muscles[choiceMuscles[i]]-muscles[choiceMuscles[i]-1]);
	}

	std::vector<int> particleMap(N, -1);
	std::vector<int> revParticleMap;
	std::vector<Particle> particlesNew;
	for (int i = 0; i < numMusclesNew; ++i) {
		for (int part = (choiceMuscles[i] == 0) ? 0 : muscles[choiceMuscles[i]-1]; part < muscles[choiceMuscles[i]]; ++part) {
			particleMap[part] = particlesNew.size();
			particlesNew.push_back(particles[part]);
			revParticleMap.push_back(part);
			// musclesNew[i] = part+1;
		}
	}

	std::vector<std::vector<int> > groupsNew(particlesNew.size());
	for (int i = 0; i < particlesNew.size(); ++i) {
		for (int gr = 0; gr < groups[revParticleMap[i]].size(); ++gr){
		 	groupsNew[i].push_back(particleMap[groups[revParticleMap[i]][gr]]);
		}
	}

	// std::vector< std::vector<std::vector<int> > > Muscle_interPNew;
	for (int i = 0; i < numMusclesNew; ++i) {
		std::vector<std::vector<int> > interP;
		interP = Muscle_interP[choiceMuscles[i]];
		for (int node = 0; node < interP.size(); ++node) {
			for (int pr = 0; pr < interP[node].size(); ++pr) {
				interP[node][pr] = particleMap[interP[node][pr]];
			}
		}
		// Muscle_interPNew.push_back(interP);
		Muscle_interP[choiceMuscles[i]] = interP;
	}

	// std::cout <<  << std::endl;

	std::cout << numMusclesNew << std::endl;
	for (int i = 0; i < musclesNew.size(); ++i) {
		std::cout << musclesNew[i] << " ";
	} std::cout << std::endl;
	for (int i = 0; i < particlesNew.size(); ++i) {
		std::cout << i << " ";
		for (int inf = 0; inf < 7; ++inf) {
			std::cout << manyInfo[revParticleMap[i]][inf] << " ";
		} std::cout << std::endl;
	}
	std::cout << numNeighbours << std::endl;
	for (int i = 0; i < groupsNew.size(); ++i) {
		std::cout << i << " ";
		for (int j = 0; j < groupsNew[i].size(); ++j) {
			std::cout << groupsNew[i][j] << " ";
		} std::cout << std::endl;
	}
	for (int i = 0; i < numMusclesNew; ++i) {
		int muscleIndex = choiceMuscles[i]; char trans = (onlyTranslation) ? 't' : 'r';
		std::cout << Muscle_interP[muscleIndex].size() << " " << Muscle_interP[muscleIndex][0].size() << " " << trans << std::endl;
		for (int mp = 0; mp < Muscle_interP[muscleIndex].size(); ++mp) {
			for (int p = 0; p < Muscle_interP[muscleIndex][mp].size(); ++p) {
				std::cout << Muscle_interP[muscleIndex][mp][p] << " " << 
							 Muscle_interPw[muscleIndex][mp][p] << " " << 
							 Muscle_interVec[muscleIndex][mp][p].x() << " " << 
							 Muscle_interVec[muscleIndex][mp][p].y() << " " << 
							 Muscle_interVec[muscleIndex][mp][p].z() << " ";//std::endl;
				// std::cin >> temp[p] >> tempw[p] >> x >> y >> z;
				// tempdVec[p] = Vector3d(x,y,z);
			} std::cout << std::endl;
			// interP.push_back(temp); interPw.push_back(tempw); interVec.push_back(tempdVec);
		}
		// std::cin >> numTris; int t1,t2,t3;
		std::cout << Muscle_interTris[muscleIndex].size() << std::endl;
		for (int mt = 0; mt < Muscle_interTris[muscleIndex].size(); ++mt) {
			std::cout << Muscle_interTris[muscleIndex][mt].x() << " " <<
						 Muscle_interTris[muscleIndex][mt].y() << " " <<
						 Muscle_interTris[muscleIndex][mt].z() << " " <<
						 Muscle_interTrisLen[muscleIndex][mt].x() << " " <<
						 Muscle_interTrisLen[muscleIndex][mt].y() << " " <<
						 Muscle_interTrisLen[muscleIndex][mt].z() << std::endl;
			// std::cin >> t1 >> t2 >> t3 >> x >> y >> z;
			// interTris.push_back(Vector3i(t1,t2,t3));
			// interTrisLen.push_back(Vector3d(x,y,z));
		}
	}
	std::cout << particleMap.size() << " ";
	for (int i = 0; i < particleMap.size(); ++i) {
		std::cout << particleMap[i] << " ";
	} std::cout << std::endl;

	return 0;
}