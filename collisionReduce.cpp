#include <cstdio>
#include <iostream>
#include <set>
#include "Particle.h"
// #include "omp.h"

std::vector< std::vector<std::vector<int> > > Muscle_interP;
std::vector< std::vector<std::vector<double> > > Muscle_interPw;
std::vector< std::vector<std::vector<Vector3d> > > Muscle_interVec;
std::vector< std::vector<Vector3i> > Muscle_interTris;
std::vector< std::vector<Vector3d> > Muscle_interTrisLen;
bool onlyTranslation;
std::vector<int> muscles;

std::vector<std::vector<Vector3d> > boneVerticies;
std::vector<std::vector<Vector3i> > boneTriangles;

int main(int argc, char **argv) {
	bool remapReq = false, useBones = false;
	if (argc == 2) {
		if (argv[1][0] == 'l') remapReq = true;
		else if (argv[1][0] == 'b') useBones = true;
	} else if (argc == 3) {
		useBones = true; remapReq = true;
	}

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
	
	if (remapReq) {
		int wholeN; std::cin >> wholeN;
		std::vector<int> particleMap(wholeN, -1);
		for (int i = 0; i < wholeN; ++i) std::cin >> particleMap[i];
	}

	if (useBones) {
		int numBones,numVerticies,numTriangles,t1,t2,t3; 
		double x,y,z; char vORf;
		std::cin >> numBones;
		while(numBones--) {
			std::vector<Vector3d> verticies;
			std::cin >> numVerticies;
			while(numVerticies--) {
				std::cin >> vORf >> x >> y >> z; verticies.push_back(Vector3d(x,y,z));
			}
			boneVerticies.push_back(verticies);

			std::vector<Vector3i> triangles;
			std::cin >> numTriangles;
			while(numTriangles--) {
				std::cin >> vORf >> t1 >> t2 >> t3; triangles.push_back(Vector3i(t1-1,t2-1,t3-1));
			}
			boneTriangles.push_back(triangles);
		}
	}

	double factor = 10, defFactor = 5;
	std::vector<std::vector<int> > probableMuscleCollIdx;
	std::vector<std::vector<int> > probableTriangleCollIdx;
	std::vector<std::vector<int> > probableBoneCollIdx;
	std::vector<std::vector<int> > probableTriBoneCollIdx;
	// #pragma omp parallel for
	for (int s = 0; s < particles.size(); ++s) {
		int muscleIndex;
		for (int i = 0; i < muscles.size(); ++i) {
			if (s < muscles[i]) {muscleIndex = i; break;}
		}
		if (muscleIndex == 8) factor = 15; else factor = defFactor;
		// if (muscleIndex == muscles.size()-1) factor = 15; else factor = defFactor;
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

		if (useBones){
			std::vector<int> boneIdxs;
			std::vector<int> triBoneIdxs;
			for (int boneIdx = 0; boneIdx < boneTriangles.size(); ++boneIdx) {
				for (int triaIdx = 0; triaIdx < boneTriangles[boneIdx].size(); ++triaIdx) {
					std::vector<Vector3d> tria(3);
					tria[0] = boneVerticies[boneIdx][boneTriangles[boneIdx][triaIdx].x()];
					tria[1] = boneVerticies[boneIdx][boneTriangles[boneIdx][triaIdx].y()];
					tria[2] = boneVerticies[boneIdx][boneTriangles[boneIdx][triaIdx].z()];
					if ((tria[0]-particles[s].pos).norm() < factor*particles[s].size.x() ||
						(tria[1]-particles[s].pos).norm() < factor*particles[s].size.x() ||
						(tria[2]-particles[s].pos).norm() < factor*particles[s].size.x() ) {
						boneIdxs.push_back(boneIdx); triBoneIdxs.push_back(triaIdx);
					}
				}
			}
			probableBoneCollIdx.push_back(boneIdxs); probableTriBoneCollIdx.push_back(triBoneIdxs);
		}
	}
	for (int pmci = 0; pmci < probableMuscleCollIdx.size(); ++pmci) {
		std::cout << pmci << std::endl;
		std::cout << probableMuscleCollIdx[pmci].size() << " ";
		for (int indx = 0; indx < probableMuscleCollIdx[pmci].size(); ++indx) {
			std::cout << probableMuscleCollIdx[pmci][indx] << " " << probableTriangleCollIdx[pmci][indx] << " ";
		// } if (probableMuscleCollIdx[pmci].size() != 0) std::cout << std::endl;
		} std::cout << std::endl;
		if (useBones){
			std::cout << probableBoneCollIdx[pmci].size() << " ";
			for (int indx = 0; indx < probableBoneCollIdx[pmci].size(); ++indx) {
				std::cout << probableBoneCollIdx[pmci][indx] << " " << probableTriBoneCollIdx[pmci][indx] << " ";
			} std::cout << std::endl;
		}
	}

	return 0;
}