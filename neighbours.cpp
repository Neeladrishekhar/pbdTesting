/*
g++ neighbours.cpp -I/usr/include/eigen3 -o neighbours
./neighbours < tempfile.txt > tempNeigh.txt && cat tempfile.txt tempNeigh.txt > particles

*
*
./neighbours < tempfile.txt > tempNeigh.txt && cat tempfile.txt tempNeigh.txt > particles && cat particles lowObjs/allObjs > particleIntermediate && ./objInter < particleIntermediate > particleInterpolInfo && cat particles particleInterpolInfo > interpolatedParticles && ./stressTest < interpolatedParticles

*/

#include <cstdio>
#include <iostream>
#include <utility>
#include <set>

#include "Particle.h"

int main(int argc, char **argv) {
	int N, n, numMuscles;
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

	int numNeighbours = 6;
	std::cout << numNeighbours << std::endl;
	int prevMuscle = 0;
	for (int i = 0; i < numMuscles; ++i) {
		for (int j = prevMuscle; j < muscles[i]; ++j) {
			std::vector<std::pair<double,int> >distances;
			// std::vector<double> distances(muscles[i]-prevMuscle);
			for (int k = prevMuscle; k < muscles[i]; ++k) {
				distances.push_back(std::make_pair((particles[k].pos - particles[j].pos).norm(),k));
				// distances[k-prevMuscle] = (particles[k].pos - particles[j].pos).norm();
			}
			sort(distances.begin(),distances.end());
			std::cout << j;
			for (int k = 1; k < numNeighbours+1; ++k) {
				std::cout << " " << distances[k].second;
			}
			std::cout << std::endl;
		}
		prevMuscle = muscles[i];
	}

	return 0;
}