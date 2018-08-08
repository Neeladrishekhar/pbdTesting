#include <iostream>
#include <ctime>
// shapeTime: 1.68151
// collisionTime: 25.2939
// meshTime: 4.7047
// stressTime: 2.76962
// renderTime: 17.5294
// totalTime: 51.9792
/*
graph of comparision
gaussian smoothing
rotate movie
*/
using namespace std;

int main(int argc, char const *argv[])
{
	int count = 50;
	char t; int shapeTime,collisionTime,meshTime,stressTime,renderTime,idle_num;
	double shape=0.0,collision=0.0,mesh=0.0,stress=0.0,render=0.0;
	for (int i = 0; i < count; ++i) {
		cin >> t >> shapeTime >> collisionTime >> meshTime >> stressTime >> renderTime >> idle_num;
		shape += ((1000.0*shapeTime/idle_num)/CLOCKS_PER_SEC);
		collision += ((1000.0*collisionTime/idle_num)/CLOCKS_PER_SEC);
		mesh += ((1000.0*meshTime/idle_num)/CLOCKS_PER_SEC);
		stress += ((1000.0*stressTime/idle_num)/CLOCKS_PER_SEC);
		render += ((1000.0*renderTime/idle_num)/CLOCKS_PER_SEC);
	}
	cout << "shapeTime: " << (1.0*shape/count) << "\n"
		<< "collisionTime: " << (1.0*collision/count) << "\n"
		<< "meshTime: " << (1.0*mesh/count) << "\n"
		<< "stressTime: " << (1.0*stress/count) << "\n"
		<< "renderTime: " << (1.0*render/count) << "\n"
		<< "totalTime: " << (1.0*(shape+collision+mesh+stress+render)/count)
		<< endl;
	return 0;
}