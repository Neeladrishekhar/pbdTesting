#include <vector>
#include <Eigen/Dense>
using namespace Eigen;

class Particle
{
public:
	Vector3d pos;
	Matrix3d mat3_rot;
	Matrix3d A; //moment matrix
	double mass;
	Vector3d vel;
	Matrix3d omega;
	Vector3d size;

	Particle(){}
	Particle(std::vector<double> info);
	~Particle(){}
	void set_rotation(Matrix3d rot);
};