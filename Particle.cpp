#include <vector>
#include <Eigen/Dense>

#include "Particle.h"

Particle::Particle(std::vector<double> info) {
	pos = Vector3d(info[0],info[1],info[2]);
	mat3_rot = Quaterniond(info[3],info[4],info[5],info[6]).toRotationMatrix();
	mass = 1.0;
	vel = Vector3d::Zero();
	omega = Matrix3d::Identity();
	size = Vector3d(0.25, 0.18, 0.125);
	A = mass * 0.2 * DiagonalMatrix<double,3>(size.x()*size.x(), size.y()*size.y(), size.z()*size.z()) * mat3_rot;
}

void Particle::set_rotation(Matrix3d rot) {
	mat3_rot = rot;
	A = mass * 0.2 * DiagonalMatrix<double,3>(size.x()*size.x(), size.y()*size.y(), size.z()*size.z()) * mat3_rot;
}