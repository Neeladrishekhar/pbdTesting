/**
* ./tetgen -pRzf ../upperArmModels/arm/tricep.ply
* 
* usage : g++ stressTest.cpp -lglut -lGLU -lGL -lm -fPIC -I/usr/include/eigen3 -o stressTest
* ./stressTest < interpolatedParticles
*/

#include <cstdio>
#include <iostream>
#include <set>
#include "Particle.h"
#include "omp.h"

#include <GL/glut.h>

void tripa(float a, float b, float c) {
	Vector3d nor = ((Vector3d(0,0,c) - Vector3d(0,b,0)).cross(Vector3d(a,0,0) - Vector3d(0,b,0))).normalized();
	if (a*b*c < 0) nor *= -1;
	glBegin(GL_TRIANGLES);
	glNormal3f(nor.x(),nor.y(),nor.z());
	glVertex3f(a,0,0);glVertex3f(0,b,0);glVertex3f(0,0,c); 
	glEnd();
}

void render_particle(Vector3d abc, GLfloat the_col[]) {
	// glPushMatrix();	

	float scale = 2.0f, a = abc.x(), b = abc.y(), c = abc.z();
	a *= scale; b *= scale; c *= scale;
	glMaterialfv(GL_FRONT, GL_DIFFUSE, the_col);
	// glColor3f(0,1,0);
	tripa(a,b,c); tripa(a,b,-c); tripa(a,-b,c); tripa(a,-b,-c); 
	tripa(-a,b,c); tripa(-a,b,-c); tripa(-a,-b,c); tripa(-a,-b,-c);

	// glPopMatrix();
}

std::vector< std::vector<std::vector<int> > > Muscle_interP;
std::vector< std::vector<std::vector<double> > > Muscle_interPw;
std::vector< std::vector<std::vector<Vector3d> > > Muscle_interVec;
std::vector< std::vector<Vector3i> > Muscle_interTris;
std::vector< std::vector<Vector3d> > Muscle_interTrisLen;
bool onlyTranslation;

std::vector<bool> controlParticles;
std::vector<std::vector<int> > boneParticles(3);
int clavicle_arr[] = {105,106,107,108};
std::vector<int> clavicle(clavicle_arr, clavicle_arr + sizeof(clavicle_arr) / sizeof(int) );
int scapula_arr[] = {0,1,2,43,44,45,46,74,75,76,183,216,217,218,219,250};
std::vector<int> scapula(scapula_arr, scapula_arr + sizeof(scapula_arr) / sizeof(int) );
Vector3d humerusPivot(7.5779,-16.2654,131.5013);
Matrix3d humerusRot = Matrix3d::Identity();
int humerus_arr[] = {47,48,49,98,101,102,104,136,139,140,141,142,143,144,145,148,149,157,159,160,165,169,172,174,186,193,194,195,197,
				199,200,201,204,208,209,246,247,248,249,251,252,253,254,281,282,283,284,285,286,288,289,290,291,317,318,361,383,
				384,386,404,405};
std::vector<int> humerus(humerus_arr, humerus_arr + sizeof(humerus_arr) / sizeof(int) );
Vector3d ulnaPivot(6.4562,-20.8208,103.8195);
Matrix3d ulnaRot = Matrix3d::Identity();
int ulna_arr[] = {40,41,180,181,182,196,214,215,313,314,315,316,349,362,365,370,371,372,377,378,380,381,382,399,401,402,403,408};
std::vector<int> ulna(ulna_arr, ulna_arr + sizeof(ulna_arr) / sizeof(int) );
int radius_arr[] = {278,279,280,345,346,347,348,354,357,364,369,392,393,394};
std::vector<int> radius(radius_arr, radius_arr + sizeof(radius_arr) / sizeof(int) );
// int stuck_arr[] = {0,1,2,43,44,45,46,74,75,76,105,106,107,108,143,144,145,146,147,183,194,201,202,216,217,218,219,220,250,251,252,253,254,281,317,349,359,360,361,370,383,384,386,402,404,405};
// std::vector<int> stuck(stuck_arr, stuck_arr + sizeof(stuck_arr) / sizeof(int) );
// int moved_arr[] = {40,41,47,48,49,101,102,103,104,136,139,140,141,142,180,181,182,193,195,196,197,198,246,247,248,249,278,279,280,313,314,315,316,344,345,346,347,348,354,369,382,394,399};
// std::vector<int> moved(moved_arr, moved_arr + sizeof(moved_arr) / sizeof(int) );
std::vector<std::vector<Vector3d> > boneVerticies;
std::vector<std::vector<Vector3i> > boneTriangles;

std::vector<std::vector<bool> > collisionTriangles;
std::vector<std::vector<bool> > collBoneTriangles;
bool handleCollision = false;
bool useBones = false;
bool pauseManual = false;

class Shape
{
public:	
	std::vector<Particle> particles;
	std::vector<std::vector<int> > groups;

	Shape(){}
	~Shape(){}

	void render() {
		// std::cout << "shape render taken!" << std::endl;
		GLfloat red[] = {.8f, 0.f, 0.f, 1.f};
		GLfloat green[] = {0.f, .8f, 0.f, 1.f};
		GLfloat yellow[] = {.8f, .8f, 0.f, 1.f};
		GLfloat tBlue[] = {0.f, 0.f, .8f, 1.f};
		GLfloat tMusc[] = {.8f, .2f, .1f, 1.f};
		// GLfloat blue[] = {0.f, 0.f, .8f, 1.f};
		// int stu = 0,mov = 0; 
		Matrix4d mat4 = Matrix4d::Identity();
		for (int i = 0; i < particles.size(); ++i) {
			glTranslatef(particles[i].pos.x(), particles[i].pos.y(), particles[i].pos.z());
			mat4.block(0,0,3,3) = particles[i].mat3_rot;
			glMultMatrixd(mat4.data());

			if (controlParticles[i]) {
				render_particle(particles[i].size, red);
			} else {
				render_particle(particles[i].size, green);
			}
			glPopMatrix();
			glPushMatrix();
		}

		// glColor3f(1,0,0);
		for (int i = 0; i < groups.size(); ++i) {
			for (int j = 0; j < groups[i].size(); ++j) {
				glBegin(GL_LINES);
				glLineWidth(5.0);
				glMaterialfv(GL_FRONT, GL_DIFFUSE, red);
				glVertex3f(particles[i].pos.x(),particles[i].pos.y(),particles[i].pos.z());
				glVertex3f(particles[groups[i][j]].pos.x(),particles[groups[i][j]].pos.y(),particles[groups[i][j]].pos.z());
				glEnd();
			}
		}
		
		GLfloat tColor[] = {0.f, 0.f, .8f, .6f};
		for (int muscleCount = 0; muscleCount < Muscle_interTris.size(); ++muscleCount) { //if (muscleCount > 0) continue;
			// std::cout << "Hello World! " << omp_get_thread_num() << " of " << omp_get_num_threads() << std::endl;
			
			std::vector<std::vector<int> > interP = Muscle_interP[muscleCount];
			std::vector<std::vector<double> > interPw = Muscle_interPw[muscleCount];
			std::vector<std::vector<Vector3d> > interVec = Muscle_interVec[muscleCount];
			std::vector<Vector3i> interTris = Muscle_interTris[muscleCount];
			std::vector<Vector3d> interTrisLen = Muscle_interTrisLen[muscleCount];

			for (int i = 0; i < interTris.size(); ++i) {
				std::vector<Vector3d> tria(3);
				// std::cout << "render_func taken! " << i << " " << interTris.size() << std::endl;
				for (int t = 0; t < 3; ++t) { tria[t] = Vector3d::Zero();
					for (int p = 0; p < interP[interTris[i][t]].size(); ++p) {
						if (onlyTranslation) {
							tria[t] += interPw[interTris[i][t]][p] * (particles[interP[interTris[i][t]][p]].pos + interVec[interTris[i][t]][p]);
						} else {
							int tem = interP[interTris[i][t]][p];
							// if (i == 84) std::cout << "84 index tria : " << t << " " << p  << " " << interTris[i][t] << std::endl;
							tria[t] += interPw[interTris[i][t]][p] * (particles[interP[interTris[i][t]][p]].pos + (particles[interP[interTris[i][t]][p]].mat3_rot*interVec[interTris[i][t]][p]));
						}
					}
				}
				Vector3d nor = ((tria[1]-tria[0]).cross(tria[2]-tria[0])).normalized();

				if (handleCollision && collisionTriangles[muscleCount][i]) {
					if (!pauseManual) collisionTriangles[muscleCount][i] = false;
					glBegin(GL_TRIANGLES);
					glNormal3f(nor.x(),nor.y(),nor.z());
					glMaterialfv(GL_FRONT, GL_DIFFUSE, yellow);
					glVertex3f(tria[0].x(),tria[0].y(),tria[0].z());
					glVertex3f(tria[1].x(),tria[1].y(),tria[1].z());
					glVertex3f(tria[2].x(),tria[2].y(),tria[2].z());
					glEnd();
					continue;
				}
				
				Vector3d deltaLen((tria[1]-tria[0]).norm(),(tria[2]-tria[1]).norm(),(tria[0]-tria[2]).norm());
				deltaLen -= interTrisLen[i];
				Vector3d stress(deltaLen.x()/interTrisLen[i].x(), deltaLen.y()/interTrisLen[i].y(), deltaLen.z()/interTrisLen[i].z());
				stress *= 4;
				Vector3d colorScale(std::max(std::min(stress.x(),1.0),-1.0), std::max(std::min(stress.y(),1.0),-1.0), std::max(std::min(stress.z(),1.0),-1.0));
				std::vector<Vector3d> colors(3);
				// for (int t = 0; t < 3; ++t) { colors[t][0] = std::max(0.0, colorScale[t]); colors[t][2] = -1.0 * std::min(0.0, colorScale[t]); colors[t][1] = 1.0 - std::max(colorScale[t], -1*colorScale[t]); }
				double grey = 0.2;
				for (int t = 0; t < 3; ++t) { colors[t][0] = 1-((1 - std::max(0.0, colorScale[t]))*(1-grey)); colors[t][2] = 1-((1 - (-1.0 * std::min(0.0, colorScale[t])))*(1-grey)); colors[t][1] = grey; }

				// glBegin(GL_TRIANGLES);
				glBegin(GL_POLYGON);
				glNormal3f(nor.x(),nor.y(),nor.z()); Vector3d mid;
				for (int t = 0; t < 3; ++t) {
					mid = 0.5*(colors[(t+2)%3] + colors[t]);
					tColor[0] = mid[0]; tColor[2] = mid[2]; tColor[1] = mid[1];
					glMaterialfv(GL_FRONT, GL_DIFFUSE, tColor);
					glVertex3f(tria[t].x(),tria[t].y(),tria[t].z());
					tColor[0] = colors[t][0]; tColor[2] = colors[t][2]; tColor[1] = colors[t][1];
					glMaterialfv(GL_FRONT, GL_DIFFUSE, tColor);
					mid = 0.5*(tria[(t+1)%3] + tria[t]);
					glVertex3f(mid.x(),mid.y(),mid.z());
				}
				// glVertex3f(tria[0].x(),tria[0].y(),tria[0].z());glVertex3f(tria[1].x(),tria[1].y(),tria[1].z());
				// glVertex3f(tria[2].x(),tria[2].y(),tria[2].z());
				glEnd();
			}
		}

		if (useBones) {
			GLfloat boneColor[] = {.8f, .8f, .8f, 1.f};
			for (int boneIdx = 0; boneIdx < boneTriangles.size(); ++boneIdx) {
				for (int triaIdx = 0; triaIdx < boneTriangles[boneIdx].size(); ++triaIdx) {
					std::vector<Vector3d> tria(3);
					
					// scapula and clavicle
					tria[0] = boneVerticies[boneIdx][boneTriangles[boneIdx][triaIdx].x()];
					tria[1] = boneVerticies[boneIdx][boneTriangles[boneIdx][triaIdx].y()];
					tria[2] = boneVerticies[boneIdx][boneTriangles[boneIdx][triaIdx].z()];
					// humerus
					if (boneIdx == 2) {
						tria[0] = (humerusRot * (tria[0] - humerusPivot)) + humerusPivot;
						tria[1] = (humerusRot * (tria[1] - humerusPivot)) + humerusPivot;
						tria[2] = (humerusRot * (tria[2] - humerusPivot)) + humerusPivot;
					// ulna and radius
					} else if (boneIdx > 2) {
						tria[0] = (ulnaRot * (tria[0] - ulnaPivot)) + ulnaPivot; tria[0] = (humerusRot * (tria[0] - humerusPivot)) + humerusPivot;
						tria[1] = (ulnaRot * (tria[1] - ulnaPivot)) + ulnaPivot; tria[1] = (humerusRot * (tria[1] - humerusPivot)) + humerusPivot;
						tria[2] = (ulnaRot * (tria[2] - ulnaPivot)) + ulnaPivot; tria[2] = (humerusRot * (tria[2] - humerusPivot)) + humerusPivot;
					}

					Vector3d nor = ((tria[1]-tria[0]).cross(tria[2]-tria[0])).normalized();
					glBegin(GL_TRIANGLES);
					glNormal3f(nor.x(),nor.y(),nor.z());
					if (handleCollision && collBoneTriangles[boneIdx][triaIdx]) {
						if (!pauseManual) collBoneTriangles[boneIdx][triaIdx] = false;
						glMaterialfv(GL_FRONT, GL_DIFFUSE, yellow);
					} else {
						glMaterialfv(GL_FRONT, GL_DIFFUSE, boneColor);
					}
					glVertex3f(tria[0].x(),tria[0].y(),tria[0].z());
					glVertex3f(tria[1].x(),tria[1].y(),tria[1].z());
					glVertex3f(tria[2].x(),tria[2].y(),tria[2].z());
					glEnd();
				}
			}
		}
	}
};

Shape model, model_base;
Vector3d lookAt(3,-15,120);
// Vector3d eyePos(20,20,120);
Vector3d eyePos(7,24,120);

void render(){
	// std::cout << "Render taken!" << std::endl;
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(eyePos.x(),eyePos.y(),eyePos.z(), lookAt.x(),lookAt.y(),lookAt.z(), 0,0,1);// gluLookAt(eye, center, up);
	// gluLookAt(eyePos.x(),eyePos.y(),eyePos.z(), 3,-15,120, 0,0,1);
	// gluLookAt(eyePos.x(),eyePos.y(),eyePos.z(), 0,0,0, 0,0,1);

	glPushMatrix();

	model.render();

	glPopMatrix();

	glutSwapBuffers();
}

void reshape_func(int w, int h){
	glViewport(0, 0, (GLsizei) w, (GLsizei) h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0, float(1.0*w/h), 1.0, 1000.0);
	// glOrtho(0.0f, w, h, 0.0f, 0.0f, 1.0f);

	glutPostRedisplay();
}

void keyboard_func(unsigned char key, int x, int y){
	double cam_speed = 2.5;
	Matrix3d m;
	switch (key)
	{
	case 27:
		exit(0);
		break;
	case 'p':
		pauseManual = !pauseManual;
		break;
	case 'q':
		eyePos += Vector3d::UnitZ();lookAt += cam_speed * Vector3d::UnitZ(); glutPostRedisplay();
		break;
	case 'e':
		eyePos -= Vector3d::UnitZ();lookAt += -cam_speed * Vector3d::UnitZ(); glutPostRedisplay();
		break;
	case 'w':
		eyePos += cam_speed*(lookAt - eyePos).normalized(); glutPostRedisplay();
		break;
	case 's':
		eyePos -= cam_speed*(lookAt - eyePos).normalized(); glutPostRedisplay();
		break;
	case 'd':
		m = AngleAxisd(M_PI/36, Vector3d::UnitZ());
		eyePos = (m * (eyePos - lookAt)) + lookAt; glutPostRedisplay();
		// std::cout << "eyePos : " << eyePos << std::endl;
		break;
	case 'a':
		m = AngleAxisd(-M_PI/36, Vector3d::UnitZ());
		eyePos = (m * (eyePos - lookAt)) + lookAt; glutPostRedisplay();
		// std::cout << "eyePos : " << eyePos << std::endl;
		break;
	default:
		break;
	}
}

void mouse_func(int button, int state, int x, int y){

}

std::vector<int> muscles;
int idle_num = 0;
std::vector<std::vector<int> > probableMuscleCollIdx;
std::vector<std::vector<int> > probableTriangleCollIdx;
std::vector<std::vector<int> > probableBoneCollIdx;
std::vector<std::vector<int> > probableTriBoneCollIdx;

// void idle_func() {}
void idle_func() {
	// std::cout << "idle_func taken!" << std::endl;
	if (pauseManual) return;
	Shape model_cur; double delta_t = 0.004;
	model_cur.particles = model.particles;
	// model_cur.particles = model_base.particles;
	for (int i = 0; i < model.particles.size(); ++i) { // update position from velocity
		model_cur.particles[i].pos += 0.9*model.particles[i].vel * delta_t;
		model_cur.particles[i].set_rotation(model.particles[i].omega * model.particles[i].mat3_rot);
	}

	// while (constraint(model_cur.particles)) {
	int step = 2;
	while (step--) {

		// apply_constraint(model_cur.particles);
		int maxVal = 100;
		double vis = 1.0*idle_num; int vist = vis/(2*maxVal); vis -= vist*2*maxVal; if (vis > maxVal) vis = maxVal-(vis-maxVal);// vis/=2.0;
		
		for (int i = 0; i < boneParticles[0].size(); ++i) {
			model_cur.particles[boneParticles[0][i]].pos = model_base.particles[boneParticles[0][i]].pos;
		}
		// humerusRot = AngleAxisd(vis*M_PI/180, -1*Vector3d::UnitY());
		for (int i = 0; i < boneParticles[1].size(); ++i) {
			model_cur.particles[boneParticles[1][i]].pos = (humerusRot * (model_base.particles[boneParticles[1][i]].pos - humerusPivot)) + humerusPivot;
		}
		ulnaRot = AngleAxisd(vis*M_PI/180, -1*Vector3d::UnitY());
		for (int i = 0; i < boneParticles[2].size(); ++i) {
			model_cur.particles[boneParticles[2][i]].pos = (ulnaRot * (model_base.particles[boneParticles[2][i]].pos - ulnaPivot)) + ulnaPivot;
			model_cur.particles[boneParticles[2][i]].pos = (humerusRot * (model_cur.particles[boneParticles[2][i]].pos - humerusPivot)) + humerusPivot;
		}
		/*
		Matrix3d Rotate;
		Rotate = AngleAxisd(vis*M_PI/180, -1*Vector3d::UnitY());
		// int cut = 184, pivot = 188;
		double cutZ = 127.242; Vector3d pivot(9.09645,-13.9159,130.049);
		// if (idle_num == 0) {
		// 	std::cout << idle_num << " " << model_base.particles[cut].pos.z() << std::endl;
		// 	std::cout << model_base.particles[pivot].pos << std::endl;
		// }
		for (int i = 0; i < stuck.size(); ++i) { if (stuck[i] >= model_base.particles.size()) continue;
			if (model_base.particles[stuck[i]].pos.z() < cutZ) {
				model_cur.particles[stuck[i]].pos = (Rotate*(model_base.particles[stuck[i]].pos - pivot)) + pivot;
				model_cur.particles[stuck[i]].set_rotation(Rotate*model_base.particles[stuck[i]].mat3_rot);
			} else {
				model_cur.particles[stuck[i]].pos = model_base.particles[stuck[i]].pos;
				model_cur.particles[stuck[i]].set_rotation(model_base.particles[stuck[i]].mat3_rot);
			}
		}
		for (int i = 0; i < moved.size(); ++i){ if (moved[i] >= model_base.particles.size()) continue;
			if (model_base.particles[moved[i]].pos.z() < cutZ) {
				model_cur.particles[moved[i]].pos = (Rotate*(model_base.particles[moved[i]].pos - pivot)) + pivot;
				model_cur.particles[moved[i]].set_rotation(Rotate*model_base.particles[moved[i]].mat3_rot);
			} else {
				model_cur.particles[moved[i]].pos = model_base.particles[moved[i]].pos;
				model_cur.particles[moved[i]].set_rotation(model_base.particles[moved[i]].mat3_rot);
			}
		}
		*/
		
		std::vector<Vector3d> goals; std::vector<Matrix3d> Rs;
		for (int i = 0; i < model.groups.size(); ++i) goals.push_back(Vector3d::Zero());
		for (int i = 0; i < model.groups.size(); ++i) {
			Matrix3d A = model_cur.particles[i].A;
			A += model_cur.particles[i].mass * model_cur.particles[i].pos * model_base.particles[i].pos.transpose();
			double M = model_cur.particles[i].mass;
			Vector3d c = model_cur.particles[i].mass * model_cur.particles[i].pos;
			Vector3d c_base = model_base.particles[i].mass * model_base.particles[i].pos;
			for (int j = 0; j < model.groups[i].size(); ++j) {
				A += model_cur.particles[model.groups[i][j]].A;
				A += model_cur.particles[model.groups[i][j]].mass * model_cur.particles[model.groups[i][j]].pos * model_base.particles[model.groups[i][j]].pos.transpose();
				M += model_cur.particles[model.groups[i][j]].mass;
				c += model_cur.particles[model.groups[i][j]].mass * model_cur.particles[model.groups[i][j]].pos;
				c_base += model_base.particles[model.groups[i][j]].mass * model_base.particles[model.groups[i][j]].pos;
			}
			c /= M; c_base /= M;
			A -= M * c * c_base.transpose();

			Affine3d affA; affA = A;
			Matrix3d R = affA.rotation();
			// if (i == 10) std::cout << R << std::endl;

			Rs.push_back(R);
			// model_cur.particles[i].set_rotation(R*model_base.particles[i].mat3_rot);
			Vector3d goal = ( R * (model_base.particles[i].pos - c_base) ) + c;
			// goals[i] = goal;
			// goals[i] += 0.8 * goal / (model.groups[i].size() + 1.0);
			goals[i] += 1.0 * goal / (model.groups[i].size() + 1.0);
			for (int j = 0; j < model.groups[i].size(); ++j) {
				Vector3d goal = ( R * (model_base.particles[model.groups[i][j]].pos - c_base) ) + c;
				goals[model.groups[i][j]] += goal / (model.groups[model.groups[i][j]].size() + 1.0);
			}
		}

		for (int i = 0; i < model_cur.particles.size(); ++i) {
			model_cur.particles[i].pos += 1*(goals[i] - model_cur.particles[i].pos);
			model_cur.particles[i].set_rotation(Rs[i]*model_base.particles[i].mat3_rot); /// this might be wrong
			// model_cur.particles[i].set_rotation(Rs[i]); /// this might be the correct way
		}

		// break; // just for observations

		// try to detect collision (sphere to plane/triangle)
		if (handleCollision) {
			for (int s = 0; s < model.particles.size(); ++s) {
				// std::cout << "Hello World! " << omp_get_thread_num() << " of " << omp_get_num_threads() << std::endl;
				// std::cout << idle_num << " " << step << " " << s << std::endl;
				for (int closeTrias = 0; closeTrias < probableMuscleCollIdx[s].size(); ++closeTrias){
					int i = probableTriangleCollIdx[s][closeTrias], muscleCount = probableMuscleCollIdx[s][closeTrias];
					std::vector<Vector3d> tria(3);
					for (int t = 0; t < 3; ++t) { tria[t] = Vector3d::Zero();
						for (int p = 0; p < Muscle_interP[muscleCount][Muscle_interTris[muscleCount][i][t]].size(); ++p) {
							if (onlyTranslation) {
								tria[t] += Muscle_interPw[muscleCount][Muscle_interTris[muscleCount][i][t]][p] * (model_cur.particles[Muscle_interP[muscleCount][Muscle_interTris[muscleCount][i][t]][p]].pos + Muscle_interVec[muscleCount][Muscle_interTris[muscleCount][i][t]][p]);
							} else {
								tria[t] += Muscle_interPw[muscleCount][Muscle_interTris[muscleCount][i][t]][p] * (model_cur.particles[Muscle_interP[muscleCount][Muscle_interTris[muscleCount][i][t]][p]].pos + (model_cur.particles[Muscle_interP[muscleCount][Muscle_interTris[muscleCount][i][t]][p]].mat3_rot*Muscle_interVec[muscleCount][Muscle_interTris[muscleCount][i][t]][p]));
							}
						}
					}						
					Vector3d side1 = tria[1]-tria[0];
					Vector3d side2 = tria[2]-tria[0];
					Vector3d nor = (side1.cross(side2)).normalized();
					// using normal and triangle verticies for collision
		
					if ((nor.dot(model_cur.particles[s].pos+(model_cur.particles[s].size.x()*nor)) - nor.dot(tria[0])) * (nor.dot(model_cur.particles[s].pos-(model_cur.particles[s].size.x()*nor)) - nor.dot(tria[0])) >= 0) continue;
		
					Vector3d planeIntersectionPoint = model_cur.particles[s].pos + (((nor.dot(tria[0]-model_cur.particles[s].pos))/nor.squaredNorm()) * nor);
					planeIntersectionPoint -= tria[0];
					double alpha = ((planeIntersectionPoint.x()*side2.y())-(planeIntersectionPoint.y()*side2.x()))/((side1.x()*side2.y())-(side1.y()*side2.x()));
					double beta = ((planeIntersectionPoint.x()*side1.y())-(planeIntersectionPoint.y()*side1.x()))/((side2.x()*side1.y())-(side2.y()*side1.x()));
					if (alpha < 0.0 || beta < 0.0 || alpha+beta > 1.0) {
						continue;
					} else { // there is a collision
						collisionTriangles[muscleCount][i] = true;
						model_cur.particles[s].pos = planeIntersectionPoint + tria[0] + (model_cur.particles[s].size.x()+0.0001)*nor;
					}
				}
				if (useBones) {
					for (int closeTrias = 0; closeTrias < probableBoneCollIdx[s].size(); ++closeTrias){
						int triaCount = probableTriBoneCollIdx[s][closeTrias], boneCount = probableBoneCollIdx[s][closeTrias];
						std::vector<Vector3d> tria(3);
						tria[0] = boneVerticies[boneCount][boneTriangles[boneCount][triaCount].x()];
						tria[1] = boneVerticies[boneCount][boneTriangles[boneCount][triaCount].y()];
						tria[2] = boneVerticies[boneCount][boneTriangles[boneCount][triaCount].z()];
						Vector3d side1 = tria[1]-tria[0];
						Vector3d side2 = tria[2]-tria[0];
						Vector3d nor = (side1.cross(side2)).normalized();
						// using normal and triangle verticies for collision
			
						if ((nor.dot(model_cur.particles[s].pos+(model_cur.particles[s].size.x()*nor)) - nor.dot(tria[0])) * (nor.dot(model_cur.particles[s].pos-(model_cur.particles[s].size.x()*nor)) - nor.dot(tria[0])) >= 0) continue;
			
						Vector3d planeIntersectionPoint = model_cur.particles[s].pos + (((nor.dot(tria[0]-model_cur.particles[s].pos))/nor.squaredNorm()) * nor);
						planeIntersectionPoint -= tria[0];
						double alpha = ((planeIntersectionPoint.x()*side2.y())-(planeIntersectionPoint.y()*side2.x()))/((side1.x()*side2.y())-(side1.y()*side2.x()));
						double beta = ((planeIntersectionPoint.x()*side1.y())-(planeIntersectionPoint.y()*side1.x()))/((side2.x()*side1.y())-(side2.y()*side1.x()));
						if (alpha < 0.0 || beta < 0.0 || alpha+beta > 1.0) {
							continue;
						} else { // there is a collision
							collBoneTriangles[boneCount][triaCount] = true;
							model_cur.particles[s].pos = planeIntersectionPoint + tria[0] + (model_cur.particles[s].size.x()+0.0001)*nor;
						}
					}	
				}
			}
		}
	}
	
	for (int i = 0; i < model.particles.size(); ++i) {
		// if (i == 25) std::cout << model_cur.particles[i].quat_rot.w() << " " << model_cur.particles[i].quat_rot.x() << " " << model_cur.particles[i].quat_rot.y() << " " << model_cur.particles[i].quat_rot.z() << std::endl;
		
		// model.particles[i].vel = Vector3d::Zero();
		model.particles[i].vel = (model_cur.particles[i].pos - model.particles[i].pos) / delta_t;
		model.particles[i].pos = model_cur.particles[i].pos;

		Matrix3d omega_new = model_cur.particles[i].mat3_rot * model.particles[i].mat3_rot.inverse();
		// if (AngleAxisd(w_new).axis().dot(AngleAxisd(model.particles[i].wQ).axis()) < 0) {
		// 	model.particles[i].wQ = w_new.conjugate();
		// 	model.particles[i].vel = Vector3d::Zero();
		// } else {
		// 	model.particles[i].wQ = w_new;
		// }
		model.particles[i].omega = omega_new;
		model.particles[i].set_rotation(model_cur.particles[i].mat3_rot);
	}
	// if (idle_num > 1000) sleep(1);
	idle_num += 1;
	glutPostRedisplay();
}

int main(int argc, char **argv) {
	for (int i = 0; i < clavicle.size() ; ++i) boneParticles[0].push_back(clavicle[i]);
	for (int i = 0; i < scapula.size() ; ++i) boneParticles[0].push_back(scapula[i]);
	for (int i = 0; i < humerus.size() ; ++i) boneParticles[1].push_back(humerus[i]);
	for (int i = 0; i < ulna.size() ; ++i) boneParticles[2].push_back(ulna[i]);
	for (int i = 0; i < radius.size() ; ++i) boneParticles[2].push_back(radius[i]);
	
	bool remapReq = false;
	if (argc == 2) {
		if (argv[1][0] == 'c') handleCollision = true;
		else if (argv[1][0] == 'l') remapReq = true;
		else if (argv[1][0] == 'b') useBones = true;
	} else if (argc == 3) {
		if (argv[1][0] == 'b') {
			useBones = true;
			if (argv[2][0] == 'c') handleCollision = true;
			else if (argv[2][0] == 'l') remapReq = true;
		} else {handleCollision = true; remapReq = true;}
	} else if (argc == 4) {
		handleCollision = true; remapReq = true; useBones = true;
	}
	if (remapReq) std::cout << "we will reduce muscles" << std::endl;
	if (handleCollision) std::cout << "we will handle collision" << std::endl;

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
	
	model.particles = particles;
	model.groups = groups;

	model_base.particles = particles;
	model_base.groups = groups;
	
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
			// if (i==55){ std::cout << "55 " << temp[0] << std::endl; }
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

		if (handleCollision) {
			std::vector<bool> trias(numTris, false);
			collisionTriangles.push_back(trias);
		}
	}

	std::cout << "All inputs taken" << std::endl;

	if (remapReq) {
		int wholeN; std::cin >> wholeN;
		std::vector<int> particleMap(wholeN, -1);
		for (int i = 0; i < wholeN; ++i) std::cin >> particleMap[i];
		for (int i = 0; i < boneParticles.size(); ++i) {
			std::vector<int> bonePartNew;
			for (int m = 0; m < boneParticles[i].size(); ++m) {
				if (particleMap[boneParticles[i][m]] >= 0) bonePartNew.push_back(particleMap[boneParticles[i][m]]);
			} boneParticles[i] = bonePartNew;
		}
	}

	std::vector<bool> tempo(particles.size(), false);
	for (int i = 0; i < boneParticles.size(); ++i){
		for (int m = 0; m < boneParticles[i].size(); ++m) {
			tempo[boneParticles[i][m]] = true;
		}
	} controlParticles = tempo;

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
			if (handleCollision) {
				std::vector<bool> trias(numTriangles, false);
				collBoneTriangles.push_back(trias);
			}
			while(numTriangles--) {
				std::cin >> vORf >> t1 >> t2 >> t3; triangles.push_back(Vector3i(t1-1,t2-1,t3-1));
			}
			boneTriangles.push_back(triangles);
		}
	}

	if (handleCollision) {
	// if (false){
		int c,probNum;
		for (int i = 0; i < particles.size(); ++i) {
			std::cin >> c >> probNum;
			std::vector<int> muscleIdxs(probNum);
			std::vector<int> triangleIdxs(probNum);
			for (int j = 0; j < probNum; ++j) {
				std::cin >> muscleIdxs[j] >> triangleIdxs[j];
			}
			probableMuscleCollIdx.push_back(muscleIdxs); probableTriangleCollIdx.push_back(triangleIdxs);
			if (useBones) {
				std::cin >> probNum;
				std::vector<int> boneIdxs(probNum);
				std::vector<int> triBoneIdxs(probNum);
				for (int j = 0; j < probNum; ++j) {
					std::cin >> boneIdxs[j] >> triBoneIdxs[j];
				}
				probableBoneCollIdx.push_back(boneIdxs); probableTriBoneCollIdx.push_back(triBoneIdxs);
			}
		}
	}
	std::cout << "All alpha inputs taken!" << std::endl;

	// std::cout << "All alpha inputs taken : " << probableMuscleCollIdx[particles.size()-1][probableMuscleCollIdx[particles.size()-1].size()-1] << std::endl;
	// std::cout << "All alpha inputs taken : " << probableMuscleCollIdx[particles.size()-1].size() << std::endl;

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(600, 50);
	glutInitWindowSize(640, 640);
	glutCreateWindow("Oriented Particles");

	glClearColor(0.3,0.3,0.3,0.1);

	// GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	// GLfloat mat_shininess[] = { 50.0 };
	// GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 };
	// glClearColor (0.0, 0.0, 0.0, 0.0);
	// glShadeModel (GL_SMOOTH);

	// glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	// glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
	// glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	// glEnable(GL_CULL_FACE);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_DEPTH_TEST);

	GLfloat lightpos[] = {10, 10, 10, 1.0};
	glLightfv(GL_LIGHT0, GL_POSITION, lightpos);

	glutDisplayFunc(render);
	glutReshapeFunc(reshape_func);
	glutKeyboardFunc(keyboard_func);
 	glutMouseFunc(mouse_func);
 	glutIdleFunc(idle_func);

	glutMainLoop();

	return 0;
}
