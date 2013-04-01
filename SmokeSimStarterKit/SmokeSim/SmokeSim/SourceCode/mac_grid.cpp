// Modified by Peter Kutz, 2011 and 2012.

#include "mac_grid.h"
#include "open_gl_headers.h"
#include "camera.h"
#include "custom_output.h"
#include "constants.h"
#include <math.h>
#include <map>
#include <stdio.h>
#undef max
#undef min
#include <fstream>


// Globals:
MACGrid target;


// NOTE: x -> cols, z -> rows, y -> stacks
MACGrid::RenderMode MACGrid::theRenderMode = SHEETS;
bool MACGrid::theDisplayVel = false;

#define FOR_EACH_CELL \
   for(int k = 0; k < theDim[MACGrid::Z]; k++)  \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_CELL_REVERSE \
   for(int k = theDim[MACGrid::Z] - 1; k >= 0; k--)  \
      for(int j = theDim[MACGrid::Y] - 1; j >= 0; j--) \
         for(int i = theDim[MACGrid::X] - 1; i >= 0; i--) 

#define FOR_EACH_FACE \
   for(int k = 0; k < theDim[MACGrid::Z]+1; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]+1; j++) \
         for(int i = 0; i < theDim[MACGrid::X]+1; i++) 

#define FOR_EACH_FACE_X \
	for(int k = 0; k < theDim[MACGrid::Z]; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X]+1; i++) 

#define FOR_EACH_FACE_X1 \
	for(int k = 0; k < theDim[MACGrid::Z]; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 1; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_FACE_Y \
	for(int k = 0; k < theDim[MACGrid::Z]; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]+1; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_FACE_Y1 \
	for(int k = 0; k < theDim[MACGrid::Z]; k++) \
      for(int j = 1; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_FACE_Z \
	for(int k = 0; k < theDim[MACGrid::Z]+1; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_FACE_Z1 \
	for(int k = 1; k < theDim[MACGrid::Z]; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 



MACGrid::MACGrid()
{
   initialize();
}

MACGrid::MACGrid(const MACGrid& orig)
{
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;
}

MACGrid& MACGrid::operator=(const MACGrid& orig)
{
   if (&orig == this)
   {
      return *this;
   }
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;   

   return *this;
}

MACGrid::~MACGrid()
{
}
//level set
void MACGrid::set_liquid_boundary(double (*phi)(const vec3&)) {
   FOR_EACH_CELL
   {
	   vec3 pos = this->getCenter(i,j,k);
		levelset_phi(i,j,k) = phi(pos);
		
	}
}

void MACGrid::initialize_mparticles(double (*phi)(const vec3&)) {
	int seed = 0;
	int seedi = 0; // for determining whether the particle is inside
	FOR_EACH_CELL
	{
		double corner[8];
		vec3 tmp;
		vec3 tmp2;
		float halfsize = 0.5f*theCellSize;
		tmp = this->getCenter(i,j,k);
		tmp2 = vec3(tmp[0] - halfsize,tmp[1] - halfsize,tmp[2] - halfsize);
		corner[0] = interpolate_value(tmp2, levelset_phi);

		tmp2 = vec3(tmp[0] - halfsize,tmp[1] - halfsize,tmp[2]+halfsize);
		corner[1] = interpolate_value(tmp2, levelset_phi);

		tmp2 = vec3(tmp[0] - halfsize,tmp[1]+halfsize,tmp[2] - halfsize);
		corner[2] = interpolate_value(tmp2, levelset_phi);

		tmp2 = vec3(tmp[0] - halfsize,tmp[1] + halfsize,tmp[2] + halfsize);
		corner[3] = interpolate_value(tmp2, levelset_phi);

		tmp2 = vec3(tmp[0] + halfsize,tmp[1] - halfsize,tmp[2] - halfsize);
		corner[4] = interpolate_value(tmp2, levelset_phi);

		tmp2 = vec3(tmp[0] + halfsize,tmp[1] - halfsize, tmp[2] + halfsize);
		corner[5] = interpolate_value(tmp2, levelset_phi);

		tmp2 = vec3(tmp[0] + halfsize,tmp[1] + halfsize,tmp[2] - halfsize);
		corner[6] = interpolate_value(tmp2, levelset_phi);

		tmp2 = vec3(tmp[0] + halfsize,tmp[1] + halfsize,tmp[2] + halfsize);
		corner[7] = interpolate_value(tmp2, levelset_phi);
				
		bool particleCell = false;
		for(int n=0; n<8; n++) {
			if (corner[n] < 3.0 * theCellSize)
				particleCell = true;
		}

		tmp2 = vec3(tmp[0] - halfsize,tmp[1] - halfsize,tmp[2] - halfsize);
		vec3 pos = tmp2;
		for (int p=0; p<8; p++) {	
			float a = randhashf(seed++); float b = randhashf(seed++); float c = randhashf(seed++);
			vec3 mpos = pos + theCellSize * vec3(a,b,c);
			Particle* mp = new Particle(mpos);
			mparticles.push_back(mp);
		}

		float bmin = rmin;
		float bmax = theCellSize;
		float db = bmax - bmin;
		for (int p=0; p<mparticles.size(); p++) {
			Particle* mp = mparticles[p];
			float phi_goal;
			float inside  = randhashf(seedi++);
			if (inside <= 0.5)
				phi_goal = -bmin - db * randhashf(seed++);
			else
				phi_goal = bmin + db * randhash(seed++);
					
			float phi_curr = interpolate_value(mp->pos, levelset_phi);
			vec3 norm;
			interpolate_gradient(norm, mp->pos, levelset_phi);
			norm = norm.Normalize();
			float l = 1.0; // lambda
			mp->pos += l * (phi_goal - phi_curr) * norm;

			// sign of the particle
			phi_curr = interpolate_value(mp->pos, levelset_phi);
			if (phi_curr <= 0)
				mp->s = -1;
			else
				mp->s = 1;

			// radius of the particle
			if (mp->s * phi_curr > rmax)
				mp->r = rmax;
			else if (mp->s * phi_curr < rmin)
				mp->r = rmin;
			else
				mp->r = mp->s * phi_curr;

			//TODO: check if the particle is inside solid boundary
		}
	}
}
void MACGrid::advect_mparticles(float dt) { 
   for(unsigned int p = 0; p < mparticles.size(); ++p) {
		mparticles[p]->pos = trace_rk2(mparticles[p]->pos, dt);
	}
}

vec3 MACGrid::trace_rk2(const vec3& position, float dt) {
   vec3 input = position;
   vec3 velocity = this->getVelocity(input);
   velocity = this->getVelocity(input + 0.5f*dt*velocity);
   input += dt*velocity;
   return input;
}
//end of level set
void MACGrid::reset()
{
   mU.initialize();
   mV.initialize();
   mW.initialize();
   mP.initialize();
   mD.initialize();
   mT.initialize(0.0);

   setUpAMatrix();
}

void MACGrid::initialize()
{
   reset();
}

void MACGrid::updateSources()
{
    // TODO:Set initial values for density, temperature, and velocity.

	//this->mV(0,1,0) = 10;
	this->mV(0,1,0) = 10;
	//this->mW(0,0,1) = 10;

	this->mT(0,0,0) = 10;

	this->mD(0,0,0) = 2;
}

void MACGrid::advectVelocity(double dt)
{
    // TODO: Calculate new velocities and store in target.
	target.mU = mU;
    target.mV = mV;
    target.mW = mW;
	FOR_EACH_FACE_X
	{
		if(i==0||i==theDim[MACGrid::X]) {target.mU(i,j,k) = 0;continue;}
		vec3 p = this->getCenter(i,j,k);
		p[0] = p[0] - theCellSize/2.0f - mU(i,j,k)*dt;
		target.mU(i,j,k) = target.mU.interpolate(p);
		
	}
	FOR_EACH_FACE_Y
	{
		if(j==0||j==theDim[MACGrid::Y]) {target.mV(i,j,k) = 0;continue;}
		vec3 p = this->getCenter(i,j,k);
		p[1] = p[1] - theCellSize/2.0f - mV(i,j,k)*dt;
		target.mV(i,j,k) = target.mV.interpolate(p);
		
	}

	FOR_EACH_FACE_Z
	{
		if(k==0||k==theDim[MACGrid::X]) {target.mW(i,j,k) = 0;continue;}
		vec3 p = this->getCenter(i,j,k);
		p[2] = p[2] - theCellSize/2.0f - mW(i,j,k)*dt;
		target.mW(i,j,k) = target.mW.interpolate(p);
		
	}

	

    // Then save the result to our object.
    mU = target.mU;
    mV = target.mV;
    mW = target.mW;
}

void MACGrid::advectTemperature(double dt)
{
    // TODO: Calculate new temp and store in target.
	target.mT = mT;
	FOR_EACH_CELL
	{
		vec3 p = this->getCenter(i,j,k);
		p -= this->getVelocity(p)*dt;
		target.mT(i,j,k) = target.mT.interpolate(p);
	}
    // Then save the result to our object.
    mT = target.mT;
}

void MACGrid::advectDensity(double dt)
{
    // TODO: Calculate new densitities and store in target.
	target.mD = mD;
	FOR_EACH_CELL
	{
		vec3 p = this->getCenter(i,j,k);
		p -= this->getVelocity(p)*dt;
		target.mD(i,j,k) = target.mD.interpolate(p);
	}
    // Then save the result to our object.
    mD = target.mD;
}

void MACGrid::computeBouyancy(double dt)
{
	// TODO: Calculate bouyancy and store in target.
	target.mV = mV;
	float alpha = 0;
	float beta = 0;
	float Tamb = 0;
	FOR_EACH_FACE_Y
	{
		if(j==0) continue;
		else if (j==theDim[MACGrid::Y]) continue;
		vec3 p = this->getCenter(i,j,k);
		//vec3 tmp = p;
		p[1] = p[1] - theCellSize/2.0f;
		double T;
		//if(p[0]==0||p[0]==theDim[MACGrid::X]*theCellSize) T = Tamb;
		//else if(p[1]==0||p[1] == theDim[MACGrid::Y]*theCellSize) T = Tamb;
		//else if(p[2] == 0 || p[2] == theDim[MACGrid::Z]*theCellSize) T = Tamb;
		T = this->getTemperature(p);
		double D = this->getDensity(p);
		double BouyF = -alpha*D + beta*(T - Tamb);
		
		target.mV(i,j,k) = target.mV.interpolate(p) + BouyF*dt;

	}
   // Then save the result to our object.
   mV = target.mV;
}
double length(vec3 p)
{
	return sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
}
vec3 cross(vec3 v1,vec3 v2)
{
	return vec3(v1[1]*v2[2] - v1[2]*v2[1],v1[2]*v2[0]-v1[0]*v2[2],v1[0]*v2[1]-v1[1]*v2[0]);
}
void MACGrid::computeVorticityConfinement(double dt)
{
	// TODO: Calculate vorticity confinement forces.
	
	GridData omega_i;
	omega_i.initialize();
	GridData omega_j;
	omega_j.initialize();
	GridData omega_k;
	omega_k.initialize();
	FOR_EACH_CELL
	{
		double wjp1 = this->getVelocityZ(this->getCenter(i,j+1,k));
		double wjn1 = this->getVelocityZ(this->getCenter(i,j-1,k));
		double vkp1 = this->getVelocityY(this->getCenter(i,j,k+1));
		double vkn1 = this->getVelocityY(this->getCenter(i,j,k-1));
		double ukp1 = this->getVelocityX(this->getCenter(i,j,k+1));
		double ukn1 = this->getVelocityX(this->getCenter(i,j,k-1));
		double wip1 = this->getVelocityZ(this->getCenter(i+1,j,k));
		double win1 = this->getVelocityZ(this->getCenter(i-1,j,k));
		double vip1 = this->getVelocityY(this->getCenter(i+1,j,k));
		double vin1 = this->getVelocityY(this->getCenter(i-1,j,k));
		double ujp1 = this->getVelocityX(this->getCenter(i,j+1,k));
		double ujn1 = this->getVelocityX(this->getCenter(i,j-1,k));
		omega_i(i,j,k) = (wjp1 - wjn1 - vkp1 + vkn1)/(2*theCellSize);
		omega_j(i,j,k) = (ukp1 - ukn1 - wip1 + win1)/(2*theCellSize);
		omega_k(i,j,k) = (vip1 - vin1 - ujp1 + ujn1)/(2*theCellSize);
		
	}
	GridData F_i;
	F_i.initialize();
	GridData F_j;
	F_j.initialize();
	GridData F_k;
	F_k.initialize();
	///coefficient here!!
	float ita = 2;
	FOR_EACH_CELL
	{
		vec3 wip1 = vec3(omega_i(i+1,j,k),omega_j(i+1,j,k),omega_k(i+1,j,k));
		vec3 win1 = vec3(omega_i(i-1,j,k),omega_j(i-1,j,k),omega_k(i-1,j,k));
		vec3 wjp1 = vec3(omega_i(i,j+1,k),omega_j(i,j+1,k),omega_k(i,j+1,k));
		vec3 wjn1 = vec3(omega_i(i,j-1,k),omega_j(i,j-1,k),omega_k(i,j-1,k));
		vec3 wkp1 = vec3(omega_i(i,j,k+1),omega_j(i,j,k+1),omega_k(i,j,k+1));
		vec3 wkn1 = vec3(omega_i(i,j,k-1),omega_j(i,j,k-1),omega_k(i,j,k-1));
		double deta_omega_i = (length(wip1) - length(win1))/(2*theCellSize);
		double deta_omega_j = (length(wjp1) - length(wjn1))/(2*theCellSize);
		double deta_omega_k = (length(wkp1) - length(wkn1))/(2*theCellSize);
		vec3 deta_omega = vec3(deta_omega_i,deta_omega_j,deta_omega_k);
		vec3 N = deta_omega/(length(deta_omega)+0.0000001);
		
		if(N.Length()<=0.01)
		{
			int i = 0;
		}
		if(deta_omega.Length()>10)
		{
			int i = 0;
		}
		//std::cout<<N[0]<<" "<<N[1]<<" "<<N[2]<<"\n";
		//system("pause");
		vec3 omega = vec3(omega_i(i,j,k),omega_j(i,j,k),omega_k(i,j,k));
		vec3 F = ita*theCellSize*N.Cross(omega);
		//F.Print("force");
	
		//system("pause");
		F_i(i,j,k) = F[0];
		F_j(i,j,k) = F[1];
		F_k(i,j,k) = F[2];
	}
	
	
	
	// Apply the forces to the current velocity and store the result in target.
	target.mU = mU;
	target.mV = mV;
	target.mW = mW;
	FOR_EACH_FACE_X
	{
		if(i==0||i==theDim[MACGrid::X]) {target.mU(i,j,k) = 0;continue;}
		vec3 p = this->getCenter(i,j,k) - vec3(0.5*theCellSize,0,0);
		double Fx = F_i.interpolate(p);
		target.mU(i,j,k) = this->getVelocityX(p) + Fx * dt;
	}

	FOR_EACH_FACE_Y
	{
		if(j==0||j==theDim[MACGrid::Y]) {target.mV(i,j,k) = 0;continue;}
		vec3 p = this->getCenter(i,j,k) - vec3(0,0.5*theCellSize,0);
		double Fy = F_i.interpolate(p);
		target.mV(i,j,k) = this->getVelocityY(p) + Fy * dt;
	}

	FOR_EACH_FACE_Z
	{
		if(k==0||k==theDim[MACGrid::Z]) {target.mW(i,j,k) = 0;continue;}
		vec3 p = this->getCenter(i,j,k) - vec3(0,0,0.5*theCellSize);
		double Fz = F_i.interpolate(p);
		target.mW(i,j,k) = this->getVelocityZ(p) + Fz * dt;
	}

	// Then save the result to our object.
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
}

void MACGrid::addExternalForces(double dt)
{
   computeBouyancy(dt);
   computeVorticityConfinement(dt);
}
bool MACGrid::check_boundary(vec3 p)
{
	if(p[0]<=0||p[0]>=theDim[MACGrid::X]*theCellSize||p[1]<=0||p[1]>=theDim[MACGrid::Y]*theCellSize||p[2]<=0||p[2]>=theDim[MACGrid::Z]*theCellSize)
	return true;
	else return false;
}

bool MACGrid::checkDivergence(GridData d)
{
	double sum = 0;
	for(int i = 0;i<d.data().size();i++)
	{
		sum += d.data()[i];
	}
	//std::cout<<sum<<"\n";
	if(abs(sum)<0.0001) return true;
	else return false;
	//return true;
}
void MACGrid::project(double dt)
{
	// TODO: Solve Ap = d for pressure.
	// 1. Construct d
	GridData d;
	d.initialize();
	std::vector<double> test2;
	/*
	FOR_EACH_FACE_X
	{
		if(i==0||i==theDim[MACGrid::X]) {target.mU(i,j,k) = 0;continue;}
		
		
	}
	FOR_EACH_FACE_Y
	{
		if(j==0||j==theDim[MACGrid::Y]) {target.mV(i,j,k) = 0;continue;}
		
		
	}

	FOR_EACH_FACE_Z
	{
		if(k==0||k==theDim[MACGrid::X]) {target.mW(i,j,k) = 0;continue;}
		
		
	}*/
	
	FOR_EACH_CELL
	{
		double tmp;
		vec3 p = this->getCenter(i,j,k);
		double half = theCellSize/2.0f;
		vec3 xl = vec3(p[0] - half,p[1],p[2]);
		vec3 xh = vec3(p[0] + half,p[1],p[2]);
		vec3 yl = vec3(p[0],p[1] - half,p[2]);
		vec3 yh = vec3(p[0],p[1] + half,p[2]);
		vec3 zl = vec3(p[0],p[1],p[2] - half);
		vec3 zh = vec3(p[0],p[1],p[2] + half);

		double xvl,xvh,yvl,yvh,zvl,zvh;
		if(xl[0]<=0||xl[0]>=theDim[MACGrid::X]*theCellSize)xvl = 0;
		else xvl = this->getVelocityX(xl);
		if(xh[0]<=0||xh[0]>=theDim[MACGrid::X]*theCellSize)xvh = 0;
		else xvh = this->getVelocityX(xh);

		if(yl[0]<=0||yl[0]>=theDim[MACGrid::Y]*theCellSize)yvl = 0;
		else yvl = this->getVelocityY(yl);
		if(yh[0]<=0||yh[0]>=theDim[MACGrid::Y]*theCellSize)yvh = 0;
		else yvh = this->getVelocityY(yh);
	
		if(zl[0]<=0||zl[0]>=theDim[MACGrid::Z]*theCellSize)zvl = 0;
		else zvl = this->getVelocityZ(zl);
		if(zh[0]<=0||zh[0]>=theDim[MACGrid::Z]*theCellSize)zvh = 0;
		else zvh = this->getVelocityZ(zh);

		
		
		double dx_reverse = 1/theCellSize;
		double constant = -theCellSize*theCellSize/dt;
		tmp = constant*dx_reverse*(xvh - xvl + yvh - yvl + zvh - zvl);
		
		d(i,j,k) = tmp;
		
	}
	double sum = 0;
	for(int i = 0;i<d.data().size();i++)
	{
		sum += d.data()[i];
	}
	// 2. Construct A
	setUpAMatrix();
	// 3. Solve for p
	GridData p;
	p.initialize();
	bool re = conjugateGradient(this->AMatrix,p,d,100,0.01);
	if(!re)
	{
		int i = 0;
	}
	// Subtract pressure from our velocity and save in target.
	target.mP = mP;
	target.mU = mU;
	target.mV = mV;
	target.mW = mW;
	int ii = 0;
	FOR_EACH_CELL
	{
		target.mP(i,j,k) = p(i,j,k);
	}
	FOR_EACH_FACE_X1
	{
		if(i==0||i==theDim[MACGrid::X]) {target.mU(i,j,k) = 0;continue;}
		double px1 = target.mP(i,j,k);
		double px0 = target.mP(i-1,j,k);
		target.mU(i,j,k) -= dt*(1/theCellSize)*(px1 - px0); 
	}
	FOR_EACH_FACE_Y1
	{
		if(j==0||j==theDim[MACGrid::Y]) {target.mV(i,j,k) = 0;continue;}
		double py1 = target.mP(i,j,k);
		double py0 = target.mP(i,j-1,k);
 		target.mV(i,j,k) -= dt*(1/theCellSize)*(py1 - py0);
	}

	FOR_EACH_FACE_Z1
	{
		if(k==0||k==theDim[MACGrid::X]) {target.mW(i,j,k) = 0;continue;}
		double pz1 = target.mP(i,j,k);
		double pz0 = target.mP(i,j,k-1);
		target.mW(i,j,k) -= dt*(1/theCellSize)*(pz1 - pz0);
	}
	// Then save the result to our object

	mP = target.mP;
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
	// IMPLEMENT THIS AS A SANITY CHECK: assert (checkDivergence());
	assert (this->checkDivergence(d));
	
}

vec3 MACGrid::getVelocity(const vec3& pt)
{
   vec3 vel;
   vel[0] = getVelocityX(pt); 
   vel[1] = getVelocityY(pt); 
   vel[2] = getVelocityZ(pt); 
   return vel;
}

double MACGrid::getVelocityX(const vec3& pt)
{
   return mU.interpolate(pt);
}

double MACGrid::getVelocityY(const vec3& pt)
{
   return mV.interpolate(pt);
}

double MACGrid::getVelocityZ(const vec3& pt)
{
   return mW.interpolate(pt);
}

double MACGrid::getTemperature(const vec3& pt)
{
   return mT.interpolate(pt);
}

double MACGrid::getDensity(const vec3& pt)
{
   return mD.interpolate(pt);
}

vec3 MACGrid::getCenter(int i, int j, int k)
{
   double xstart = theCellSize/2.0;
   double ystart = theCellSize/2.0;
   double zstart = theCellSize/2.0;

   double x = xstart + i*theCellSize;
   double y = ystart + j*theCellSize;
   double z = zstart + k*theCellSize;
   return vec3(x, y, z);
}

bool MACGrid::isValidCell(int i, int j, int k)
{
	if (i >= theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
		return false;
	}

	if (i < 0 || j < 0 || k < 0) {
		return false;
	}

	return true;
}

void MACGrid::setUpAMatrix() {

	FOR_EACH_CELL {

		int numFluidNeighbors = 0;
		if (i-1 >= 0) {
			AMatrix.plusI(i-1,j,k) = -1;
			numFluidNeighbors++;
		}
		if (i+1 < theDim[MACGrid::X]) {
			AMatrix.plusI(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (j-1 >= 0) {
			AMatrix.plusJ(i,j-1,k) = -1;
			numFluidNeighbors++;
		}
		if (j+1 < theDim[MACGrid::Y]) {
			AMatrix.plusJ(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (k-1 >= 0) {
			AMatrix.plusK(i,j,k-1) = -1;
			numFluidNeighbors++;
		}
		if (k+1 < theDim[MACGrid::Z]) {
			AMatrix.plusK(i,j,k) = -1;
			numFluidNeighbors++;
		}
		// Set the diagonal:
		AMatrix.diag(i,j,k) = numFluidNeighbors;
	}
}







/////////////////////////////////////////////////////////////////////

bool MACGrid::conjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance) {
	// Solves Ap = d for p.

	FOR_EACH_CELL {
		p(i,j,k) = 0.0; // Initial guess p = 0.	
	}

	GridData r = d; // Residual vector.

	GridData z; z.initialize();
	// TODO: Apply a preconditioner here.
	// For now, just bypass the preconditioner:
	z = r;

	GridData s = z; // Search vector;

	double sigma = dotProduct(z, r);

	for (int iteration = 0; iteration < maxIterations; iteration++) {

		double rho = sigma;

		apply(A, s, z);

		double alpha = rho/dotProduct(z, s);

		GridData alphaTimesS; alphaTimesS.initialize();
		multiply(alpha, s, alphaTimesS);
		add(p, alphaTimesS, p);

		GridData alphaTimesZ; alphaTimesZ.initialize();
		multiply(alpha, z, alphaTimesZ);
		subtract(r, alphaTimesZ, r);

		if (maxMagnitude(r) <= tolerance) {
			//PRINT_LINE("PCG converged in " << (iteration + 1) << " iterations.");
			return true;
		}

		// TODO: Apply a preconditioner here.
		// For now, just bypass the preconditioner:
		z = r;

		double sigmaNew = dotProduct(z, r);

		double beta = sigmaNew / rho;

		GridData betaTimesS; betaTimesS.initialize();
		multiply(beta, s, betaTimesS);
		add(z, betaTimesS, s);
		//s = z + beta * s;

		sigma = sigmaNew;
	}

	PRINT_LINE( "PCG didn't converge!" );
	return false;

}

double MACGrid::dotProduct(const GridData & vector1, const GridData & vector2) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		result += vector1(i,j,k) * vector2(i,j,k);
	}

	return result;
}

void MACGrid::add(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) + vector2(i,j,k);
	}

}

void MACGrid::subtract(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) - vector2(i,j,k);
	}

}

void MACGrid::multiply(const double scalar, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = scalar * vector(i,j,k);
	}

}

double MACGrid::maxMagnitude(const GridData & vector) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		if (abs(vector(i,j,k)) > result) result = abs(vector(i,j,k));
	}

	return result;
}

void MACGrid::apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL { // For each row of the matrix.

		double diag = 0;
		double plusI = 0;
		double plusJ = 0;
		double plusK = 0;
		double minusI = 0;
		double minusJ = 0;
		double minusK = 0;

		diag = matrix.diag(i,j,k) * vector(i,j,k);
		if (isValidCell(i+1,j,k)) plusI = matrix.plusI(i,j,k) * vector(i+1,j,k);
		if (isValidCell(i,j+1,k)) plusJ = matrix.plusJ(i,j,k) * vector(i,j+1,k);
		if (isValidCell(i,j,k+1)) plusK = matrix.plusK(i,j,k) * vector(i,j,k+1);
		if (isValidCell(i-1,j,k)) minusI = matrix.plusI(i-1,j,k) * vector(i-1,j,k);
		if (isValidCell(i,j-1,k)) minusJ = matrix.plusJ(i,j-1,k) * vector(i,j-1,k);
		if (isValidCell(i,j,k-1)) minusK = matrix.plusK(i,j,k-1) * vector(i,j,k-1);

		result(i,j,k) = diag + plusI + plusJ + plusK + minusI + minusJ + minusK;
	}

}




/////////////////////////////////////////////////////////////////////

void MACGrid::saveSmoke(const char* fileName) {
	std::ofstream fileOut(fileName);
	if (fileOut.is_open()) {
		FOR_EACH_CELL {
			fileOut << mD(i,j,k) << std::endl;
		}
		fileOut.close();
	}
}





/////////////////////////////////////////////////////////////////////

void MACGrid::draw(const Camera& c)
{   
   drawWireGrid();
   if (theDisplayVel) drawVelocities();   
   if (theRenderMode == CUBES) drawSmokeCubes(c);
   else drawSmoke(c);
}

void MACGrid::drawVelocities()
{
   // Draw line at each center
   //glColor4f(0.0, 1.0, 0.0, 1.0); // Use this if you want the lines to be a single color.
   glBegin(GL_LINES);
      FOR_EACH_CELL
      {
         vec3 pos = getCenter(i,j,k);
         vec3 vel = getVelocity(pos);
		 vel = vel.Normalize();
         if (vel.Length() > 0.0001)
         {
		   // Un-comment the line below if you want all of the velocity lines to be the same length.
           //vel.Normalize();
           vel *= theCellSize/2.0;
           vel += pos;
		   glColor4f(1.0, 1.0, 0.0, 1.0);
           glVertex3dv(pos.n);
		   glColor4f(0.0, 1.0, 0.0, 1.0);
           glVertex3dv(vel.n);
         }
      }
   glEnd();
}

vec4 MACGrid::getRenderColor(int i, int j, int k)
{

	// Modify this if you want to change the smoke color, or modify it based on other smoke properties.
    double value = mD(i, j, k); 
    return vec4(0, 1.0, 1.0, value);

}

vec4 MACGrid::getRenderColor(const vec3& pt)
{

	// TODO: Modify this if you want to change the smoke color, or modify it based on other smoke properties.
    double value = getDensity(pt); 
    return vec4(0, 1.0, 1.0, value);

}

void MACGrid::drawZSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double startk = back - stepsize;
   double endk = 0;
   double stepk = -theCellSize;

   if (!backToFront)
   {
      startk = 0;
      endk = back;   
      stepk = theCellSize;
   }

   for (double k = startk; backToFront? k > endk : k < endk; k += stepk)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double i = 0.0; i <= right; i += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double i = right; i >= 0.0; i -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawXSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double starti = right - stepsize;
   double endi = 0;
   double stepi = -theCellSize;

   if (!backToFront)
   {
      starti = 0;
      endi = right;   
      stepi = theCellSize;
   }

   for (double i = starti; backToFront? i > endi : i < endi; i += stepi)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double k = 0.0; k <= back; k += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double k = back; k >= 0.0; k -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawSmoke(const Camera& c)
{
   vec3 eyeDir = c.getBackward();
   double zresult = fabs(Dot(eyeDir, vec3(1,0,0)));
   double xresult = fabs(Dot(eyeDir, vec3(0,0,1)));
   //double yresult = fabs(Dot(eyeDir, vec3(0,1,0)));

   if (zresult < xresult)
   {      
      drawZSheets(c.getPosition()[2] < 0);
   }
   else 
   {
      drawXSheets(c.getPosition()[0] < 0);
   }
}

void MACGrid::drawSmokeCubes(const Camera& c)
{
   std::multimap<double, MACGrid::Cube, std::greater<double> > cubes;
   FOR_EACH_CELL
   {
      MACGrid::Cube cube;
      cube.color = getRenderColor(i,j,k);
      cube.pos = getCenter(i,j,k);
      cube.dist = DistanceSqr(cube.pos, c.getPosition());
      cubes.insert(make_pair(cube.dist, cube));
   } 

   // Draw cubes from back to front
   std::multimap<double, MACGrid::Cube, std::greater<double> >::const_iterator it;
   for (it = cubes.begin(); it != cubes.end(); ++it)
   {
      drawCube(it->second);
   }
}

void MACGrid::drawWireGrid()
{
   // Display grid in light grey, draw top & bottom

   double xstart = 0.0;
   double ystart = 0.0;
   double zstart = 0.0;
   double xend = theDim[0]*theCellSize;
   double yend = theDim[1]*theCellSize;
   double zend = theDim[2]*theCellSize;

   glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
      glDisable(GL_LIGHTING);
      glColor3f(0.25, 0.25, 0.25);

      glBegin(GL_LINES);
      for (int i = 0; i <= theDim[0]; i++)
      {
         double x = xstart + i*theCellSize;
         glVertex3d(x, ystart, zstart);
         glVertex3d(x, ystart, zend);

         glVertex3d(x, yend, zstart);
         glVertex3d(x, yend, zend);
      }

      for (int i = 0; i <= theDim[2]; i++)
      {
         double z = zstart + i*theCellSize;
         glVertex3d(xstart, ystart, z);
         glVertex3d(xend, ystart, z);

         glVertex3d(xstart, yend, z);
         glVertex3d(xend, yend, z);
      }

      glVertex3d(xstart, ystart, zstart);
      glVertex3d(xstart, yend, zstart);

      glVertex3d(xend, ystart, zstart);
      glVertex3d(xend, yend, zstart);

      glVertex3d(xstart, ystart, zend);
      glVertex3d(xstart, yend, zend);

      glVertex3d(xend, ystart, zend);
      glVertex3d(xend, yend, zend);
      glEnd();
   glPopAttrib();

   glEnd();
}

#define LEN 0.5
void MACGrid::drawFace(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);
      glEnd();
   glPopMatrix();
}

void MACGrid::drawCube(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0, -1.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN, -LEN);         

         glNormal3d( 0.0,  0.0, -0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN,  LEN, -LEN);
         glVertex3d( LEN,  LEN, -LEN);
         glVertex3d( LEN, -LEN, -LEN);

         glNormal3d(-1.0,  0.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d(-LEN,  LEN,  LEN);
         glVertex3d(-LEN,  LEN, -LEN);

         glNormal3d( 0.0, 1.0,  0.0);
         glVertex3d(-LEN, LEN, -LEN);
         glVertex3d(-LEN, LEN,  LEN);
         glVertex3d( LEN, LEN,  LEN);
         glVertex3d( LEN, LEN, -LEN);

         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);

         glNormal3d(1.0,  0.0,  0.0);
         glVertex3d(LEN, -LEN, -LEN);
         glVertex3d(LEN, -LEN,  LEN);
         glVertex3d(LEN,  LEN,  LEN);
         glVertex3d(LEN,  LEN, -LEN);
      glEnd();
   glPopMatrix();
}
