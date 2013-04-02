#include "grid_data.h"


GridData::GridData() :
   mDfltValue(0.0), mMax(0.0,0.0,0.0)
{
}

GridData::GridData(const GridData& orig) :
   mDfltValue(orig.mDfltValue)
{
   mData = orig.mData;
   mMax = orig.mMax;
}

GridData::~GridData() 
{
}

std::vector<double>& GridData::data()
{
   return mData;
}

GridData& GridData::operator=(const GridData& orig)
{
   if (this == &orig)
   {
      return *this;
   }
   mDfltValue = orig.mDfltValue;
   mData = orig.mData;
   mMax = orig.mMax;
   return *this;
}

void GridData::initialize(double dfltValue)
{
   mDfltValue = dfltValue;
   mMax[0] = theCellSize*theDim[0];
   mMax[1] = theCellSize*theDim[1];
   mMax[2] = theCellSize*theDim[2];
   mData.resize(theDim[0]*theDim[1]*theDim[2], false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridData::operator()(int i, int j, int k)
{
   static double dflt = 0;
   dflt = mDfltValue;  // HACK: Protect against setting the default value

   if (i< 0 || j<0 || k<0 || 
       i > theDim[0]-1 || 
       j > theDim[1]-1 || 
       k > theDim[2]-1) return dflt;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];

   return mData[col+row+stack];
}

const double GridData::operator()(int i, int j, int k) const
{
   static double dflt = 0;
   dflt = mDfltValue;  // HACK: Protect against setting the default value

   if (i< 0 || j<0 || k<0 || 
       i > theDim[0]-1 || 
       j > theDim[1]-1 || 
       k > theDim[2]-1) return dflt;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];

   return mData[col+row+stack];
}

void GridData::getCell(const vec3& pt, int& i, int& j, int& k)
{
   vec3 pos = worldToSelf(pt); 
   i = (int) (pos[0]/theCellSize);
   j = (int) (pos[1]/theCellSize);
   k = (int) (pos[2]/theCellSize);   
}

double GridData::interpolate(const vec3& pt)
{

	// TODO: Implement sharper cubic interpolation here.
	
	vec3 pos = worldToSelf(pt);

	int i = (int) (pos[0]/theCellSize);
	int j = (int) (pos[1]/theCellSize);
	int k = (int) (pos[2]/theCellSize);

	double scale = 1.0/theCellSize;  
	double fractx = scale*(pos[0] - i*theCellSize);
	double fracty = scale*(pos[1] - j*theCellSize);
	double fractz = scale*(pos[2] - k*theCellSize);
	assert (fractx < 1.0 && fractx >= 0);
	assert (fracty < 1.0 && fracty >= 0);
	assert (fractz < 1.0 && fractz >= 0);
	//cubic
	/*bool t = (i==0||i==theDim[0]||j==0||j==theDim[1]||k==0||k==theDim[2]);
	if(!t)
	{
		//Y @ low X low Z:
		double yll_tmp1 = (*this)(i-1,j-1,k-1);
		double yll_tmp2 = (*this)(i-1,j,k-1);
		double yll_tmp3 = (*this)(i-1,j+1,k-1);
		double yll_tmp4 = (*this)(i-1,j+2,k-1);
		double delta_yll = yll_tmp3 - yll_tmp2;
		double yll_tmp13 = 0.5f*(yll_tmp3 - yll_tmp1)/theCellSize;
		double yll_tmp24 = 0.5f*(yll_tmp4 - yll_tmp2)/theCellSize;
		double yll_1234 = CUBIC(yll_tmp2,yll_tmp13,yll_tmp24,delta_yll,fracty);
		//Y @ low X middle Z:
		double ylm_tmp1 = (*this)(i-1,j-1,k);
		double ylm_tmp2 = (*this)(i-1,j,k);
		double ylm_tmp3 = (*this)(i-1,j+1,k);
		double ylm_tmp4 = (*this)(i-1,j+2,k);
		double delta_ylm = ylm_tmp3 - ylm_tmp2;
		double ylm_tmp13 = 0.5f*(ylm_tmp3 - ylm_tmp1)/theCellSize;
		double ylm_tmp24 = 0.5f*(ylm_tmp4 - ylm_tmp2)/theCellSize;
		double ylm_1234 = CUBIC(ylm_tmp2,ylm_tmp13,ylm_tmp24,delta_ylm,fracty);
		//Y @ low X high Z:
		double ylh_tmp1 = (*this)(i-1,j-1,k+1);
		double ylh_tmp2 = (*this)(i-1,j,k+1);
		double ylh_tmp3 = (*this)(i-1,j+1,k+1);
		double ylh_tmp4 = (*this)(i-1,j+2,k+1);
		double delta_ylh = ylh_tmp3 - ylh_tmp2;
		double ylh_tmp13 = 0.5f*(ylh_tmp3 - ylh_tmp1)/theCellSize;
		double ylh_tmp24 = 0.5f*(ylh_tmp4 - ylh_tmp2)/theCellSize;
		double ylh_1234 = CUBIC(ylh_tmp2,ylh_tmp13,ylh_tmp24,delta_ylh,fracty);
		//Y @ low X top Z:
		double ylt_tmp1 = (*this)(i-1,j-1,k+2);
		double ylt_tmp2 = (*this)(i-1,j,k+2);
		double ylt_tmp3 = (*this)(i-1,j+1,k+2);
		double ylt_tmp4 = (*this)(i-1,j+2,k+2);
		double delta_ylt = ylt_tmp3 - ylt_tmp2;
		double ylt_tmp13 = 0.5f*(ylt_tmp3 - ylt_tmp1)/theCellSize;
		double ylt_tmp24 = 0.5f*(ylt_tmp4 - ylt_tmp2)/theCellSize;
		double ylt_1234 = CUBIC(ylt_tmp2,ylt_tmp13,ylt_tmp24,delta_ylt,fracty);
		//X @ low X
		double dil = 0.5f*(ylh_1234 - yll_1234)/theCellSize;
		double dip1l = 0.5f*(ylt_1234 - ylm_1234)/theCellSize;
		double delta_pl = ylh_1234-ylm_1234;
		double xl = CUBIC(ylm_1234,dil,dip1l,delta_pl,fractx);

		//Y @ middle X low Z:
		double yml_tmp1 = (*this)(i,j-1,k-1);
		double yml_tmp2 = (*this)(i,j,k-1);
		double yml_tmp3 = (*this)(i,j+1,k-1);
		double yml_tmp4 = (*this)(i,j+2,k-1);
		double delta_yml = yml_tmp3 - yml_tmp2;
		double yml_tmp13 = 0.5f*(yml_tmp3 - yml_tmp1)/theCellSize;
		double yml_tmp24 = 0.5f*(yml_tmp4 - yml_tmp2)/theCellSize;
		double yml_1234 = CUBIC(yml_tmp2,yml_tmp13,yml_tmp24,delta_yml,fracty);
		//Y @ middle X middle Z:
		double ymm_tmp1 = (*this)(i,j-1,k);
		double ymm_tmp2 = (*this)(i,j,k);
		double ymm_tmp3 = (*this)(i,j+1,k);
		double ymm_tmp4 = (*this)(i,j+2,k);
		double delta_ymm = ymm_tmp3 - ymm_tmp2;
		double ymm_tmp13 = 0.5f*(ymm_tmp3 - ymm_tmp1)/theCellSize;
		double ymm_tmp24 = 0.5f*(ymm_tmp4 - ymm_tmp2)/theCellSize;
		double ymm_1234 = CUBIC(ymm_tmp2,ymm_tmp13,ymm_tmp24,delta_ymm,fracty);
		//Y @ middle X high Z:
		double ymh_tmp1 = (*this)(i,j-1,k+1);
		double ymh_tmp2 = (*this)(i,j,k+1);
		double ymh_tmp3 = (*this)(i,j+1,k+1);
		double ymh_tmp4 = (*this)(i,j+2,k+1);
		double delta_ymh = ymh_tmp3 - ymh_tmp2;
		double ymh_tmp13 = 0.5f*(ymh_tmp3 - ymh_tmp1)/theCellSize;
		double ymh_tmp24 = 0.5f*(ymh_tmp4 - ymh_tmp2)/theCellSize;
		double ymh_1234 = CUBIC(ymh_tmp2,ymh_tmp13,ymh_tmp24,delta_ymh,fracty);
		//Y @ middle X top Z:
		double ymt_tmp1 = (*this)(i,j-1,k+2);
		double ymt_tmp2 = (*this)(i,j,k+2);
		double ymt_tmp3 = (*this)(i,j+1,k+2);
		double ymt_tmp4 = (*this)(i,j+2,k+2);
		double delta_ymt = ymt_tmp3 - ymt_tmp2;
		double ymt_tmp13 = 0.5f*(ymt_tmp3 - ymt_tmp1)/theCellSize;
		double ymt_tmp24 = 0.5f*(ymt_tmp4 - ymt_tmp2)/theCellSize;
		double ymt_1234 = CUBIC(ymt_tmp2,ymt_tmp13,ymt_tmp24,delta_ymt,fracty);
		//X @ middle X
		double dim = 0.5f*(ymh_1234 - yml_1234)/theCellSize;
		double dip1m = 0.5f*(ymt_1234 - ymm_1234)/theCellSize;
		double delta_pm = ymh_1234-ymm_1234;
		double xm = CUBIC(ymm_1234,dim,dip1m,delta_pm,fractx);

		//Y @ high X low Z:
		double yhl_tmp1 = (*this)(i+1,j-1,k-1);
		double yhl_tmp2 = (*this)(i+1,j,k-1);
		double yhl_tmp3 = (*this)(i+1,j+1,k-1);
		double yhl_tmp4 = (*this)(i+1,j+2,k-1);
		double delta_yhl = yhl_tmp3 - yhl_tmp2;
		double yhl_tmp13 = 0.5f*(yhl_tmp3 - yhl_tmp1)/theCellSize;
		double yhl_tmp24 = 0.5f*(yhl_tmp4 - yhl_tmp2)/theCellSize;
		double yhl_1234 = CUBIC(yhl_tmp2,yhl_tmp13,yhl_tmp24,delta_yhl,fracty);
		//Y @ high X middle Z:
		double yhm_tmp1 = (*this)(i+1,j-1,k);
		double yhm_tmp2 = (*this)(i+1,j,k);
		double yhm_tmp3 = (*this)(i+1,j+1,k);
		double yhm_tmp4 = (*this)(i+1,j+2,k);
		double delta_yhm = yhm_tmp3 - yhm_tmp2;
		double yhm_tmp13 = 0.5f*(yhm_tmp3 - yhm_tmp1)/theCellSize;
		double yhm_tmp24 = 0.5f*(yhm_tmp4 - yhm_tmp2)/theCellSize;
		double yhm_1234 = CUBIC(yhm_tmp2,yhm_tmp13,yhm_tmp24,delta_yhm,fracty);
		//Y @ high X high Z:
		double yhh_tmp1 = (*this)(i+1,j-1,k+1);
		double yhh_tmp2 = (*this)(i+1,j,k+1);
		double yhh_tmp3 = (*this)(i+1,j+1,k+1);
		double yhh_tmp4 = (*this)(i+1,j+2,k+1);
		double delta_yhh = yhh_tmp3 - yhh_tmp2;
		double yhh_tmp13 = 0.5f*(yhh_tmp3 - yhh_tmp1)/theCellSize;
		double yhh_tmp24 = 0.5f*(yhh_tmp4 - yhh_tmp2)/theCellSize;
		double yhh_1234 = CUBIC(yhh_tmp2,yhh_tmp13,yhh_tmp24,delta_yhh,fracty);
		//Y @ high X top Z:
		double yht_tmp1 = (*this)(i+1,j-1,k+2);
		double yht_tmp2 = (*this)(i+1,j,k+2);
		double yht_tmp3 = (*this)(i+1,j+1,k+2);
		double yht_tmp4 = (*this)(i+1,j+2,k+2);
		double delta_yht = yht_tmp3 - yht_tmp2;
		double yht_tmp13 = 0.5f*(yht_tmp3 - yht_tmp1)/theCellSize;
		double yht_tmp24 = 0.5f*(yht_tmp4 - yht_tmp2)/theCellSize;
		double yht_1234 = CUBIC(yht_tmp2,yht_tmp13,yht_tmp24,delta_yht,fracty);
		//X @ high X
		double dih = 0.5f*(yhh_1234 - yhl_1234)/theCellSize;
		double dip1h = 0.5f*(yht_1234 - yhm_1234)/theCellSize;
		double delta_ph = yhh_1234-yhm_1234;
		double xh = CUBIC(yhm_1234,dih,dip1h,delta_ph,fractx);

		//Y @ top X low Z:
		double ytl_tmp1 = (*this)(i+2,j-1,k-1);
		double ytl_tmp2 = (*this)(i+2,j,k-1);
		double ytl_tmp3 = (*this)(i+2,j+1,k-1);
		double ytl_tmp4 = (*this)(i+2,j+2,k-1);
		double delta_ytl = ytl_tmp3 - ytl_tmp2;
		double ytl_tmp13 = 0.5f*(ytl_tmp3 - ytl_tmp1)/theCellSize;
		double ytl_tmp24 = 0.5f*(ytl_tmp4 - ytl_tmp2)/theCellSize;
		double ytl_1234 = CUBIC(ytl_tmp2,ytl_tmp13,ytl_tmp24,delta_ytl,fracty);
		//Y @ top X middle Z:
		double ytm_tmp1 = (*this)(i+2,j-1,k);
		double ytm_tmp2 = (*this)(i+2,j,k);
		double ytm_tmp3 = (*this)(i+2,j+1,k);
		double ytm_tmp4 = (*this)(i+2,j+2,k);
		double delta_ytm = ytm_tmp3 - ytm_tmp2;
		double ytm_tmp13 = 0.5f*(ytm_tmp3 - ytm_tmp1)/theCellSize;
		double ytm_tmp24 = 0.5f*(ytm_tmp4 - ytm_tmp2)/theCellSize;
		double ytm_1234 = CUBIC(ytm_tmp2,ytm_tmp13,ytm_tmp24,delta_ytm,fracty);
		//Y @ top X high Z:
		double yth_tmp1 = (*this)(i+2,j-1,k+1);
		double yth_tmp2 = (*this)(i+2,j,k+1);
		double yth_tmp3 = (*this)(i+2,j+1,k+1);
		double yth_tmp4 = (*this)(i+2,j+2,k+1);
		double delta_yth = yth_tmp3 - yth_tmp2;
		double yth_tmp13 = 0.5f*(yth_tmp3 - yth_tmp1)/theCellSize;
		double yth_tmp24 = 0.5f*(yth_tmp4 - yth_tmp2)/theCellSize;
		double yth_1234 = CUBIC(yth_tmp2,yth_tmp13,yth_tmp24,delta_yth,fracty);
		//Y @ top X top Z:
		double ytt_tmp1 = (*this)(i+2,j-1,k+2);
		double ytt_tmp2 = (*this)(i+2,j,k+2);
		double ytt_tmp3 = (*this)(i+2,j+1,k+2);
		double ytt_tmp4 = (*this)(i+2,j+2,k+2);
		double delta_ytt = ytt_tmp3 - yht_tmp2;
		double ytt_tmp13 = 0.5f*(ytt_tmp3 - ytt_tmp1)/theCellSize;
		double ytt_tmp24 = 0.5f*(ytt_tmp4 - ytt_tmp2)/theCellSize;
		double ytt_1234 = CUBIC(ytt_tmp2,ytt_tmp13,ytt_tmp24,delta_ytt,fracty);
		//X @ top X
		double dit = 0.5f*(yth_1234 - ytl_1234)/theCellSize;
		double dip1t = 0.5f*(ytt_1234 - ytm_1234)/theCellSize;
		double delta_pt = yth_1234-ytm_1234;
		double xt = CUBIC(ytm_1234,dit,dip1t,delta_pt,fractx);
		//Z
		double di = 0.5*(xh - xl)/theCellSize;
		double dip1 = 0.5*(xt - xm)/theCellSize;
		double delta_p = xh - xm;
		double z = CUBIC(xm,di,dip1,delta_p,fractz);
		return z;
	}*/
	/*
	// LINEAR INTERPOLATION:
	vec3 pos = worldToSelf(pt);

	int i = (int) (pos[0]/theCellSize);
	int j = (int) (pos[1]/theCellSize);
	int k = (int) (pos[2]/theCellSize);

	double scale = 1.0/theCellSize;  
	double fractx = scale*(pos[0] - i*theCellSize);
	double fracty = scale*(pos[1] - j*theCellSize);
	double fractz = scale*(pos[2] - k*theCellSize);
	assert (fractx < 1.0 && fractx >= 0);
	assert (fracty < 1.0 && fracty >= 0);
	assert (fractz < 1.0 && fractz >= 0);
	*/
	
	
	// Y @ low X, low Z:
	double tmp1 = (*this)(i,j,k);
	double tmp2 = (*this)(i,j+1,k);
	// Y @ high X, low Z:
	double tmp3 = (*this)(i+1,j,k);
	double tmp4 = (*this)(i+1,j+1,k);

	// Y @ low X, high Z:
	double tmp5 = (*this)(i,j,k+1);
	double tmp6 = (*this)(i,j+1,k+1);
	// Y @ high X, high Z:
	double tmp7 = (*this)(i+1,j,k+1);
	double tmp8 = (*this)(i+1,j+1,k+1);

	// Y @ low X, low Z
	double tmp12 = LERP(tmp1, tmp2, fracty);
	// Y @ high X, low Z
	double tmp34 = LERP(tmp3, tmp4, fracty);

	// Y @ low X, high Z
	double tmp56 = LERP(tmp5, tmp6, fracty);
	// Y @ high X, high Z
	double tmp78 = LERP(tmp7, tmp8, fracty);

	// X @ low Z
	double tmp1234 = LERP (tmp12, tmp34, fractx);
	// X @ high Z
	double tmp5678 = LERP (tmp56, tmp78, fractx);

	// Z
	double tmp = LERP(tmp1234, tmp5678, fractz);
	return tmp;


}

vec3 GridData::worldToSelf(const vec3& pt) const
{
   vec3 out;
   out[0] = min(max(0.0, pt[0] - theCellSize*0.5), mMax[0]);
   out[1] = min(max(0.0, pt[1] - theCellSize*0.5), mMax[1]);
   out[2] = min(max(0.0, pt[2] - theCellSize*0.5), mMax[2]);
   return out;
}

GridDataX::GridDataX() : GridData()
{
}

GridDataX::~GridDataX()
{
}

void GridDataX::initialize(double dfltValue)
{
   GridData::initialize(dfltValue);
   mMax[0] = theCellSize*(theDim[0]+1);
   mMax[1] = theCellSize*theDim[1];
   mMax[2] = theCellSize*theDim[2];
   mData.resize((theDim[0]+1)*theDim[1]*theDim[2], false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridDataX::operator()(int i, int j, int k)
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (i < 0 || i > theDim[0]) return dflt;

   if (j < 0) j = 0;
   if (j > theDim[1]-1) j = theDim[1]-1;
   if (k < 0) k = 0;
   if (k > theDim[2]-1) k = theDim[2]-1;

   int col = i;
   int row = k*(theDim[0]+1);
   int stack = j*(theDim[0]+1)*theDim[2];
   return mData[stack + row + col];
}

const double GridDataX::operator()(int i, int j, int k) const
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (i < 0 || i > theDim[0]) return dflt;

   if (j < 0) j = 0;
   if (j > theDim[1]-1) j = theDim[1]-1;
   if (k < 0) k = 0;
   if (k > theDim[2]-1) k = theDim[2]-1;

   int col = i;
   int row = k*(theDim[0]+1);
   int stack = j*(theDim[0]+1)*theDim[2];
   return mData[stack + row + col];
}

vec3 GridDataX::worldToSelf(const vec3& pt) const
{   
   vec3 out;
   out[0] = min(max(0.0, pt[0]), mMax[0]);
   out[1] = min(max(0.0, pt[1]-theCellSize*0.5), mMax[1]);
   out[2] = min(max(0.0, pt[2]-theCellSize*0.5), mMax[2]);
   return out;
}

GridDataY::GridDataY() : GridData()
{
}

GridDataY::~GridDataY()
{
}

void GridDataY::initialize(double dfltValue)
{
   GridData::initialize(dfltValue);
   mMax[0] = theCellSize*theDim[0];
   mMax[1] = theCellSize*(theDim[1]+1);
   mMax[2] = theCellSize*theDim[2];
   mData.resize(theDim[0]*(theDim[1]+1)*theDim[2], false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridDataY::operator()(int i, int j, int k)
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (j < 0 || j > theDim[1]) return dflt;

   if (i < 0) i = 0;
   if (i > theDim[0]-1) i = theDim[0]-1;
   if (k < 0) k = 0;
   if (k > theDim[2]-1) k = theDim[2]-1;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];
   return mData[stack + row + col];
}

const double GridDataY::operator()(int i, int j, int k) const
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (j < 0 || j > theDim[1]) return dflt;

   if (i < 0) i = 0;
   if (i > theDim[0]-1) i = theDim[0]-1;
   if (k < 0) k = 0;
   if (k > theDim[2]-1) k = theDim[2]-1;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];
   return mData[stack + row + col];
}

vec3 GridDataY::worldToSelf(const vec3& pt) const
{
   vec3 out;
   out[0] = min(max(0.0, pt[0]-theCellSize*0.5), mMax[0]);
   out[1] = min(max(0.0, pt[1]), mMax[1]);
   out[2] = min(max(0.0, pt[2]-theCellSize*0.5), mMax[2]);
   return out;
}

GridDataZ::GridDataZ() : GridData()
{
}

GridDataZ::~GridDataZ()
{
}

void GridDataZ::initialize(double dfltValue)
{
   GridData::initialize(dfltValue);
   mMax[0] = theCellSize*theDim[0];
   mMax[1] = theCellSize*theDim[1];
   mMax[2] = theCellSize*(theDim[2]+1);
   mData.resize(theDim[0]*theDim[1]*(theDim[2]+1), false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridDataZ::operator()(int i, int j, int k)
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (k < 0 || k > theDim[2]) return dflt;

   if (i < 0) i = 0;
   if (i > theDim[0]-1) i = theDim[0]-1;
   if (j < 0) j = 0;
   if (j > theDim[1]-1) j = theDim[1]-1;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*(theDim[2]+1);

   return mData[stack + row + col];
}

const double GridDataZ::operator()(int i, int j, int k) const
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (k < 0 || k > theDim[2]) return dflt;

   if (i < 0) i = 0;
   if (i > theDim[0]-1) i = theDim[0]-1;
   if (j < 0) j = 0;
   if (j > theDim[1]-1) j = theDim[1]-1;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*(theDim[2]+1);

   return mData[stack + row + col];
}

vec3 GridDataZ::worldToSelf(const vec3& pt) const
{
   vec3 out;
   out[0] = min(max(0.0, pt[0]-theCellSize*0.5), mMax[0]);
   out[1] = min(max(0.0, pt[1]-theCellSize*0.5), mMax[1]);
   out[2] = min(max(0.0, pt[2]), mMax[2]);
   return out;
}



GridDataInt::GridDataInt() :
   mDfltValue(0)
{
}

GridDataInt::~GridDataInt() 
{
}

std::vector<int>& GridDataInt::data()
{
   return mData;
}

GridDataInt& GridDataInt::operator=(const GridDataInt& orig)
{
   if (this == &orig)
   {
      return *this;
   }
   mDfltValue = orig.mDfltValue;
   mData = orig.mData;
   return *this;
}

void GridDataInt::initialize(int dfltValue)
{
   mDfltValue = dfltValue;
   mData.resize(theDim[0]*theDim[1]*theDim[2], false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

int& GridDataInt::operator()(int i, int j, int k)
{
   static int dflt = 0;
   dflt = mDfltValue;  // HACK: Protect against setting the default value

   if (i< 0 || j<0 || k<0 || 
       i > theDim[0]-1 || 
       j > theDim[1]-1 || 
       k > theDim[2]-1) return dflt;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];

   return mData[col+row+stack];
}

const int GridDataInt::operator()(int i, int j, int k) const
{
   static int dflt = 0;
   dflt = mDfltValue;  // HACK: Protect against setting the default value

   if (i< 0 || j<0 || k<0 || 
       i > theDim[0]-1 || 
       j > theDim[1]-1 || 
       k > theDim[2]-1) return dflt;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];

   return mData[col+row+stack];
}