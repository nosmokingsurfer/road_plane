#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <random>
#include <ctime>
#include <fstream>
#include <chrono>

using namespace std;

typedef vector<vector<float>> Matrix;

struct Point
{
  float x;
  float y;
  float z;

  Point& normalize()
  {
    float sum = powf(x*x + y*y + z*z, 0.5);

    this->x /= sum;
    this->y /= sum;
    this->z /= sum;

    return *this;
  }
};

struct PlaneModel
{
  float A;
  float B;
  float C;
  float D;

  PlaneModel(){}

  inline bool operator== (const PlaneModel& rhs)
  {
    float eps = 0.00001;
    return ((fabs(A - rhs.A) < eps) && (fabs(B - rhs.B) < eps) && (fabs(C - rhs.C) < eps) && (fabs(D - rhs.D) < eps));
  }

  PlaneModel(const float a, const float b, const float c, const float d)
  {
    A = a;
    B = b;
    C = c;
    D = d;
  }

  void print()
  {
    cout << A << " " << B << " " << C << " " << D << endl;
  }

  float dist(const Point& p) const
  {
    float result = {0};

    result = (fabs(A*p.x +B*p.y + C*p.z + D))/(powf(A*A + B*B + C*C, 0.5));

    return result;
  }

  bool isValid()
  {
    if((A == 0)&&(B == 0)&&(C == 0)&&(D == 0))
      return false;
    else
      return true;
  }
};

Point operator*(const Matrix& lhs, const Point& rhs)
{
  return{ rhs.x*lhs[0][0] + rhs.y*lhs[0][1] + rhs.z*lhs[0][2],
          rhs.x*lhs[1][0] + rhs.y*lhs[1][1] + rhs.z*lhs[1][2],
          rhs.x*lhs[2][0] + rhs.y*lhs[2][1] + rhs.z*lhs[2][2] };
}

void printPoint(const Point& p)
{
  cout << p.x << " " << p.y << " " << p.z << endl;
}

void printMatrix(const Matrix& m)
{
  for(const auto& row: m)
  {
    for(const auto& col: row)
    {
      cout << col << " ";
    }
    cout << endl;
  }
}

Point getCOM(const vector<Point>& points)
{
    Point result = {0,0,0};

    for(const auto& p: points)
    {
      result.x += p.x;
      result.y += p.y;
      result.z += p.z;
    }

    auto N = points.size();

    result.x /= N;
    result.y /= N;
    result.z /= N;

    return result;
}

Point getCOM(const vector<Point>& points, const set<int> indices)
{
  Point result = {0,0,0};

  if(indices.size() > 0)
  { 
    for(const auto& i: indices)
    {
      result.x += points[i].x;
      result.y += points[i].y;
      result.z += points[i].z;
    }

    auto N = indices.size();

    result.x /= N;
    result.y /= N;
    result.z /= N;
  }
  else
  {
    result = getCOM(points);
  }

  return result;
}

void getMoments(const vector<Point>& points, Point& com, Matrix& inertia)
{
  inertia = {{0,0,0},{0,0,0},{0,0,0}};

  com = getCOM(points);

  for(auto p: points)
  {
    p.x -= com.x;
    p.y -= com.y;
    p.z -= com.z;

    inertia[0][0] += pow(p.x, 2);
    inertia[0][1] += p.x*p.y;
    inertia[0][2] += p.x*p.z;

    inertia[1][0] += p.x*p.y;
    inertia[1][1] += pow(p.y, 2);
    inertia[1][2] += p.y*p.z;

    inertia[2][0] += p.x*p.z;
    inertia[2][1] += p.y*p.z;
    inertia[2][2] += pow(p.z, 2);
  }
}

void getMoments(const vector<Point>& points, const set<int> indices, Point& com, Matrix& inertia)
{
   inertia = {{0,0,0},{0,0,0},{0,0,0}};

   float x = 0;
   float y = 0;
   float z = 0;

   if(indices.size() > 0)
   { 
      com = getCOM(points, indices);
      for(const auto& i: indices)
      {
         x = points[i].x - com.x;
         y = points[i].y - com.y;
         z = points[i].z - com.z;

         inertia[0][0] += pow(x, 2);
         inertia[0][1] += x*y;
         inertia[0][2] += x*z;

         inertia[1][0] += x*y;
         inertia[1][1] += pow(y, 2);
         inertia[1][2] += y*z;

         inertia[2][0] += x*z;
         inertia[2][1] += y*z;
         inertia[2][2] += pow(z, 2);
      }
   }
   else
   {
      getMoments(points, com, inertia);
   }
}

vector<vector<float>> getAdjointMatrix(const vector<vector<float>>& m)
{
  vector<vector<float>> result = {{0,0,0},{0,0,0},{0,0,0}};

  float a = m[0][0];
  float b = m[0][1];
  float c = m[0][2];

  float d = m[1][0];
  float e = m[1][1];
  float f = m[1][2];

  float g = m[2][0];
  float h = m[2][1];
  float i = m[2][2];

  result[0][0] = e*i - h*f;
  result[0][1] = h*c - b*i;
  result[0][2] = b*f - e*c;

  result[1][0] = g*f - d*i;
  result[1][1] = a*i - g*c;
  result[1][2] = d*c - a*f;

  result[2][0] = d*h - g*e;
  result[2][1] = g*b - a*h;
  result[2][2] = a*e - d*b;

  return result;
}

float checkModelQuality(const float thresh, const vector<Point>& cloud, const PlaneModel& plane)
{
  float result = 0;

  for(const auto& p: cloud)
  {
    result += plane.dist(p);
  }

  return result; 
}

set<int> getMayBeInliers(const vector<Point>& cloud, const int num)
{
  auto size = cloud.size();

  std::mt19937 mersenne(static_cast<std::mt19937::result_type>(chrono::high_resolution_clock::now().time_since_epoch().count()));

  std::uniform_int_distribution<> die(0, static_cast<int>(size) - 1);

  set<int> result;

  if(cloud.size() < num)
    return result;

  while (result.size() < num)
  {
    result.insert(die(mersenne));
  }

  return result;
}

PlaneModel fitModel(const vector<Point>& data, const set<int> indices = {})
{
  PlaneModel model = {0,0,0,0};

  Point com;
  Matrix inertia;
  getMoments(data, indices, com, inertia);
  Matrix adjoined = getAdjointMatrix(inertia);

  Point eigenVector = { 1,2,3 };

  for (auto i = 0; i < 5; i++)
  {
    eigenVector = adjoined*eigenVector;
    eigenVector.normalize();
  }

  model.A = eigenVector.x;
  model.B = eigenVector.y;
  model.C = eigenVector.z;
  model.D = -model.A*com.x - model.B*com.y - model.C*com.z;

  return model;
}

PlaneModel RANSAC(const vector<Point>& data, const float thresh, const float n = 0.55, const int k = 10)
{
  PlaneModel bestFit = {0,0,0,0};

  bool bestFirInitialized = false;

  auto size = data.size();

  int iterations = 0;
  float bestError = 1.0f*size;
  float prevBestError = 2* bestError;

  while (iterations < k)
  {
    set<int> mayBeInliers = getMayBeInliers(data, (int)(size*n));

    PlaneModel mayBeModel = fitModel(data, mayBeInliers);

    if(!bestFirInitialized)
      bestFit = mayBeModel;
      bestFirInitialized = true;

    set<int> alsoInliers = {};

    for(auto i = 0; i < size; i++)
    {
      if(!(mayBeInliers.find(i) != mayBeInliers.end()))
      {
        if(mayBeModel.dist(data[i]) < thresh)
        {
          alsoInliers.insert(i);
        }
      }
    }



    if(alsoInliers.size() + mayBeInliers.size() >= size/2)
    {
      alsoInliers.insert(begin(mayBeInliers), end(mayBeInliers));
      PlaneModel betterModel = fitModel(data, alsoInliers);

      float thisErr = checkModelQuality(thresh, data, betterModel);

      if(thisErr < bestError)
      {
        bestFit = betterModel;
        prevBestError = bestError;
        bestError = thisErr;
      }
    }
    else if (fabs(prevBestError - bestError) < 0.0001)
      break;

    iterations++;
  }

  return bestFit;
}

void readData(float& p, int& N, vector<Point>& data)
{
  cin >> p;
  cin >> N;

  data.resize(N,{0,0,0});

  for(auto& p: data)
  {
    cin >> p.x;
    cin >> p.y;
    cin >> p.z;
  }
}

void readData(float& p, int& N, vector<Point>& data, const string& filename)
{
  ifstream myFile;
  myFile.open(filename);
  
  myFile >> p;
  myFile >> N;

  data.resize(N, { 0,0,0 });

  for (auto& p : data)
  {
    myFile >> p.x;
    myFile >> p.y;
    myFile >> p.z;
  }

  myFile.close();
}

void solve(const float p, const int N, const vector<Point> data)
{
  PlaneModel solution;

  if(data.size() >= 6)
  {
    solution = RANSAC(data, p, 0.5f, 30);  
  }
  else if(data.size() >=3)
  {
    solution = fitModel(data);
  }


  solution.print();
}

int main() {
  float p = 0;
  int N = 0;
  vector<Point> data;

  
  
  //for(auto test: {"test_0.txt", "test_1.txt", "test_2.txt", "test_3.txt", "sdc_point_cloud.txt"})
  //{ 
  //  readData(p, N, data, test);
  //  solve(p,N,data); 
  //}

  readData(p,N,data);
  solve(p,N,data);


  return 0;
}