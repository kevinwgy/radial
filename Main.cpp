#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <climits>
#include <cfloat>
#include <vector>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <KDTree.h>
#include <rbf_interp.hpp>

using namespace std;

// ------------------------------------------------------------------------
// An instantiation of "Obj" in KDTree.h (cf. PointIn3D in m2c/Vector3D.h)
// ------------------------------------------------------------------------
class Point {
public:
  int id;
  vector<double> x;
public:
  Point() {}
  Point(int i, vector<double> xin) {id = i; x = xin;}
  double val(int i) const {return x[i];}
  double width([[maybe_unused]] int i) const {return 0.0;}
  int pid() const {return id;}
};


// ------------------------------------------------------------------------
void ReadDataFile(char* filename, int dim_in, int dim_out,
                  vector<vector<double> >& Sin, vector<vector<double> >* Sout = NULL);
double norm2(vector<double>& x, vector<double>& y);
bool IsTheSame(string str1, string str2);
// ------------------------------------------------------------------------


// ------------------------------------------------------------------------
// MAIN FUNCTION
// ------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  if(argc != 8) {
    fprintf(stdout,"Usage: radial [input dim] [output dim] [data file] [prediction points file]"
                   " [#points for interpolation] [RBF (MultiQuadric or InverseMultiQuadric or"
                   " ThinPlateSpline or Gaussian)] [output file]\n");
    exit(-1);
  }

  //-------------------------------------
  // Input & output dimensions 
  //-------------------------------------
  const int dim_in  = atoi(argv[1]);
  const int dim_out = atoi(argv[2]);
  if(dim_in<=0 || dim_out<=0) {
    fprintf(stdout,"Error: Input and output dimensions must be positive.\n");
    exit(-1);
  }

  //-------------------------------------
  // Choice of radial basis function (RBF)
  //-------------------------------------
  void (*phi)(int, double[], double, double[]); //a function pointer
  string rbf(argv[6]);
  if(IsTheSame(rbf, "MULTIQUADRIC"))
    phi = MathTools::phi1;
  else if(IsTheSame(rbf, "INVERSEMULTIQUADRIC"))
    phi = MathTools::phi2;
  else if(IsTheSame(rbf, "THINPLATESPLINE"))
    phi = MathTools::phi3;
  else if(IsTheSame(rbf, "GAUSSIAN"))
    phi = MathTools::phi4;
  else {
    fprintf(stdout,"Error: Detected unknown radial basis function (RBF): %s.\n", argv[6]);
    exit(-1);
  }
  fprintf(stdout,"- Using radial basis function (RBF): %s.\n", argv[6]);

  //-------------------------------------
  // Read data file
  //-------------------------------------
  vector<vector<double> > input, output;
  ReadDataFile(argv[3], dim_in, dim_out, input, &output);
  fprintf(stdout,"- Loaded %d data points from %s.\n", (int)input.size(), argv[3]);

  //-------------------------------------
  // Build the tree
  //-------------------------------------
  if(dim_in>5) {
    fprintf(stdout,"Error: Need to generalize the code (slightly) to handle input dim > 5.\n");
    exit(-1);
  }
  int N = input.size();
  assert(N == (int)output.size());
  vector<Point> p;
  for(int i=0; i<N; i++)
    p.push_back(Point(i, input[i]));
  void* tree;
  if(dim_in == 1) //have to do this because template variable "dim" needs to be known at compile-time.
    tree = new KDTree<Point, 1>(N, p.data());
  else if(dim_in == 2)
    tree = new KDTree<Point, 2>(N, p.data());
  else if(dim_in == 3)
    tree = new KDTree<Point, 3>(N, p.data());
  else if(dim_in == 4)
    tree = new KDTree<Point, 4>(N, p.data());
  else {
    assert(dim_in == 5);
    tree = new KDTree<Point, 5>(N, p.data());
  }
  fprintf(stdout,"- Constructed a KDTree (K=%d).\n", dim_in);

  //-------------------------------------
  // Read prediction points file
  //-------------------------------------
  vector<vector<double> > input2;
  ReadDataFile(argv[4], dim_in, 0, input2);
  fprintf(stdout,"- Loaded %d data points (coordinates only) from %s.\n", (int)input2.size(), argv[4]);

  //-------------------------------------
  // Make predictions & write to file
  //-------------------------------------
  int numPoints = atoi(argv[5]); //number of points for interpolation
  assert(numPoints>0);
  int maxCand = numPoints*20;
  Point candidates[maxCand];
  double cutoff = 1.0; //tentative cutoff distance (for individual dims) --- will be adjusted
  fstream file(argv[7], fstream::out);
  if(!file.is_open()) {
    fprintf(stdout, "Error: Cannot open file %s.\n", argv[7]);
    exit(-1); 
  }

  fprintf(stdout,"- Interpolating...\n");
  vector<double> xd(dim_in*numPoints,0);
  vector<double> fd(numPoints,0);
  vector<double> rbf_weight(numPoints,0);

  int pid = 0;
  for(auto&& x : input2) {

    //Write input to file
    for(auto&& xi : x)
      file << fixed << setw(14) << setprecision(10) << xi << "  ";

    //Find points for interpolation
    int nFound = 0, counter = 0; 
    double low_cut = 0.0, high_cut = DBL_MAX;
    while(nFound<numPoints || nFound>maxCand) {
      if(++counter>1000) {
        fprintf(stdout,"*** Error: Cannot find enough sample points for "
                       "interpolation after %d iterations. "
                       "Coords:", counter);
        for(auto&& xi : x)
          fprintf(stdout," %e", xi);
        fprintf(stdout,". Candidates: %d, cutoff = %e.\n", nFound, cutoff);
        exit(-1);
      }


      if(dim_in == 1) //have to do this because template variable "dim" needs to be known at compile-time.
        nFound = static_cast<KDTree<Point, 1>*>(tree)->findCandidatesWithin(x.data(), candidates, maxCand, cutoff);
      else if(dim_in == 2)
        nFound = static_cast<KDTree<Point, 2>*>(tree)->findCandidatesWithin(x.data(), candidates, maxCand, cutoff);
      else if(dim_in == 3)
        nFound = static_cast<KDTree<Point, 3>*>(tree)->findCandidatesWithin(x.data(), candidates, maxCand, cutoff);
      else if(dim_in == 4)
        nFound = static_cast<KDTree<Point, 4>*>(tree)->findCandidatesWithin(x.data(), candidates, maxCand, cutoff);
      else {
        assert(dim_in == 5);
        nFound = static_cast<KDTree<Point, 5>*>(tree)->findCandidatesWithin(x.data(), candidates, maxCand, cutoff);
      }
 

      if(nFound<numPoints) {
        low_cut = std::max(low_cut, cutoff);
        if(high_cut>0.5*DBL_MAX)
          cutoff *= 4.0;
        else
          cutoff = 0.5*(low_cut + high_cut);
      }
      else if(nFound>maxCand) {
        high_cut = std::min(high_cut, cutoff);
        cutoff = 0.5*(low_cut + high_cut);
      }
    }

    //figure out the actual points for interpolation (numPoints)
    vector<pair<double,int> > dist2node;
    for(int i=0; i<nFound; i++)
      dist2node.push_back(make_pair(norm2(candidates[i].x, x), candidates[i].id));
    sort(dist2node.begin(), dist2node.end());
/*
    fprintf(stdout,"(%e, %e). nFound = %d.\n", x[0], x[1], nFound);
    for(auto&& d2n : dist2node) {
      int id = d2n.second;
      fprintf(stdout,"- %d (%e, %e)->%e. d = %e.\n", id, input[id][0], input[id][1], output[id][0], d2n.first);
    }
    exit(-1);
*/
    dist2node.resize(numPoints);

    //prepare to interpolate
    for(int i=0; i<numPoints; i++)
      for(int j=0; j<dim_in; j++)
        xd[dim_in*i+j] = input[dist2node[i].second][j];

    double r0 = dist2node.front().first + dist2node.back().first;  //slightly larger than maximum separation
    double xinterp = 0.0;

    //interpolation (dim per dim for function output)
    for(int i=0; i<dim_out; i++) {
      for(int j=0; j<numPoints; j++)
        fd[j] = output[dist2node[j].second][i];

      MathTools::rbf_weight(dim_in, numPoints, xd.data(), r0, phi, fd.data(), rbf_weight.data());
      MathTools::rbf_interp(dim_in, numPoints, xd.data(), r0, phi, rbf_weight.data(), 1, x.data(), &xinterp);

      file << "  " << fixed << setw(14) << setprecision(10) << xinterp;
    }
    file << endl;

    pid++;
    if(pid%1000 == 0)
      fprintf(stdout,"  o Done with %d points.\n", pid);
  }

  fprintf(stdout,"  o Done.\n");

  file.close();
  fprintf(stdout,"- Saved results in %s.\n", argv[7]);

  //delete the tree
  if(dim_in == 1) //have to do this because template variable "dim" needs to be known at compile-time.
    delete static_cast<KDTree<Point, 1>*>(tree);
  else if(dim_in == 2)
    delete static_cast<KDTree<Point, 2>*>(tree);
  else if(dim_in == 3)
    delete static_cast<KDTree<Point, 3>*>(tree);
  else if(dim_in == 4)
    delete static_cast<KDTree<Point, 4>*>(tree);
  else {
    assert(dim_in == 5);
    delete static_cast<KDTree<Point, 5>*>(tree);
  }
 
  fprintf(stdout,"- Normal termination.\n");
  return 0;
}


// ------------------------------------------------------------------------
// Read data file
// The first [dim_in] columns go into "Sin" (coordinates of points).
// The next [dim_out] columns go into "Sout" (outputs).
// Additional columns, if exist, are ignored.
// [dim_out] can be zero. In this case, "Sout" is not needed. If it is provided,
// it is not modified.
// ------------------------------------------------------------------------
void ReadDataFile(char* filename, int dim_in, int dim_out,
                  vector<vector<double> >& Sin, vector<vector<double> >* Sout)
{
  assert(dim_in>0 && dim_out>=0);
  if(dim_out>0)
    assert(Sout);

  // open file
  fstream file;
  file.open(filename, fstream::in);
  if (!file.is_open()) {
    fprintf(stdout, "Error: Cannot open file %s.\n", filename);
    exit(-1);
  }

  // Start reading the file
  std::string word, line;

  int r;
  vector<double> data_in(dim_in), data_out(dim_out);

  for(r=0; r<INT_MAX; r++) {

    getline(file, line);

    std::istringstream is(line);
    bool done = false;
    for(int i=0; i<dim_in; i++) {
      is >> data_in[i];
      if(is.fail()) {
        done = true;
        break;
      }
    }
    if(done)
      break;
    for(int i=0; i<dim_out; i++) {
      is >> data_out[i];
      if(is.fail()) {
        done = true;
        break;
      }
    }
    if(done)
      break;

    if(r>=(int)Sin.size()) {
      Sin.resize(Sin.size() + 4000, vector<double>(dim_in,0.0));
      if(dim_out>0)
        Sout->resize(Sout->size() + 4000, vector<double>(dim_out,0.0));
    }

    for(int i=0; i<dim_in; i++)
      Sin[r][i] = data_in[i];
    for(int i=0; i<dim_out; i++)
      (*Sout)[r][i] = data_out[i];
  }

  Sin.resize(r);
  if(dim_out>0)
    Sout->resize(r);

  file.close();
}

// ------------------------------------------------------------------------

double norm2(vector<double>& x, vector<double>& y)
{
  assert(x.size()==y.size());
  double d = 0.0, D = 0.0;
  for(int i=0; i<(int)x.size(); i++) {
    d = x[i]-y[i];
    D += d*d;
  }
  return sqrt(D);
}

// ------------------------------------------------------------------------

bool IsTheSame(string str1, string str2)
{
  if(str1.length() != str2.length())
    return false;
  for(unsigned int i = 0; i < str1.length(); ++i) {
    if (tolower(str1[i]) != tolower(str2[i]))
    return false;
  }
  return true;
}

// ------------------------------------------------------------------------

