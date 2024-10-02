#include <KDTree.h>
#include <rbf_interp.hpp>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <limit.h>
#include <vector>
#include <cassert>

using namespace std;

// An instantiation of "Obj" in KDTree.h (cf. PointIn3D in m2c/Vector3D.h)
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


void ReadDataFile(char* filename, int dim_in, int dim_out,
                  vector<vector<double> >& Sin, vector<vector<double> >& Sout);


int main(int argc, char* argv[])
{
  if(argc != 6) {
    fprintf(stdout,"Usage: radial [input dim] [output dim] [data file] [prediction points file] [output file]\n");
    exit(-1);
  }

  int dim_in  = atoi(argc[1]);
  int dim_out = atoi(argv[2]);
  if(dim_in<=0 || dim_out<=0) {
    fprintf(stdout,"Error: Input and output dimensions must be positive.\n");
    exit(-1);
  }

  
  //-------------------------------------
  // Read data file
  //-------------------------------------
  vector<vector<double> > input, output;
  ReadDataFile(argv[3], dim_in, dim_out, input, output);

  //-------------------------------------
  // Build the tree
  //-------------------------------------
  int N = input.size();
  assert(N == (int)output.size());
  vector<Point> p;
  for(int i=0; i<N; i++)
    p.push_back(Point(i, input[i]));

  KDTree<Point, dim_in> tree(N, p.data());

  //-------------------------------------
  // Read prediction points file
  //-------------------------------------
  vector<vector<double> > predpts, dummy; //dummy is not used
  ReadDataFile(argv[4], dim_in, 0, predpts, dummy);

  //-------------------------------------
  // Make predictions
  //-------------------------------------
  I AM HERE
  return 0;
}



void ReadDataFile(char* filename, int dim_in, int dim_out, vector<vector<double> >& Sin, vector<vector<double> >& Sout)
{
  assert(dim_in>0);

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
      Sin.resize(Sin.size() + initial_size, vector<double>(dim_in,0.0));
      if(dim_out>0)
        Sout.resize(Sout.size() + initial_size, vector<double>(dim_out,0.0));
    }

    for(int i=0; i<dim_in; i++)
      Sin[r][i] = data_in[i];
    for(int i=0; i<dim_out; i++)
      Sout[r][i] = data_out[i];
  }

  Sin.resize(r);
  if(dim_out>0)
    Sout.resize(r);
}



