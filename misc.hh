#pragma once
#include <stdio.h>
#include <string>
#include <stdint.h>
#include <zlib.h>
#include <vector>
#include <functional>

void chomp(char* line);
char* sfgets(char* p, int num, FILE* fp);
void reverseNucleotides(std::string* nucleotides);
uint64_t filesize(const char* name);
bool stringfgets(FILE* fp, std::string* line);
bool stringfgets(gzFile fp, std::string* line);

/** Rapid estimator of variance and mean of a series of doubles. 
    API compatible with a, sadly, far slower boost::accumulator_set
    doing the same thing*/
class VarMeanEstimator
{
public:
  VarMeanEstimator() : N(0), xTot(0), x2Tot(0) {}
  void operator()(double val) 
  {
    ++N;
    xTot += val;
    x2Tot += val*val;
  }
  bool valid() const
  {
    return N>0;
  }
  friend double mean(const VarMeanEstimator& vme);
  friend double variance(const VarMeanEstimator& vme);
private:
  uint64_t N;
  double xTot;
  double x2Tot;
};

//! extract 'mean' from a VarMeanEstimator
inline double mean(const VarMeanEstimator& vme)
{
  return vme.xTot/vme.N;
}

//! extract 'variance' from a VarMeanEstimator
inline double variance(const VarMeanEstimator& vme) 
{
  return (vme.x2Tot - vme.xTot*vme.xTot/vme.N)/vme.N;
}

std::string compilerVersion();
void reverseNucleotides(std::string* nucleotides);

std::vector<std::string> expandArguments(int argc, char** argv);
void visitAllNgrams(std::function<void(const std::string&)> exec, unsigned int chars, std::string start="");
