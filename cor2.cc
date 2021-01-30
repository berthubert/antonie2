#include "refgenome2.hh"

#include <iostream>
#include "misc.hh"
#include <set>
#include <algorithm>
#include "dnamisc.hh"
#include <boost/lexical_cast.hpp>
#include <thread>
#include <fstream>
#include "nucstore.hh"
#include <atomic>
#include <mutex>
#include <fstream>
#include <boost/container/small_vector.hpp>
#include <future>
#include <sstream>
#include <sys/prctl.h>
#include <boost/progress.hpp>
#include "ext/flat_hash_map/bytell_hash_map.hpp"

using namespace std;


namespace std {
    template <>
    struct hash<NucleotideStore> {
        size_t operator () (const NucleotideStore& ns) const { return ns.hash(); }
    };
}


template<typename T, size_t X, size_t Y>
struct Matrix
{
  Matrix()
  {
    values = new T[X*Y];
    for(size_t x=0 ; x < X; ++x)
      for(size_t y=0 ; y < Y; ++y)
	values[x + X*y]=T();
  }
  ~Matrix()
  {
    delete[] values;
  }
  T* values; 
  T& operator()(size_t x, size_t y) {
    return values[x +X*y];
  }

  size_t maxX() const { return X-1; } 
  size_t maxY() const { return Y-1; } 
  
};


typedef ska::bytell_hash_map<uint32_t, vector<unsigned int>> nucs_t;

std::atomic<uint32_t> g_processed{0};
void indexChromosome(nucs_t& nucs, const ReferenceGenome& rg, const ReferenceGenome::Chromosome* chromo, uint8_t val, uint32_t unitsize)
{
  cout<<"This is thread tagging with val "<< (int)val<<endl;
  auto numnucs = chromo->chromosome.size();
  for(uint32_t pos = 0 ; pos < numnucs + unitsize; ++pos) {
    auto r = chromo->chromosome.getRange(pos, unitsize);
    if(!r.isCanonical()) {
      auto& place = nucs[r.getRC().getUInt32()];
      if(place.size() < 50)
        place.push_back(pos);
    }
    else {
      auto &place = nucs[r.getUInt32()];
      if(place.size() < 50)
        place.push_back(pos);
    }

    g_processed++;
  }
}



void indexAll(nucs_t* nucs, ReferenceGenome* rg, string* name, uint8_t val, int unitsize)
{
  auto chromo = rg->getChromosome(*name);
  cout<<"Done reading genome1, have "<<rg->numChromosomes()<<" chromosomes, "<<
    rg->numNucleotides()<<" nucleotides"<<endl;
  cout<<"Chromosome '"<< *name<<"' has "<< chromo->chromosome.size()<<" nucleotides"<<endl;
  
  indexChromosome(*nucs, *rg, chromo, val, unitsize);
}

int main(int argc, char**argv)
{
  cout<<"Start reading genome"<<endl;
  if(argc < 2) {
    cerr<<"Syntax: cor2 fna1 chromo1 fna2 chromo2"<<endl;
    return EXIT_FAILURE;
  }
  constexpr int unitsize=16;

  nucs_t nucs1, nucs2;

  ReferenceGenome rg1(argv[1]);
  string name1(argv[2]);
  ReferenceGenome rg2(argv[3]);
  string name2(argv[4]);

  bool pleaseQuit=false;
  std::thread t3([&rg1, &rg2, &pleaseQuit]() {
      auto tot = rg1.numNucleotides() + rg2.numNucleotides();
      while(!pleaseQuit) {
        
        cout << "\r"<< g_processed << " / " << tot<<", " << g_processed*100.0/tot;
        cout.flush();
        sleep(1);
      }
    });

  
  //  void indexAll(nucs_t* nucs, ReferenceGenome* rg, string* name, uint8_t val, int unitsize)
  std::thread t1(indexAll, &nucs1, &rg1, &name1, 1, unitsize);
  std::thread t2(indexAll, &nucs2, &rg2, &name2, 2, unitsize);

  t1.join();
  t2.join();
  pleaseQuit=true;
  t3.join();


  cout<<"\n\n"<<nucs1.size()<<" different kmers"<<endl;
  cout<<nucs2.size()<<" different kmers"<<endl;
  //  uint32_t only1=0, only2=0, both=0, huh=0;

  constexpr int size = 1500;
  Matrix<unsigned int, size, size> m;
  auto chromo1 = rg1.getChromosome(name1), chromo2 = rg2.getChromosome(name2);
  double f1 = 1.0*size / chromo1->chromosome.size();
  double f2 = 1.0*size / chromo2->chromosome.size();

  boost::progress_display show_progress(nucs1.size());

  vector<unsigned int> xmatches(size), ymatches(size);
  
  for(const auto& n : nucs1) {
    const auto& e2 = nucs2[n.first];
    for(const auto& p1 : n.second)
      for(const auto& p2 : e2) {
        m(f1*p1, f2*p2)++;
        xmatches[f1*p1]++;
        ymatches[f2*p2]++;
      }
    ++show_progress;
  }

  ofstream plot("plot");
  for(int x =0 ; x< 1500; ++x)
    for(int y=0; y < 1500; ++y)
      plot << x << " " << y <<" " << m(x,y) << "\n";
  plot.flush();

  ofstream xplot("xplot");
  int counter=0;
  for(const auto& x : xmatches) {
    xplot<<counter<<" "<<x<<"\n";
    ++counter;
  }
  xplot.flush();
  
  ofstream yplot("yplot");
  counter=0;
  for(const auto& y : ymatches) {
    yplot<<counter<<" "<<y<<"\n";
    ++counter;
  }

  yplot.flush();
  
  _exit(0);
}

