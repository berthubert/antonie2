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
#include "bzlib.h"

using namespace std;

/* Hash all segments in their canonical form 
   Store a vector where these hashes appear
   To make the uniqueness stats, go through all those positions to check the _actual_ different k-mers they represent
   Emit their individual counts

   Memory usage: we end up storing a 4 byte position for each of 3.5 billion k-mers
*/



struct HashStat
{
  std::mutex* m;
  boost::container::small_vector<uint32_t,2> pos;
};

unsigned int g_unitsize=16;

class HashCollector
{
public:
  HashCollector()
  {
    cout<<"Creating the hashbins"<<endl;
    d_hashes.resize(d_hashsize);
    for(auto & h : d_hashes)
      h.m = new std::mutex;
    cout<<"Done"<<endl;
  }

  ~HashCollector()
  {
    for(auto & h : d_hashes)
      delete h.m;
  }
  vector<HashStat> d_hashes;
  const unsigned int d_hashsize=1<<16;

  uint32_t count(const NucleotideStore& ns, const ReferenceGenome& rg) const;
  vector<pair<uint32_t,bool>> getPositions(const NucleotideStore& ns, const ReferenceGenome& rg, uint32_t before=std::numeric_limits<uint32_t>::max()) const;
  void add(const NucleotideStore& ns, uint32_t pos);
  
} g_hashes;

NucleotideStore g_allA, g_allC, g_allG, g_allT;

void HashCollector::add(const NucleotideStore& stretch, uint32_t pos)
{
  uint32_t h;
  if(!stretch.isCanonical())
    h = stretch.getRC().hash() % d_hashsize;
  else
    h = stretch.hash() % d_hashsize;
    
  //  cout<<"Storing '"<<stretch<<"' at pos, h="<<h<<endl;
  std::lock_guard<std::mutex> l(*d_hashes[h].m);
  if(stretch==g_allA || stretch==g_allC || stretch == g_allG || stretch == g_allT)
    if(d_hashes[h].pos.size() >= 100)
      return;
  d_hashes[h].pos.push_back(pos);  
}

uint32_t HashCollector::count(const NucleotideStore& stretch, const ReferenceGenome& rg) const
{
  uint32_t h;

  if(stretch.isCanonical()) {
    h = stretch.hash() % d_hashsize;
  }
  else
    h = stretch.getRC().hash() % d_hashsize;

  uint32_t ret=0;
  std::lock_guard<std::mutex> l(*d_hashes[h].m);
  
  for(const auto& e : d_hashes[h].pos) {
    auto cmp = rg.getRange(e, g_unitsize);
    if(cmp==stretch)
      ++ret;
    else if(cmp.getRC() == stretch) {
      ++ret;
    }
  }
  return ret;
}

vector<pair<uint32_t,bool>> HashCollector::getPositions(const NucleotideStore& stretch, const ReferenceGenome& rg, uint32_t before) const
{
  vector<pair<uint32_t,bool>> ret;

  uint32_t h;

  if(stretch.isCanonical()) {
    h = stretch.hash() % d_hashsize;
  }
  else
    h = stretch.getRC().hash() % d_hashsize;


  std::lock_guard<std::mutex> l(*d_hashes[h].m);
  //  cout<<"Lookup "<<stretch<<", h="<<h<<", have "<<d_hashes[h].pos.size()<<" candidates"<<endl;
  for(const auto& e : d_hashes[h].pos) {
    if(e >= before)
      continue;
    auto cmp = rg.getRange(e, g_unitsize);
    if(cmp==stretch) 
      ret.push_back({e,false});
    else if (cmp.getRC()==stretch)
      ret.push_back({e,true});
    else {
      //      cout<<"Candidate "<<cmp<<" matched neither way"<<endl;
    }
  }
  return ret;
}

void indexChr(ReferenceGenome::Chromosome* chromosome, std::string name)
{
  if(chromosome->fullname.find("primary")==string::npos) {
    cout<<"Not indexing "<<chromosome->fullname<<endl;
    return;
  }
  
  prctl(PR_SET_NAME, string("Indexing "+name).c_str());
  auto size=chromosome->chromosome.size();
  cout<<"Starting index of '"<<name<<"' with "<<size<<" nucleotides"<<endl;
  for(size_t pos =0; pos < size - g_unitsize; ++pos) {
    auto stretch = chromosome->chromosome.getRange(pos,  g_unitsize);
    g_hashes.add(stretch, chromosome->offset+pos);

  }
  cout<<"Done with index of '"<<name<<"' with "<<size<<" nucleotides"<<endl;
}


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

set<NucleotideStore> getUniNucs(const ReferenceGenome::Chromosome& chromo, uint32_t start, uint32_t len)
{
  set<NucleotideStore> ret;
  for(uint32_t pos = start ; pos < start + len; ++pos) {
    try {
      auto r = chromo.chromosome.getRange(pos, g_unitsize);
      if(!r.isCanonical())
        ret.insert(r.getRC());
      else
        ret.insert(r);
    }
    catch(std::exception& e) {
      cerr<<e.what()<<endl;
    }
  }
  return ret;
}

template<class InputIt1, class InputIt2>
size_t set_intersection_size(InputIt1 first1, InputIt1 last1,
                               InputIt2 first2, InputIt2 last2)
{
  size_t count = 0;
  while (first1 != last1 && first2 != last2) {
    if (*first1 < *first2) {
      ++first1;
    } else  {
      if (!(*first2 < *first1)) {
        ++count;
        *first1++;
      }
      ++first2;
    }
  }
  return count;
}


void doKMERMap(const ReferenceGenome& rg1, const std::string& name1, const ReferenceGenome& rg2, const std::string& name2)
{
  struct Count
  {
    uint32_t xnucs{0}, ynucs{0}, overnucs{0};
  };
  constexpr int numchunks=1000;
  Matrix<Count, numchunks, numchunks> m;

  auto chromo1 = rg1.getChromosome(name1);
  auto chromo2 = rg2.getChromosome(name2);
  
  auto numnucs1 = chromo1->chromosome.size();
  auto numnucs2 = chromo2->chromosome.size();

  unsigned int chunksize = std::max(chromo1->chromosome.size(), chromo2->chromosome.size())/numchunks;
  auto numchunks1 = numnucs1/chunksize;
  auto numchunks2 = numnucs2/chunksize;
  // 1 goes on the x-axis
  vector<set<NucleotideStore>> xvector, yvector;
  xvector.resize(numnucs1/chunksize);
  yvector.resize(numnucs2/chunksize);

  cout<<numnucs1<<" " <<numnucs2 <<endl;
  cout<<xvector.size()<<" "<<yvector.size()<<endl;
   
  atomic<uint32_t> sofar1{0}, sofar2{0};
  auto f1 = [&sofar1, &rg1, &chromo1, &numchunks1, &chunksize, &xvector]() {
    for(uint32_t chunk = sofar1++ ; chunk < numchunks1; chunk = sofar1++) {
      xvector[chunk]=getUniNucs(*chromo1, chunk*chunksize, chunksize);
      cout<<"1: "<<chunk<<" / " <<numchunks1<<" -> " <<xvector[chunk].size() << endl;
    }
  };
  auto f2 = [&sofar2, &rg2, &chromo2, &numchunks2, &chunksize, &yvector]() {
    for(uint32_t chunk = sofar2++ ; chunk < numchunks2; chunk = sofar2++) {
      yvector[chunk]=getUniNucs(*chromo2, chunk*chunksize, chunksize);
      cout<<"2: "<<chunk<<" / " <<numchunks2<<" -> " <<yvector[chunk].size() << endl;
    }
  };

  vector<std::thread> running;
  for(int n=0; n < 4; ++n) {
    running.emplace_back(f1);
    running.emplace_back(f2);
  }
  
  for(auto& r : running) {
    r.join();
    cout<<"Thread done"<<endl;
  }
  
  running.clear();
  sofar1=0;

  
  auto f3=[&sofar1, &numchunks1, &numchunks2, &chunksize, &xvector, &yvector, &m]() {
    
    
    for(uint32_t xchunk = sofar1++ ; xchunk < numchunks1; xchunk = sofar1++) {
      cout<<xchunk<<endl;
      for(uint32_t ychunk = 0; ychunk < numchunks2; ++ychunk) {
        auto c = set_intersection_size(xvector[xchunk].begin(), xvector[xchunk].end(), yvector[ychunk].begin(), yvector[ychunk].end());
        m(xchunk, ychunk)={xvector[xchunk].size(), yvector[ychunk].size(), c};
      }
    }
  };
  for(int n=0; n < 8; ++n)
    running.emplace_back(f3);
  
  while(sofar1 < numchunks1) {
    cout<<sofar1<<endl;
    sleep(30);
    ofstream plot("plot");
    for(size_t x=0; x <= m.maxX(); ++x)
      for(size_t y=0; y <= m.maxY(); ++y) {
        auto res=m(x,y);
        plot<<x<<"\t"<<y<<"\t"<<res.xnucs<<"\t" << res.ynucs<<"\t"<<res.overnucs<<"\n";
      }
  }
  for(auto& r : running) {
    r.join();
    cout<<"Thread done"<<endl;
  }

  ofstream plot("plot");
  for(size_t x=0; x <= m.maxX(); ++x) {
    for(size_t y=0; y <= m.maxY(); ++y) {
      auto res=m(x,y);
      plot<<x<<"\t"<<y<<"\t"<<res.xnucs<<"\t" << res.ynucs<<"\t"<<res.overnucs<<"\n";
    }
  }
}


uint32_t measureBZ2(const std::string& cmp)
{
  unsigned int dstlen=cmp.size()*1.04;
  char* dest=new char[dstlen];
  int err;
  if((err=BZ2_bzBuffToBuffCompress(dest, &dstlen, (char*)cmp.c_str(), cmp.size(), 9, 0, 30)) != BZ_OK) {
    delete[] dest;
    throw runtime_error("bz2 error: "+std::to_string(err));
  }

  delete[] dest;
  return dstlen;
}

void doCompressMap(const ReferenceGenome& rg)
{
  struct Count
  {
    uint32_t xlen, ylen, totlen;
  };
  Matrix<Count, 4000, 4000> m;
  auto numnucs= rg.numNucleotides();
  unsigned int chunksize = rg.numNucleotides()/8000;
  unsigned int numchunks = numnucs/chunksize - 1; // we have some edge cases

  cout<<"Chunksize: "<<chunksize<<" nucleotides"<<endl;
  
  atomic<uint32_t> sofar{0};
  auto f=[&sofar,&rg,&numchunks,&m,&chunksize]() {
    for(uint32_t xchunk = sofar++ ; xchunk < numchunks; xchunk = sofar++) {
      try {
	auto xstretch=rg.getRange(xchunk*chunksize, chunksize);
        string xascii=xstretch.toASCII();
        auto xlen=measureBZ2(xascii);
        for(uint32_t ychunk = 0 ; ychunk < numchunks; ++ychunk) {
          auto ystretch=rg.getRange(ychunk*chunksize, chunksize);

          string yascii=ystretch.toASCII();
          m(xchunk, ychunk)={
            xlen,
            measureBZ2(yascii),
            measureBZ2(xascii+yascii)};
        }
        
      }
      catch(std::exception& e) {
	cerr<<e.what()<<endl;
      }
    }
  };

  vector<std::thread> running;
  for(int n=0; n < 8; ++n)
    running.emplace_back(f);

  while(sofar < numnucs) {
    cout<<sofar<<endl;
    sleep(30);
    {
      ofstream plot("plot.tmp");
      plot <<"# "<<m.maxX()<<" "<< m.maxY()<<endl;
      for(size_t x=0; x <= m.maxX(); ++x) {
        for(size_t y=0; y <= m.maxY(); ++y) {
          auto res=m(x,y);
          plot<<x<<"\t"<<y<<"\t"<<res.xlen<<"\t" << res.ylen<<"\t"<<res.totlen<<"\t";
          if(res.xlen + res.ylen)
            plot<<1.0*res.totlen/(res.xlen+res.ylen)<<"\n";
          else
            plot<<"0\n";
        }
      }
    }
    rename("plot.tmp", "plot");
  }
  
  for(auto& r : running)
    r.join();

}


int main(int argc, char**argv)
{
  for(unsigned int n=0; n < g_unitsize;++n) {
    g_allA.append('A');
    g_allC.append('C');
    g_allG.append('G');
    g_allT.append('T');
  }
  
  cout<<"Start reading genome"<<endl;
    
  if(argc < 2) {
    cerr<<"Syntax: genex reference.fasta"<<endl;
    return EXIT_FAILURE;
  }
  
  ReferenceGenome rg1(argv[1]);
  string chromo1(argv[2]);
  cout<<"Done reading genome1, have "<<rg1.numChromosomes()<<" chromosomes, "<<
    rg1.numNucleotides()<<" nucleotides"<<endl;

  
  ReferenceGenome rg2(argv[3]);
  string chromo2(argv[4]);

  cout<<"Done reading genome 2, have "<<rg2.numChromosomes()<<" chromosomes, "<<
    rg2.numNucleotides()<<" nucleotides"<<endl;

  
  /*
  for(const auto& a: rg.getAllChromosomes()) {
    ofstream of(a.first);
    cout<<a.first<<" " <<a.second.fullname<< " "<<a.second.chromosome.getString().size()<<endl;
    of<<a.second.chromosome.getString();
  }
  return 0;
  */
  

  // doCompressMap(rg);
  doKMERMap(rg1, chromo1, rg2, chromo2);
  
}

