#include "csv-parser/csv.hpp"
#include <iostream>
#include <map>
#include <fstream>
#include <optional>
#include <fmt/core.h>
#include <fmt/os.h>
#include "refgenome2.hh"
#include "geneannotated.hh"

#include <thread>
#include <atomic>

using namespace std;

/* idea:
   Go through a sequence, and count homopolymers.
   Per sequence emit number of 1-5, 5-10, 10-15, 15-20 etc per ntd
*/

struct ChromoStats
{
  unsigned int acount{0}, ccount{0}, gcount{0}, tcount{0};
  unsigned int size;
  double afrac, cfrac, gfrac, tfrac;
};

ChromoStats getChromoStats(const ReferenceGenome::Chromosome& ch)
{
  const auto& chr = ch.chromosome;
  struct ChromoStats ret={};
  ret.size = chr.size();
  for(uint64_t s = 0 ; s < ret.size; ++s) {
    char c = chr.get(s);
    if(c=='A')
      ++ret.acount;
    else if(c=='C')
      ++ret.ccount;
    else if(c=='G')
      ++ret.gcount;
    else if(c=='T')
      ++ret.tcount;
  }
  ret.afrac = 1.0*ret.acount/ret.size;
  ret.cfrac = 1.0*ret.ccount/ret.size;
  ret.gfrac = 1.0*ret.gcount/ret.size;
  ret.tfrac = 1.0*ret.tcount/ret.size;
  return ret;
}

int main(int argc, char** argv)
{
  std::atomic<int> ctr=1;
  std::mutex iolock;
  ofstream hpcsv("hpcount.csv");
  hpcsv<<"name";
  for(const auto& c : "acgt") {
    if(!c)continue;
    for(int n=1;n<=24;++n) 
      hpcsv<<';'<<c<<"hp"<<n;
    hpcsv<<';'<<c<<"hp25+";
    hpcsv<<';'<<c<<"hpmax";
    hpcsv<<';'<<c<<"hpmean";
    
  }
  hpcsv<<endl;
  
  auto doFunc=[&ctr, &iolock, &argc, &argv, &hpcsv]() {
    for(int n=ctr++; n < argc; n=ctr++) {
      try {
        string fname=argv[n];
        cout<<"Reading "<<fname<<".. ";
        cout.flush();
        ReferenceGenome rg(fname);
        cout<<"Done reading genome from "<<argv[n]<<", have "<<rg.numChromosomes()<<" chromosomes, "<<
          rg.numNucleotides()<<" nucleotides"<<endl;
        
        for(const auto& c : rg.getAllChromosomes()) {
          auto cstats = getChromoStats(c.second);
          
          const auto& chr = c.second.chromosome;
          if(chr.size() < 150000) {
            cout<<c.first<<" too small, "<<chr.size()<<" ntds, skipping: "<<c.second.fullname<<endl;
            continue;
          }
          cout << c.first<<": "<<chr.size() << " ntds ";
          cout<<c.second.fullname<<endl;

          // we skip the first one
          int hpcount=1;
          map<unsigned int, unsigned int> acounts, ccounts, gcounts, tcounts;
          for(uint64_t s = 1 ; s < chr.size(); ++s) {
            char prevchar = chr.get(s-1);
            if(chr.get(s) == prevchar)
              ++hpcount;
            else {
              if(prevchar=='A')
                acounts[hpcount]++;
              else if(prevchar=='C')
                ccounts[hpcount]++;
              else if(prevchar=='G')
                gcounts[hpcount]++;
              else if(prevchar=='T')
                tcounts[hpcount]++;
              hpcount=1;
            }
          }

          std::lock_guard<std::mutex> m(iolock);
          hpcsv<<c.first;
          auto emit=[&](char c, auto& m) {
                      cout<<c<<" homopolymers: ";
                      double frac;
                      if(c=='A') frac = cstats.afrac;
                      else if(c=='C') frac = cstats.cfrac;
                      else if(c=='G') frac = cstats.gfrac;
                      else if(c=='T') frac = cstats.tfrac;
                      cout<<"("<<frac<<")\n";
                      int ex=0;
                      double sum=0;
                      unsigned int count=0;
                      for(const auto& as : m) {
                        count += as.second;
                        sum += as.first * as.second;
                        cout<<as.first<<": "<<as.second<<" (";

                        double est = pow(frac, as.first);
                        for(int k = as.first +1 ; k < 20; ++k)
                          est -= pow(frac, k);
                        
                        cout <<chr.size() * est<<")\n";
                        if(as.first > 24)
                          ex += as.second;
                      }
                      for(int n=1;n<=24;++n) {
                        if(m.find(n) != m.end())
                          hpcsv<<';'<<m[n];
                        else hpcsv<<";0";
                      }
                      
                      hpcsv<<';'<<ex;
                      hpcsv<<';'<< m.rbegin()->first; // highest homopolymer count
                      hpcsv<<';'<< ( sum/count); // weighted average length
                    };
          emit('A', acounts);
          emit('C', ccounts);
          emit('G', gcounts);
          emit('T', tcounts);
          hpcsv<<'\n';
        }
      }
      catch(...){}
    }
    
              };

    thread t1(doFunc);
  thread t2(doFunc);
  thread t3(doFunc);
  thread t4(doFunc);
  t1.join();
  t2.join();
  t3.join();
  t4.join();


}
