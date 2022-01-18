#include "csv-parser/csv.hpp"
#include "refgenome2.hh"
#include "geneannotated.hh"
#include <errno.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include "dnamisc.hh"
#include <map>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include "misc.hh"
#include <boost/dynamic_bitset.hpp>
#include <unordered_map>

using namespace std;


int main(int argc, char **argv)
{
  if(argc != 2) {
    cerr<<"Syntax: cluster intronsfile"<<endl;
    return 0;
  }
  csv::CSVReader reader(argv[1]);

  auto safeindex = [&reader](const char* name) {
		     int ret = reader.index_of(name);
		     if(ret < 0)
		       throw runtime_error("Could not find index for "+string(name));
		     return ret;
		   };
    
  
  int chrpos = safeindex("chr");
  int startpos = safeindex("start"); //
  int stoppos = safeindex("stop");
  int genepos = safeindex("gene");
  int sensepos = safeindex("sense"); //
  int intronpos = safeindex("intron"); // 


  struct Intron
  {
    string chr;
    string gene;
    int startpos;
    int stoppos;
    bool sense;
    NucleotideStore dna;
  };
  vector<Intron> introns;
  for (csv::CSVRow& row: reader) { // Input iterator
    introns.push_back({
        row[chrpos].get<string>(),
        row[genepos].get<string>(),
        row[startpos].get<int>(),
        row[stoppos].get<int>(),
        row[sensepos].get<bool>(),
        NucleotideStore(row[intronpos].get<string>())});
  }
  cout<<"Have "<<introns.size()<<" introns"<<endl;

  sort(introns.begin(), introns.end(), [](const auto& a, const auto& b) {
    return a.dna.size() < b.dna.size();
  });


  auto thr=[&introns]() {
    for(unsigned int n = 0;; ++n) {
      int rnd = random() % 70000;
      int rnd2 = rnd + 1800 - (random() % 1800);
      const auto& a = introns[rnd];
      const auto& b = introns[rnd2];
      
      auto delta = a.dna.getDelta(b.dna);
      if(delta.size() < 0.25 * a.dna.size() && delta.size() > 0.15 * a.dna.size()) {
        cout<<"A: "<<a.chr<<", "<<a.gene<<", start: "<<a.startpos<<", len: "<<a.dna.size()<<endl;
        cout<<"B: "<<b.chr<<", "<<b.gene<<", start: "<<b.startpos<<", len: "<<b.dna.size()<<", delta: "<<delta.size()<<endl;
        cout<<"link: "<<a.gene<<"-"<<b.gene<<endl;
      }
    }
  };

  srandom(time(0));
  std::thread t1(thr), t2(thr), t3(thr), t4(thr);
  t1.join();
  t2.join();
  t3.join();
  t4.join();

  #if 0
  std::unordered_map<string, set<decltype(introns)::const_iterator>> iIndex;
  unsigned int slen=20;
  unsigned int ctr=0;
  for(auto iter = introns.cbegin() ; iter != introns.cend(); ++iter) {
    if(iter->size() < slen)
      continue;
    for(unsigned int p = 0 ; p < iter->size() - slen; ++p) {
      auto r = iter->getRange(p, slen).toASCII();
      iIndex[r].insert(iter);
    }
    if(!((++ctr)%256)) {
      cout<<"\r"<<ctr<<", "<<iIndex.size()<<" different "<<slen<<"-mers so far";
      cout.flush();
    }
  }
  cout<<"\n";
  cout<<"There are "<<iIndex.size()<<" different "<<slen<<"-mers"<<endl;
  vector<pair<int, NucleotideStore>> rindex;
  for(const auto& i : iIndex) {
    rindex.emplace_back(i.second.size(), i.first);
  }
  sort(rindex.begin(), rindex.end());
  for(auto iter = rindex.rbegin(); iter != rindex.rend(); ++iter) {
    cout<<iter->first<<' '<<iter->second<<'\n';
  }
  #endif
}
