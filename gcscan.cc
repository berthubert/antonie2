#include "refgenome2.hh"
#include "geneannotated.hh"

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

// reads one chromosome and bunches statistics about them together in two csv files
// scan.csv:
// codons.csv: 

int main(int argc, char**argv)
{
  if(argc < 2) {
    cerr<<"Syntax: gcscan reference.fasta"<<endl;
    return EXIT_FAILURE;
  }

  ReferenceGenome rg(argv[1]);
  
  cout<<"Done reading genome from "<<argv[1]<<", have "<<rg.numChromosomes()<<" chromosomes, "<<
    rg.numNucleotides()<<" nucleotides"<<endl;
  
  for(const auto& c : rg.getAllChromosomes()) {
    const auto& chr = c.second.chromosome;
    cout << c.first<<": "<<chr.size() << " ntds ";
    cout<<c.second.fullname<<endl;

    int aCount{0}, cCount{0}, gCount{0}, tCount{0};
    ofstream scan("scan.csv");
    ofstream codons("codons.csv");
    
    scan<<"name,relpos,abspos,atCount,gcCount,c\n";
    codons<<"name,abspos,codon"<<endl;
    for(uint64_t s = 0 ; s < chr.size(); ++s) {
      char n = chr.get(s);
      if(n=='A') ++aCount;
      else if(n=='C') ++cCount;
      else if(n=='G') ++gCount;
      else if(n=='T') ++tCount;
      scan<<c.first<<","<<1.0*s/chr.size()<<","<<s<<","<<aCount+tCount<<","<<gCount+cCount<<","<<n<<"\n";
      if(!(s%3))
	codons << c.first<<","<<s<<","<<chr.get(s)<<chr.get(s+1)<<chr.get(s+2)<<endl;
    }
    scan<<endl;
  }
}

