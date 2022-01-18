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
#include <boost/algorithm/string.hpp>
#include <future>
#include <sstream>
#include <sys/prctl.h>
#include "bzlib.h"
#include "taxoreader.hh"

using namespace std;

/* grep a fasta file, per chromosome
 *
 */

int main(int argc, char**argv)
{
  if(argc < 2) {
    cerr<<"Syntax: fagrep nucleotides fasta [nucleotides]"<<endl;
    return EXIT_FAILURE;
  }
  NucleotideStore srch(argv[1]);

  NucleotideStore srch2;
  if(argc > 3)
    srch2 = NucleotideStore(argv[3]);
  
  string fnaname(argv[2]);
  ReferenceGenome rg(fnaname);

  // this is a pretty weak effort, full of bugs
  
  cout<<"Done reading reference genome from '"<<fnaname<<"', got chromosomes: "<<endl;
  for(const auto& c : rg.getAllChromosomes()) {
    const auto& chromo = c.second;    
    auto ascii=chromo.chromosome.toASCII();
    
    auto iter=ascii.find(srch.toASCII());
    if(iter != string::npos)
      cout<<c.first<<" "<<iter<<": Found at least one copy in forward direction"<<endl;
    iter=ascii.find(srch.getRC().toASCII());
    if(iter != string::npos)
      cout<<c.first<<" "<<iter<<": Found at least one copy in reverse direction"<<endl;

    if(!srch2.size())
      continue;

    auto iter2=ascii.find(srch2.toASCII());
    if(iter2 != string::npos)
      cout<<"secondary "<<c.first<<" "<<iter2<<" ("<< ((int)iter2-(int)iter) <<"): Found at least one copy in forward direction"<<endl;
    iter2=ascii.find(srch2.getRC().toASCII());
    if(iter2 != string::npos)
      cout<<"secondary "<<c.first<<" "<<iter2<<" ("<<((int)iter2-(int)iter)<<"): Found at least one copy in reverse direction"<<endl;

    
  }
}


