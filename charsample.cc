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
#include "taxoreader.hh"

using namespace std;

/* reads a genome and samples it for c/g deltas
 * This code ignores chromosomes smaller than 1 million bp, which mostly rids us of confusing plasmids, viruses etc
 */

int main(int argc, char**argv)
{
  if(argc < 2) {
    cerr<<"Syntax: charsample reference.fasta"<<endl;
    return EXIT_FAILURE;
  }

  ofstream cgcsv("cgsample.csv");
  cgcsv<<"chromo,startpos,siz,gcount1,ccount1,gcount2,ccount2,delta"<<endl;
  for(int n=1; n < argc; ++n) {
    try {
      ReferenceGenome rg(argv[n]);
      //    string garname = argv[n];
      //    garname.replace(garname.size()-3, 3, "gff");
  
    
      cout<<"Done reading genome from "<<argv[n]<<", have "<<rg.numChromosomes()<<" chromosomes, "<<
	rg.numNucleotides()<<" nucleotides"<<endl;
    
      //    GeneAnnotationReader gar(garname);
      //    cout<<"Done with annotations from "<<garname<<" (taxid "<<gar.d_taxonID<<"), got "<<gar.size()<<" of them, for "<<gar.getChromosomes().size()<<" chromosomes: ";

    
      for(const auto& c : rg.getAllChromosomes()) {
	const auto& chr = c.second.chromosome;
	if(chr.size() < 1000000) {
	  cout<<c.first<<" too small, "<<chr.size()<<" ntds, skipping"<<endl;
	  continue;
	}
	for(int siz = 25 ; siz < 5000; siz += 25) {
	  int iters=2000;
	  if(siz <= 500)
	    iters=20000;
	  for(int n=0; n < iters; ++n) {
	    int startpos= random() % (chr.size() - 2*siz);
	    int cCount1{0}, gCount1{0};
	    int cCount2{0}, gCount2{0};
	    int pos = startpos;
	    for(; pos < startpos + siz; ++pos) {
	      char c = chr.get(pos);
	      if(c=='G')
		++gCount1;
	      else if(c=='C')
		++cCount1;
	    }
	    startpos += siz;
	    for(; pos < startpos + siz; ++pos) {
	      char c = chr.get(pos);
	      if(c=='G')
		++gCount2;
	      else if(c=='C')
		++cCount2;
	    }

	    
	    cgcsv << c.first<<","<<startpos<<","<<siz<<","<<gCount1<<","<<cCount1<<","<<gCount2<<","<<cCount2<<",";
	    cgcsv << (gCount1-cCount1) << "\n";
	  }
	}
      }
    }
    catch(std::exception& e) {
      cerr<<"Exception: "<<e.what()<<endl;
    }
  }
  
}


