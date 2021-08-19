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
#include "csv-parser/csv.hpp"

using namespace std;


int main(int argc, char**argv)
{
  if(argc < 2) {
    cerr<<"Syntax: chromopic reference.fasta"<<endl;
    return EXIT_FAILURE;
  }
  for(int n=1; n < argc; ++n) {
    try {

      ReferenceGenome rg(argv[n]);
      string garname = argv[n];
      garname.replace(garname.size()-3, 3, "gff");
  
    
      cout<<"Done reading genome from "<<argv[n]<<", have "<<rg.numChromosomes()<<" chromosomes, "<<
	rg.numNucleotides()<<" nucleotides"<<endl;
    
      GeneAnnotationReader gar(garname);
      cout<<"Done with annotations from "<<garname<<" (taxid "<<gar.d_taxonID<<"), got "<<gar.size()<<" of them, for "<<gar.getChromosomes().size()<<" chromosomes: ";

    
      for(const auto& c : rg.getAllChromosomes()) {
	cout<<c.first<<endl;
	const auto& chr = c.second.chromosome;

	ofstream piccsv(c.first+"_pic.csv");
	piccsv<<"pos,sense\n";
	if(chr.size() < 1000000) {
	  cout<<c.first<<" too small, "<<chr.size()<<" ntds, skipping"<<endl;
	  continue;
	}

	auto annotations = gar.getAll(c.first);
	cout<<"Got "<<annotations.size()<<" annotations"<<endl;
	for(const auto& a : annotations) {
	  if(a.type=="gene" || a.type=="pseudogene") {
	    if(a.strand) {
	      for(auto pos = a.startPos; pos < a.stopPos; pos++) {
		if(!(pos % 16))
		  piccsv<<pos<<","<<"+"<<"\n";
	      }

	    }
	    else {
 	      for(auto pos = a.stopPos; pos > a.startPos; pos--) {
		if(!(pos % 16))
		  piccsv<<pos<<","<<"-"<<"\n";

	      }
	    }
	  }
	}
      }
    }
    catch(std::exception& e) {
      cerr<<"Error in chromosome "<<e.what()<<endl;
    }
  }
}
