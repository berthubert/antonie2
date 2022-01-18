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

/* reads a whole bunch of genomes, and it construes a filename with gff3 for them too
 * based on this it emits a big cpgstats.csv file with statistics on frequency of CpG
 * 
 * The GFF is used to determine if we are in a gene or not
 *
 */

int main(int argc, char**argv)
{
  if(argc < 2) {
    cerr<<"Syntax: cpgstats reference.fasta"<<endl;
    return EXIT_FAILURE;
  }

  ofstream cpgstats("cpgstats.csv");
  cpgstats<<"name,pos,cgcount,pgenecount,ngenecount";
  cpgstats<<endl;

  
  for(int n=1; n < argc; ++n) {
    try {
      string garname=argv[n+1];
      GeneAnnotationReader gar(garname);
      cout<<"Done with annotations from "<<garname<<" (taxid "<<gar.d_taxonID<<"), got "<<gar.size()<<" of them, for "<<gar.getChromosomes().size()<<" chromosomes: ";

      
      ReferenceGenome rg(argv[n]);

      cout<<"Done reading genome from "<<argv[n]<<", have "<<rg.numChromosomes()<<" chromosomes, "<<
	rg.numNucleotides()<<" nucleotides"<<endl;
      
    
      for(const auto& c : rg.getAllChromosomes()) {
	const auto& chr = c.second.chromosome;

	if(chr.size() < 150000) {
	  cout<<c.first<<" too small, "<<chr.size()<<" ntds, skipping"<<endl;
	  continue;
	}
        auto annos = gar.getAll(c.first);
                
        cout<<c.first<<endl;
        map<int, int> gpos, gneg;

        vector<int> sgpos, sgneg;
        for(const auto& a: annos) {
          if(a.type=="gene" && a.strand == 1)
            sgpos.push_back(a.startPos);
          if(a.type=="gene" && a.strand == 0)
            sgneg.push_back(a.stopPos);

        }
        sort(sgpos.begin(), sgpos.end());
        sort(sgneg.begin(), sgneg.end());
        int gcount=0;
        for(const auto& sgp : sgpos)
          gpos[sgp]=gcount++;
        gcount=0;
        for(const auto& sgn : sgneg)
          gneg[sgn]=gcount++;
        
        int cgcount=0, gccount=0;
        char prev = chr.get(0);
        
	for(uint64_t s = 1 ; s < chr.size(); ++s) {
          
          char cur = chr.get(s);
          if(prev=='C' && cur=='G')
            cgcount++;

          prev = cur;
          
          if(!(s%256)) {
            int gposcount=gpos.rbegin()->second, gnegcount=gneg.rbegin()->second;
            if(auto iter = gpos.lower_bound(s); iter != gpos.end())
              gposcount=iter->second;
            if(auto iter = gneg.lower_bound(s); iter != gneg.end())
              gnegcount=iter->second;            
            cpgstats<<c.first<<','<<s<<','<<cgcount<<','<< gposcount<<','<<gnegcount<<'\n';
          }
	}
      }
    }
    catch(std::exception& e) {
      cerr<<"Exception: "<<e.what()<<endl;
    }
  }
  
}


