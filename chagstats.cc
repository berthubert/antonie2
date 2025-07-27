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

/* reads a whole bunch of genomes, and it construes a filename with gff3 for them too
 * based on this it emits a big chagraff.csv file with statistics for nucleotide
 * triplets
 * 
 * The GFF is used to determine the codon position, or if we are in a gene or not
 *
 * This code ignores chromosomes smaller than 1 million bp, which mostly rids us of confusing plasmids, viruses etc
 */

int main(int argc, char**argv)
{
  if(argc < 2) {
    cerr<<"Syntax: chagstats reference.fasta"<<endl;
    return EXIT_FAILURE;
  }

  ofstream chagstats("chagstats.csv");
  chagstats<<"name";
  char nucs[]="ACGT";
  char triplet[4];
  triplet[3]=0;
  for(int i=0; i < 4; ++i) {
    triplet[0]=nucs[i];
    for(int j=0; j < 4; ++j) {
      triplet[1]=nucs[j];
      for(int k=0; k < 4; ++k) {
	triplet[2]=nucs[k];
	chagstats<<","<<triplet;
      }
    }
  }
  chagstats<<endl;

  const int ngramlen=7;
  char ngram[ngramlen+1];
  ngram[ngramlen]=0;
  
  map<string, int64_t> ngramcount;
  visitAllNgrams([&ngramcount](const auto& str) { ngramcount[str]=0; }, ngramlen );
  
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
	if(chr.size() < 150000) {
	  cout<<c.first<<" too small, "<<chr.size()<<" ntds, skipping"<<endl;
	  continue;
	}
	triplet[1]=chr.get(1);
	triplet[2]=chr.get(2);
	map<string, unsigned int> tripcount;
	for(uint64_t s = 2 ; s < chr.size(); ++s) {
	  triplet[0]=triplet[1];
	  triplet[1]=triplet[2];
	  triplet[2] = chr.get(s);
	  
	  tripcount[triplet]++;
	}
	chagstats<<c.first;
	for(const auto& p : tripcount)
	  chagstats<<','<<p.second;
	chagstats<<"\n";


	for(int n=1; n < ngramlen  ; ++n)
	  ngram[n]=chr.get(n);

	for(uint64_t s = ngramlen-1 ; s < chr.size(); ++s) {
	  for(int n=0; n < ngramlen  -1; ++n)
	    ngram[n]=ngram[n+1];

	  ngram[ngramlen-1] = chr.get(s);
	  
	  ngramcount[ngram]++;
	}
      }
    }
    catch(std::exception& e) {
      cerr<<"Exception: "<<e.what()<<endl;
    }
  }
  ofstream ngcsv("ngrams.csv");
  ngcsv<<"ngram,rcngram,count,rccount,procdiff"<<endl;
  for(const auto& nc : ngramcount) {
    NucleotideStore ns(nc.first);
    if(!ns.isCanonical())
      continue;
    string rev = ns.getRC().toASCII();
    ngcsv<<nc.first<<","<<rev<<","<<nc.second<<','<<ngramcount[rev]<<',';
    ngcsv<< (100.0*(1.0*ngramcount[rev]-1.0*nc.second)/nc.second) << '\n';
  }
  
}


