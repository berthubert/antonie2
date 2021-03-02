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
 * based on this it emits a big skplot.csv file which has all kinds of skews for all chromosomes, posted every n loci
 * also emits a genomes.csv summary file (should be 'chromosomes' really)
 *
 * The GFF is used to determine the codon position, or if we are in a gene or not
 *
 * The GFF has a link to a species ID and it would be super nice if we could get data from there
 * "species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=869304" - this is html though
 *
 * This code ignores chromosomes smaller than 1 million bp, which mostly rids us of confusing plasmids, viruses etc
 */

int main(int argc, char**argv)
{
  if(argc < 2) {
    cerr<<"Syntax: gcstats reference.fasta"<<endl;
    return EXIT_FAILURE;
  }

  ofstream genomescsv("genomes.csv");
  genomescsv<<"name;fullname;acount;ccount;gcount;tcount;realm1;realm2;realm3;realm4;realm5\n";
  ofstream skplot("skplot.csv");
  skplot<<"name,relpos,abspos,gcskew,taskew,gcskew0,gcskew1,gcskew2,gcskewNG,taskew0,taskew1,taskew2,taskewNG,pospos,gccount";
  for(int cpos = 0 ; cpos < 4 ; ++cpos) {
    skplot << ",acounts"<<cpos<<",ccounts"<<cpos<<",gcounts"<<cpos<<",tcounts"<<cpos;
  }
  skplot<<endl;


  cout<<"Reading taxonomies..";
  cout.flush();
  TaxoReader tr("/home/ahu/git/antonie/taxonomy/new/fullnamelineage.dmp");
  cout<<" got "<<tr.size()<<" entries"<<endl;
  
  for(int n=1; n < argc; ++n) {
    try {
    ReferenceGenome rg(argv[n]);
    string garname = argv[n];
    garname.replace(garname.size()-3, 3, "gff");
  
    
    cout<<"Done reading genome from "<<argv[n]<<", have "<<rg.numChromosomes()<<" chromosomes, "<<
      rg.numNucleotides()<<" nucleotides"<<endl;
    
    GeneAnnotationReader gar(garname);
    cout<<"Done with annotations from "<<garname<<" (taxid "<<gar.d_taxonID<<"), got "<<gar.size()<<" of them, for "<<gar.getChromosomes().size()<<" chromosomes: ";

    for(const auto& c : gar.getChromosomes())
      cout<<"'"<<c<<"' ";
    cout<<endl;
    
    for(const auto& c : rg.getAllChromosomes()) {
      const auto& chr = c.second.chromosome;
      if(chr.size() < 1000000) {
	cout<<c.first<<" too small, "<<chr.size()<<" ntds, skipping"<<endl;
	continue;
      }
      cout << c.first<<": "<<chr.size() << " ntds ";
      cout<<c.second.fullname<<endl;


      
      int aCount{0}, cCount{0}, gCount{0}, tCount{0};
      for(uint64_t s = 0 ; s < chr.size(); ++s) {
	char n = chr.get(s);
	if(n=='A') ++aCount;
	else if(n=='C') ++cCount;
	else if(n=='G') ++gCount;
	else if(n=='T') ++tCount;
      }

      vector<string> taxo = tr.get(gar.d_taxonID);
      taxo.resize(6); // so the printing is safe
      genomescsv<<c.first<<";"<<c.second.fullname<<";"<<aCount<<";"<<cCount<<";"<<gCount<<";"<<tCount<<";"<<taxo[1]<<";"<<taxo[2]<<";"<<taxo[3]<<";"<<taxo[4]<<";"<<taxo[5]<<endl;
      
      cout<<"A: "<<aCount<<", ";
      cout<<"C: "<<cCount<<", ";
      cout<<"G: "<<gCount<<", ";
      cout<<"T: "<<tCount<<", ";
      double gc = 1.0*(cCount + gCount)/(aCount + cCount + gCount +tCount);
      cout<<"gc: "<<gc<<endl;
      double rat = 1.0*(cCount + gCount)/(aCount + tCount);
      double gcskew = 0, taskew=0, skew=0, pospos=0;
      double gcskews[4] = {0,0,0,0};
      double taskews[4] = {0,0,0,0};
      int acounts[4]={}, ccounts[4]={}, gcounts[4]={}, tcounts[4]={};
      
      GeneAnnotation last;
      last.gene=false;
      unsigned int gccount=0;
      for(uint64_t s = 0 ; s < chr.size(); ++s) {
	unsigned genepos = s+1;
	if(!(last.gene && last.startPos <= genepos && genepos < last.stopPos)) {
	  last.gene=false;
	  auto ans= gar.lookup(c.first, genepos);
	  for(const auto& a: ans) {
	    if(a.gene==true && (!a.tag.empty() && a.tag[0]=='g')) {
	      //	      cout<<"At position "<< genepos <<" of "<<c.first<<" we are in gene '"<<a.name<<"' tag '"<<a.tag<<"' which goes from " << a.startPos <<" to "<< a.stopPos <<", strand: "<<a.strand<<". Offset: "<< genepos - a.startPos<<", codon offset: ";
	      
	      if(a.strand) { // should be ATG
		int codonpos=(genepos - a.startPos) % 3;
		//		cout << codonpos << " ";

		/*
		if(!codonpos)
		  cout<<chr.get(s-1)<<" "<<chr.get(s)<<chr.get(s+1)<<chr.get(s+2)<<" "<<chr.get(s+3);
		cout<<"\n";
		*/
	      }
	      else {
		int codonpos = (a.stopPos - genepos) % 3;
		//		cout << codonpos << " ";
		if(!codonpos) {
		  // CAT on the positive strand
		  // cout<<chr.get(s-2)<<chr.get(s-1)<<chr.get(s);
		}
		//		cout<<"\n";
	      }
	      last = a;
	      break;
	    }
	  }
	}
	int codonpos=3; // this means "not in a gene"
	if(last.gene) {
	  if(last.strand) {
	    codonpos = ((genepos - last.startPos) % 3);
	    pospos++;
	  }
	  else {
	    codonpos = ((last.stopPos - genepos) % 3);
	    pospos--;
	  }
	}
	    
	char n = chr.get(s);
	if(n=='G') {
	  gcskew += 1.0;
	  ++gccount;
	}
	else if(n=='C') {
	  gcskew -= 1.0;
	  ++gccount;
	}
	else if(n=='T')
	  taskew += 1.0;
	else if(n=='A')
	  taskew -= 1.0;

	if(n=='G') {
	  gcskews[codonpos] += 1.0;
	  gcounts[codonpos]++;
	} else if(n=='C') {
	  gcskews[codonpos] -= 1.0;
	  ccounts[codonpos]++;
	}
	else if(n=='T') {
	  taskews[codonpos] += 1.0;
	  tcounts[codonpos]++;
	}
	else if(n=='A') {
	  taskews[codonpos] -= 1.0;
	  acounts[codonpos]++;
	}
	
	if(!(s%4096) || s == chr.size() -1 ) {
	  skplot<<c.first<<","<<1.0*s/chr.size()<<","<<s<<","<<gcskew<<","<<taskew<<","<<gcskews[0]<<","<<gcskews[1]<<","<<gcskews[2]<<"," <<gcskews[3]<<","<<taskews[0]<<","<<taskews[1]<<","<<taskews[2]<<"," <<taskews[3]<<","<<pospos<<","<<gccount;
	  for(int cpos = 0 ; cpos < 4 ; ++cpos) {
	    skplot << "," << acounts[cpos]<<","<<ccounts[cpos]<<","<<gcounts[cpos]<<","<<tcounts[cpos];
	  }
	  skplot<<endl;
	}
      }
      cout<<"skew: "<<skew<<", rat: "<<rat<<endl;

    
    }
    }
    catch(std::exception& e) {
      cerr<<"Exception: "<<e.what()<<endl;
    }
  }
}

