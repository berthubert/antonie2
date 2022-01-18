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
    cerr<<"Syntax: simsearch 1.gff 2.gff"<<endl;
    return EXIT_FAILURE;
  }
#if 0
  ofstream cpgstats("cpgstats.csv");
  cpgstats<<"name,pos,cgcount,pgenecount,ngenecount";
  cpgstats<<endl;
#endif


  struct Gene
  {
    int idx;
    string chromo;
    int startPos;
    int stopPos;
    bool strand;
    NucleotideStore ntds;
    bool operator<(const Gene& b) const
    {
      return idx < b.idx;
    }
  };
  map<string, vector<Gene> > cnts;

  /* 
     Read all chromosomes, find common gene names

     For every common gene name, gather the sequences

     Iterate over each sequence, making 18-mers and note their positions
     
     Then determine all inter-k-mer distances on a sequence, and note those between 100 and 1000
     Do the same for the other sequences
     Determine if there are inter-k-mer distances that exist on all sequence, but differ by at least 100
  */
  map<pair<int, string>, int> exolens;
  for(int n=1; n < argc; ++n) {
    try {
      string garname=argv[n];
      GeneAnnotationReader gar(garname);
      cout<<"Done with annotations from "<<garname<<" (taxid "<<gar.d_taxonID<<"), got "<<gar.size()<<" of them, for "<<gar.getChromosomes().size()<< " sequences\n";

      string fnaname = argv[n];
      boost::replace_all(fnaname, ".gff", ".fna");

      ReferenceGenome rg(fnaname);
      cout<<"Done reading reference genome from '"<<fnaname<<"', got chromosomes: "<<endl;
      for(const auto& c : rg.getAllChromosomes()) {
        cout<<c.first<<endl;
      }
      
      for(const auto& c : gar.getChromosomes()) {
        auto annos = gar.getAll(c);
        const auto& chromo = rg.getChromosome(c);
        if(!chromo) {
          cout<<"Could not get nucleotides for '"<<c<<"'"<<endl;
          continue;
        }
        for(const auto& a: annos) {
          if(a.type=="gene") {
            //            cout<<a.name<<endl;
            auto nc = chromo->chromosome.getRange(a.startPos, a.stopPos - a.startPos);
            if(!a.strand)
              nc = nc.getRC();
            if(nc.size() < 200000)
              cnts[a.name].push_back({n-1, c, (int)a.startPos, (int)a.stopPos, a.strand, nc});
            else {
              NucleotideStore emp;
              cnts[a.name].push_back({n-1, c, (int)a.startPos, (int)a.stopPos, a.strand, emp});
            }
          }
          if(a.type=="exon") {
            exolens[{n-1,(string)a.enclosing_gene}] += a.stopPos - a.startPos;
          }
        }
      }
    }
    catch(std::exception& e) {
      cerr<<"Exception: "<<e.what()<<endl;
    }
  }
  unsigned int nsets=argc-1;
  ofstream gcmp("genecmp.csv");
  gcmp<<"gene;len1;len2;exolen1;exolen2"<<endl;
  ofstream sbs("sbs");
  for(const auto& c : cnts) {

    set<Gene> u;
    for(const auto& i : c.second) {
      u.insert(i);
    }
    if(u.size() ==nsets && u.size() == c.second.size()) {
      cout<<c.first<<":";
      gcmp<<c.first;
      for(const auto& i : u) {
        cout<<" "<<i.idx << " ("<<i.stopPos - i.startPos<<")";
        gcmp<<";"<<(i.stopPos - i.startPos);
      }
      gcmp<<";"<<exolens[{1,c.first}]<<";"<<exolens[{2,c.first}];
      gcmp<<"\n";
      cout<<"\n";
      map<NucleotideStore, set<int>> present;
      int primerlen=20;
      for(const auto& i : u) {

        sbs<<c.first<< " " << i.idx <<" "<<i.ntds.toASCII()<<"\n";
        for(int n=0; n < (int)i.ntds.size() - primerlen; ++n) {
          auto p = i.ntds.getRange(n, primerlen);
          present[p].insert(i.idx);
        }
      }
      set<NucleotideStore> candidates;
      for(const auto& p : present) {
        if(p.second.size() == nsets) {
          //          cout<<p.first<<endl;
          candidates.insert(p.first);
        }
      }
      cout<<"Gene has "<<candidates.size()<<" primer candidates"<<endl;
      


      vector<map<pair<NucleotideStore,NucleotideStore>, int>> dists;
      dists.resize(nsets);
      for(const auto& au : u) {
        auto a = au.ntds.toASCII();
        for(auto iter1 = candidates.cbegin(); iter1 != candidates.cend(); ++iter1) {
          auto c1ascii = iter1->toASCII();
          int pos1 = a.find(c1ascii);
          for(auto iter2 = candidates.cbegin(); iter2 != candidates.cend(); ++iter2) {
            auto c2ascii = iter2->toASCII();
            int pos2 = a.find(c2ascii);
            int dist = pos2-pos1;
            if(dist > 75 && dist < 1300)
              dists[au.idx][{*iter1, *iter2}]=dist;
          }
        }
      }

      map<pair<NucleotideStore, NucleotideStore>, vector<unsigned int>> distances;
      cout<<"Have";
      for(const auto& d : dists) {
        for(const auto& dp : d)
          distances[dp.first].push_back(dp.second);
        cout<<" "<<d.size();
      }
      cout<<" well-spaced pairs, of which ";
      int allcounter=0;


      for(auto& ac : distances) {
        if(ac.second.size() != nsets)
          continue;
        allcounter++;
        auto unsorted=ac.second;
        sort(ac.second.begin(), ac.second.end());
        auto iter = ac.second.begin() + 1;
        for(; iter != ac.second.end(); ++iter) {
          if((*iter - *(iter-1)) < 75)
            break;
        }
        if(iter == ac.second.end()) {
          cout<<"Candidate: "<<ac.first.first <<" - "<<ac.first.second<<":";
          for(const auto& dd : unsorted)
            cout<<" "<<dd;
          cout<<endl;
        }
        
      }
      cout<<allcounter<<" present in all "<<nsets<<" genes"<<endl;

      
      /*
      for(const auto& d : dists1) {
        auto iter = dists2.find(d.first);
        if(iter == dists2.end())
          continue;
        if(abs(iter->second - d.second) > 100) {
          cout<<"Found a sane candidate: "<<d.first.first<<" - " <<d.first.second;
          cout<<": delta "<< (iter->second - d.second) <<", absolutes "<<iter->second<<" - "<<d.second<<endl;
        }

      }
      */
    }
  }
}


