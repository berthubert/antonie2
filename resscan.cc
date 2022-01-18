#include "refgenome2.hh"
#include "geneannotated.hh"

#include <iostream>
#include "misc.hh"
#include <set>
#include <algorithm>
#include "dnamisc.hh"
#include <boost/algorithm/string.hpp>
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

// reads a fasta file & associated GFF and counts k-mers for all genes

static int countStrings(const std::string& s, const std::string& tocount)
{
  int count = 0;
  size_t nPos = s.find(tocount, 0); // first occurrence
  while(nPos != string::npos) {
    count++;
    nPos = s.find(tocount, nPos + 1);
  }
  return count;
}

int main(int argc, char**argv)
{
  if(argc < 2) {
    cerr<<"Syntax: resscan reference.fasta"<<endl;
    return EXIT_FAILURE;
  }
  ofstream rescsv("rescount.csv"), sizescsv("ressizes.csv");
  rescsv<<"chromo,sequence,genecount,precount,totcount\n";
  sizescsv<<"chromo,genesize,presize,totsize\n";
  
  std::atomic<int> ctr=1;
  std::mutex iolock;

  int a0, a1, a2, a3, a4, a5;
  char acgt[]="ACGT";
  vector<string> all;
  string str;
  str.resize(6);
   
  for(a0=0; a0 < 4 ; a0++) {
    str[0]=acgt[a0];
    for(a1=0; a1 < 4 ; a1++) {
      str[1]=acgt[a1];
      for(a2=0; a2 < 4 ; a2++) {
        str[2]=acgt[a2];
        for(a3=0; a3 < 4 ; a3++) {
          str[3]=acgt[a3];
          for(a4=0; a4 < 4 ; a4++) {
            str[4]=acgt[a4];
            for(a5=0; a5 < 4 ; a5++) {
              str[5]=acgt[a5];
              all.push_back(str);
            }
          }
        }
      }
    }
  }
  
  auto doFunc=[&ctr, &iolock, &argc, &argv, &rescsv, &sizescsv, &all]() {
    for(int n=ctr++; n < argc; n=ctr++) {
      ReferenceGenome rg(argv[n]);
      string garname = argv[n];
      boost::replace_all(garname, ".fna", ".gff");
      
      GeneAnnotationReader gar(garname);
      cout<<"Done with annotations from "<<garname<<" (taxid "<<gar.d_taxonID<<"), got "<<gar.size()<<" of them, for "<<gar.getChromosomes().size()<<" chromosomes\n";
      
      cout<<"Done reading genome from "<<argv[n]<<", have "<<rg.numChromosomes()<<" chromosomes, "<<
        rg.numNucleotides()<<" nucleotides"<<endl;
      for(const auto& c : rg.getAllChromosomes()) {
        const auto& chr = c.second.chromosome;
        if(chr.size() < 150000) {
          cout<<c.first<<" too small, "<<chr.size()<<" ntds, skipping: "<<c.second.fullname<<endl;
          continue;
        }

        
        cout << c.first<<": "<<chr.size() << " ntds ";
        cout<<c.second.fullname<<endl;
        
        auto anno = gar.getAll(c.first);
        map<string, int> genecounts, precounts, totcounts;
        unsigned int genesize{0}, presize{0}, totsize{0};
        for(const auto& a : all) {
          genecounts[a];
          precounts[a];
        }
        for(const auto& a: anno) {
          if(a.type=="gene") {
            //                    cout<<a.name<<": "<<a.stopPos - a.startPos<<" "<<a.strand<<": ";
            NucleotideStore ntds, pre;
            if(a.strand) {
              ntds=chr.getRange(a.startPos - 1, a.stopPos - a.startPos);
              if(a.startPos > 100)
                pre = chr.getRange(a.startPos - 100, 100);
            }
            else {
              ntds=chr.getRange(a.startPos, a.stopPos - a.startPos).getRC();
              if(a.stopPos + 100 < chr.size())
                pre = chr.getRange(a.stopPos, 100).getRC();
            }
            
            auto a= ntds.toASCII();
            genesize += a.size();
            for(string::size_type pos = 6; pos < a.size(); ++pos) {
              genecounts[a.substr(pos-6, 6)]++;
            }
            
            a=pre.toASCII();
            presize += a.size();
            for(string::size_type pos = 6; pos < a.size(); ++pos) {
              precounts[a.substr(pos-6, 6)]++;
            }
            
          }
        }

        totsize = chr.size();
        for(string::size_type pos = 6; pos < chr.size(); ++pos) {
          NucleotideStore ns =chr.getRange(pos-6, 6);
          if(!ns.isCanonical())
            ns = ns.getRC();
          totcounts[ns.toASCII()]++;
        }

        ostringstream tmp, tmpsizes;
        tmpsizes<<c.first<<","<<genesize<<","<<presize<<","<<totsize; 
        for(auto& gc : genecounts) {
          //          cout<<gc.first<<": Totcount: "<<gc.second<<", precount: "<<precounts[gc.first]<<endl;
          NucleotideStore ns(gc.first);
          if(!ns.isCanonical())
            ns = ns.getRC();
          tmp<<c.first<<","<<gc.first<<","<<gc.second<<","<<precounts[gc.first]<<","<<totcounts[ns.toASCII()]<<"\n";

        }
        std::lock_guard<std::mutex> m(iolock);
        rescsv << tmp.str();
        sizescsv << tmpsizes.str() <<endl;
                
      }
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

