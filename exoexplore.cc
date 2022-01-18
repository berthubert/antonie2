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
  if(argc < 2) {
    cerr<<"Syntax: exoexplore annotations.gff"<<endl;
    return EXIT_FAILURE;
  }

  GeneAnnotationReader gar(argv[1]);
  cout<<"Done with annotations, got "<<gar.size()<<" of them"<<endl;

  ReferenceGenome rg(argv[2]);
  cout<<"Done reading genome from "<<argv[2]<<", have "<<rg.numChromosomes()<<" chromosomes, "<<
    rg.numNucleotides()<<" nucleotides"<<endl;

  map<string,string> chrotrans;
  for(const auto& c : rg.getAllChromosomes()) {
    // NC_0000001.x
    int num=atoi(c.first.substr(3).c_str());
    chrotrans["chr"+to_string(num)]=c.first;

    cout<<c.first<<": "<<num<<", "<<("chr"+to_string(num))<<endl;
        
  }

  struct GeneData
  {
    string chromo;
    bool sense;
  };
  unordered_map<string, GeneData> chromomap;
  map<string, vector<pair<int,int> > > genex;
  ofstream exoncsv("exons.csv");
  exoncsv<<"chr,start,stop,gene,sense,index"<<endl;
  for(const auto& c : gar.getChromosomes()) {

    auto annos = gar.getAll(c);
    for(const auto& an : annos) {
      if(an.type=="exon") {
        genex[an.enclosing_gene].push_back({an.startPos, an.stopPos});
        chromomap[an.enclosing_gene]={an.chromosome, an.strand};
      }
    }
  }

  ofstream introncsv("introns.csv");
  introncsv<<"chr,start,stop,gene,sense,index"<<endl;
  set<tuple<string,int,int>> seen;
  ofstream allintronscsv("intronsfull.csv");
  allintronscsv<<"chr,start,stop,gene,sense,intron,index"<<endl;
  
  for(auto& g : genex) {
    bool sense = chromomap[g.first].sense;

            
    if(g.second.size() == 1) {
      introncsv<<chromomap[g.first].chromo<<",0,0,"<< sense <<"\n";
      continue; // no introns
    }

    sort(g.second.begin(), g.second.end());
    int count=1;
    exoncsv<<chromomap[g.first].chromo<<","<<g.second.begin()->first<<","<<g.second.begin()->second<<","<<g.first<<","<<sense;
    if(sense)
      exoncsv<<",1\n";
    else
      exoncsv<<","<<g.second.size()<<"\n";
    for(auto iter = g.second.begin() + 1; iter != g.second.end(); ++iter) {
      int relcount = sense ? count : (g.second.size() - count);
      exoncsv<<chromomap[g.first].chromo<<","<<iter->first<<","<<iter->second<<","<<g.first<<","<<sense<<","<<relcount<<"\n";

      
      ++count;
      introncsv<<chromomap[g.first].chromo<<","<<(iter - 1)->second << ","<<iter->first<<","<<g.first<<","<<chromomap[g.first].sense<<","<<relcount<<'\n';
      allintronscsv<<chromomap[g.first].chromo<<","<<(iter - 1)->second << ","<<iter->first<<","<<g.first<<","<<chromomap[g.first].sense<<","<<relcount;

      if(auto iter2 = chrotrans.find(chromomap[g.first].chromo); iter2 != chrotrans.end()) {
        const auto& chr = rg.getChromosome(iter2->second);
        if(!chr)
          continue;
        auto len = iter->first - (iter-1)->second;
        if(len >=0 && len < 100000) {
          tuple<string,int,int> token={iter2->second, (iter-1)->second, len-1};
          if(!seen.count(token)) {
            auto r=chr->chromosome.getRange((iter-1)->second, len-1);
            if(!chromomap[g.first].sense)
              r=r.getRC();
            allintronscsv<<','<<r<<"\n";
            seen.insert(token);
          }
        }
      }
    }
  }
}
