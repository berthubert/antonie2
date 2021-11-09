#include "refgenome2.hh"
#include "geneannotated.hh"
#include <errno.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include "dnamisc.hh"
#include <map>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include "misc.hh"
#include <boost/dynamic_bitset.hpp>
#include <unordered_map>

using namespace std;

struct Node
{
  std::string parent;
  vector<std::string> children;
  std::string type;
  std::string tag;
};

std::unordered_map<string, Node> nodes;

int main(int argc, char **argv)
{
  if(argc < 3) {
    cerr<<"Syntax: gfftool annotations.gff refgenome.fna offset1 [offset2]"<<endl;
    return EXIT_FAILURE;
  }

  GeneAnnotationReader gar(argv[1]);
  cout<<"Done with annotations, got "<<gar.size()<<" of them"<<endl;
  auto names=gar.getChromosomes();
  cout<<"Got annotations for "<<names.size()<<" DNA sequences"<<endl;
  cout<<"Reading reference genome"<<endl;
  ReferenceGenome rg(argv[2]);
  
  cout<<"Done reading genome from "<<argv[2]<<", have "<<rg.numChromosomes()<<" sequences, "<<
    rg.numNucleotides()<<" nucleotides"<<endl;

  
  for(const auto&n : names) {
    auto anns = gar.getAll(n);
    cout << n<<"\t"<<anns.size()<<" annotations,";

    auto ptr = rg.getChromosome(n);
    if(ptr)
      cout<<" "<<ptr->fullname<<" have "<<ptr->chromosome.size()<<" nucleotides,";
    map<string, int> tcounts;
    unordered_set<string> genes;
    for(const auto& a : anns) {
      if(a.type=="gene") {
        genes.insert(a.name);
        cout<<' ';
        if(a.strand)
          cout<<'+';
        cout<<a.name<<"("<<a.gene_biotype<<")";
        tcounts[a.gene_biotype]++;
      }
    }
    cout<<": "<<genes.size()<<" genes.";
    for(const auto& t: tcounts)
      cout<<" "<<t.first<<"("<<t.second<<")";
    cout<<endl;

  }

}
