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

/* 
 * Reads many genomes and outputs codon histograms for genes, split out
 * by if they live on the positive or negative strand.
 *
 * Note that this is not the leading or lagging strand!
 *
 * This code ignores chromosomes smaller than 1 million bp, which mostly rids us of confusing plasmids, viruses etc
 */

int main(int argc, char**argv)
{
  if(argc < 2) {
    cerr<<"Syntax: genehisto reference.fasta"<<endl;
    return EXIT_FAILURE;
  }

  csv::CSVReader rescsv("results.csv");
  auto safeindex = [&rescsv](const char* name) {
		     int ret = rescsv.index_of(name);
		     if(ret < 0)
		       throw runtime_error("Could not find index for "+string(name));
		     return ret;
		   };
  int name = safeindex("name");
  int siz = safeindex("siz");
  int shift = safeindex("shift");

  struct Results
  {
    int siz;
    int shift;
  };
  map<string, Results> results;
  for (csv::CSVRow& row: rescsv) { // Input iterator
    results.insert({row[name].get<string>(), {row[siz].get<int>(), row[shift].get<int>()}});
  }
  cout<<"Got "<<results.size()<<" fit results"<<endl;
  
  ofstream chagstats("genehisto.csv");
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

  int codcount[4][4][4]={};
  
  ofstream codgccsv("codongc.csv");

  // ggcfrac the g/(g+c) ratio
  // gfrac = fraction of g within genes (read from their 5' to the 3')
  
  // leadgfrac is the g fraction of genes found on the leading strand, with the gene being read in the 5' - 3' direction
  // so this represents if genes themselves have the same GC fraction, dependin on which copying direction they are on
  codgccsv<<"name;ggcfrac;cgcfrac;ttafrac;atafrac;gfrac;cfrac;tfrac;afrac;leadafrac;leadcfrac;leadgfrac;leadtfrac;lagafrac;lagcfrac;laggfrac;lagtfrac"<<endl;
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
	const auto& chr = c.second.chromosome;
	
	if(chr.size() < 1000000) {
	  cout<<c.first<<" too small, "<<chr.size()<<" ntds, skipping"<<endl;
	  continue;
	}

	std::string chromoname = c.first;
	if(!results.count(chromoname)) {
	  cerr<<"Could not find fit results for "<<chromoname<<endl;
	  exit(1);
	}
	int chromolen = results[chromoname].siz;
	int shift = results[chromoname].shift;
	unsigned int leadnucs[4]={}, lagnucs[4]={};
	auto annotations = gar.getAll(c.first);
	cout<<"Got "<<annotations.size()<<" annotations"<<endl;
	// codon bias within genes, so aligned to ORF start
	// codons: positive sense genes, anticodons antisense
	map<NucleotideStore, int> codons, anticodons;

	// sliding codons within genes
	map<NucleotideStore, int> scodons, santicodons;

	// triplet ie sliding codons
	map<NucleotideStore, int> tcodons;

	// amino acid counts
	map<char, int> aacount, antiaacount;


	/*	
	char triplet[4]={};
	triplet[1]=chr.get(1);
	triplet[2]=chr.get(2);
	for(uint64_t s = 2 ; s < chr.size(); ++s) {
	  triplet[0]=triplet[1];
	  triplet[1]=triplet[2];
	  triplet[2] = chr.get(s);
	  
	  tcodons[triplet]++;
	}
	*/
	cout<<"Got "<<tcodons.size()<<" tcodons"<<endl;
	//	for(const auto& cd : tcodons) {
	//	  cout<<cd.first<<" "<<cd.second<<"\n";
	// }
	
	for(const auto& a : annotations) {
	  if(a.type=="gene") {
	    
	    //	    cout<<"name '"<<a.name<<"' strand "<<a.strand<<" gene "<<a.gene<<" type " <<a.type<<" tag " <<a.tag<< " start " <<a.startPos<<" - "<<a.stopPos<<endl;
	    if(a.strand) {
	      for(auto pos = a.startPos; pos < a.stopPos; pos++) {
		auto nuc = chr.getRange(pos-1, 3);
		//		if(pos == a.startPos)
		//		  cout << nuc << endl;
		if(!((pos - a.startPos) % 3)) {
		  //		  aacount[DNAToAminoAcid(nuc.toASCII().c_str())]++;
		  codons[nuc]++;
		  codcount[nuc.getNum(0)][nuc.getNum(1)][nuc.getNum(2)]++;
		}
		scodons[nuc]++;

		int rpos = pos - shift;
		if(rpos < 0)
		  rpos += chromolen;
		else if(rpos > chromolen)
		  rpos -= chromolen;

		if(rpos < 0.5*chromolen) // leading strand
		  leadnucs[nuc.getNum(0)]++;
		else
		  lagnucs[nuc.getNum(0)]++;
		
	      }

	    }
	    else {
 	      for(auto pos = a.stopPos; pos > a.startPos; pos--) {
		auto nuc = chr.getRange(pos-3, 3).getRC();
		//		if(pos == a.stopPos)
		//		  cout << nuc << endl;
		if(!((a.stopPos - pos)%3)) {
		  //		  antiaacount[DNAToAminoAcid(nuc.toASCII().c_str())]++;		  
		  anticodons[nuc]++;
		  codcount[nuc.getNum(0)][nuc.getNum(1)][nuc.getNum(2)]++;
		}
		santicodons[nuc]++;

		int rpos = pos - shift;
		if(rpos < 0)
		  rpos += chromolen;
		else if(rpos > chromolen)
		  rpos -= chromolen;

		if(rpos < 0.5*chromolen) // leading strand
		  leadnucs[nuc.getNum(0)]++;
		else
		  lagnucs[nuc.getNum(0)]++;
	      }
	    }
	  }
	}

	int gCodonCount=0, cCodonCount=0, tCodonCount=0, aCodonCount=0, totCodonCount=0;
	vector<decltype(codons)> poss({codons, anticodons});
	for(const auto& c : poss) {
	  for(const auto& codo : c) {
	    for(int pos = 0; pos < 3 ; ++pos) {
	      if(codo.first[pos]=='G') {
		gCodonCount += codo.second;
	      }
	      else if(codo.first[pos]=='C') {
		cCodonCount += codo.second;
	      }
	      else if(codo.first[pos]=='T') {
		tCodonCount += codo.second;
	      }
	      else if(codo.first[pos]=='A') {
		aCodonCount += codo.second;
	      }
	      totCodonCount += codo.second;
	    }
	  }
	}
	//cout<<gCodonCount<<" " << cCodonCount<<" " << gCodonCount*100.0/(gCodonCount+cCodonCount) << "% - "<< cCodonCount*100.0/(gCodonCount+cCodonCount)<<endl;
	codgccsv<<c.first<<";"<<1.0*gCodonCount/(gCodonCount+cCodonCount) << ";"<< 1.0*cCodonCount/(gCodonCount+cCodonCount)<<";"<<1.0*tCodonCount/(tCodonCount+aCodonCount) << ";"<< 1.0*aCodonCount/(tCodonCount+aCodonCount)<<";";
	codgccsv<<1.0*gCodonCount/totCodonCount<<";";
	codgccsv<<1.0*cCodonCount/totCodonCount<<";";
	codgccsv<<1.0*tCodonCount/totCodonCount<<";";
	codgccsv<<1.0*aCodonCount/totCodonCount;
	int totleadnucs=0, totlagnucs=0;
	for(int p =0 ; p < 4; ++p) {
	  totleadnucs += leadnucs[p];
	  totlagnucs += lagnucs[p];
	}
	for(int p =0 ; p < 4; ++p)
	  codgccsv<<";"<<1.0*leadnucs[p]/totleadnucs;
	for(int p =0 ; p < 4; ++p)
	  codgccsv<<";"<<1.0*lagnucs[p]/totlagnucs;

	codgccsv<<endl;

	ofstream codoncsv(c.first+"_codons.csv");
	codoncsv<<"codon,sense,antisense,ssense,santisense,total,antitotal,dup,gcscore,tascore"<<endl;
	set<NucleotideStore> done;
	for(const auto& codo : codons) {
	  codoncsv<<codo.first<<","<<codo.second<<","<<anticodons[codo.first];
	  codoncsv<<","<<scodons[codo.first]<<","<<santicodons[codo.first];
	  codoncsv<<","<<tcodons[codo.first]<<","<<tcodons[codo.first.getRC()]<<","<<done.count(codo.first.getRC());
	  int gcscore=0, tascore=0;
	  for(int n=0; n < 3; ++n) {
	    if(codo.first.get(n)=='G')
	      gcscore++;
	    else if(codo.first.get(n)=='C')
	      gcscore--;
	    else if(codo.first.get(n)=='T')
	      tascore++;
	    else if(codo.first.get(n)=='A')
	      tascore--;
	  }
	  codoncsv<<","<<gcscore<<","<<tascore;
	  codoncsv<<"\n";
	  done.insert(codo.first);
	}

	ofstream aacsv(c.first+"_aa.csv");
	aacsv<<"aa,sense,antisense"<<endl;

	for(const auto& aa : aacount) {
	  aacsv<<aa.first<<","<<aa.second<<","<<antiaacount[aa.first]<<"\n";
	}

	ofstream coddot(c.first+".dot");
	coddot <<"digraph D {\n";
	for(int x=0; x < 4; ++x) {
	  int totX=0;
	  for(int y=0; y < 4; ++y) {
	    int totY=0;
	    for(int z=0; z < 4; ++z) {
	      int totZ=codcount[x][y][z];
	      //	      coddot<< ("ACGT"[x]) << ("ACGT"[y]) << ("ACGT"[z]) <<": "<<totZ<<endl;
	      coddot<< ("ACGT"[x]) << ("ACGT"[y]) <<" -> "<<("ACGT"[x]) << ("ACGT"[y]) << ("ACGT"[z]) <<" [ label=\" "<<totZ<<"\"]\n";
	      totY+=totZ;
	      string codon({"ACGT"[x], "ACGT"[y], "ACGT"[z]});
	      coddot  << codon << " -> " << "AA_"<<DNAToAminoAcid(codon.c_str())<<endl;
	    }
	    //	    coddot<< ("ACGT"[x]) << ("ACGT"[y]) <<": "<<totY<<endl;
	    coddot<< ("ACGT"[x]) <<" -> " << ("ACGT"[x]) << ("ACGT"[y]) <<" [ label=\" "<<totY<<"\" ] \n";
	    totX += totY;
	  }
	  //  coddot<< ("ACGT"[x]) <<": "<<totX<<endl;
	  //	  coddot<< "start -> "<< ("ACGT"[x]) <<"\n";
	}

	
	
	coddot<<"}\n";
	
      }
    }
    catch(std::exception& e) {
      cerr<<"Exception: "<<e.what()<<endl;
    }
  }
}


