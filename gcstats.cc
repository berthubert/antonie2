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
 * based on this it emits a big skplot.csv file which has all kinds of skews for all chromosomes, posted every n loci
 * also emits a genomes.csv summary file (should be 'chromosomes' really)
 *
 * The GFF is used to determine the codon position, or if we are in a gene or not
 *
 * Also reads taxonomies (fullnamelineage.dmp)
 *
 * This code ignores chromosomes smaller than x
 */

int main(int argc, char**argv)
{
  if(argc < 3) {
    cerr<<"Syntax: gcstats fullnamelineage.dmp reference.fasta [reference.fasta...] "<<endl;
    return EXIT_FAILURE;
  }

  ofstream genomescsv("genomes.csv");
  genomescsv<<"name;fullname;acount;ccount;gcount;tcount;plasmid;taxonid;realm1;realm2;realm3;realm4;realm5;protgenecount;stopTAG;stopTAA;stopTGA;stopXXX;startATG;startGTG;startTTG;startXXX;dnaApos;dnaAsense"<<endl;
  ofstream skplot("skplot.csv");
  skplot<<"name,relpos,abspos,gcskew,taskew,gcskew0,gcskew1,gcskew2,gcskewNG,taskew0,taskew1,taskew2,taskewNG,pospos,gccount,ngcount";
  for(int cpos = 0 ; cpos < 4 ; ++cpos) {
    skplot << ",acounts"<<cpos<<",ccounts"<<cpos<<",gcounts"<<cpos<<",tcounts"<<cpos;
  }
  skplot<<endl;


  string taxofname(argv[1]);
  cout<<"Reading taxonomies from "<<taxofname<<"... ";
  cout.flush();
  TaxoReader tr(taxofname);
  cout<<" got "<<tr.size()<<" entries"<<endl;

  std::atomic<int> ctr=2;
  std::mutex iolock;
  auto doFunc=[&ctr, &iolock, &argc, &argv, & genomescsv, &skplot, &tr]() {
    for(int n=ctr++; n < argc; n=ctr++) {
      try {
        ReferenceGenome rg(argv[n]);

        cout<<"Done reading genome from "<<argv[n]<<", have "<<rg.numChromosomes()<<" chromosomes, "<<
          rg.numNucleotides()<<" nucleotides"<<endl;

        string garname = argv[n];
        boost::replace_all(garname, ".fna", ".gff");

        GeneAnnotationReader gar(garname);
        cout<<"Done with annotations from "<<garname<<" (taxid "<<gar.d_taxonID<<"), got "<<gar.size()<<" of them, for "<<gar.getChromosomes().size()<<" chromosomes: ";

        for(const auto& c : gar.getChromosomes())
          cout<<"'"<<c<<"' ";
        cout<<endl;
    
        for(const auto& c : rg.getAllChromosomes()) {
          const auto& chr = c.second.chromosome;
          if(chr.size() < 150000) {
            cout<<c.first<<" too small, "<<chr.size()<<" ntds, skipping: "<<c.second.fullname<<endl;
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

      
          cout<<"A: "<<aCount<<", ";
          cout<<"C: "<<cCount<<", ";
          cout<<"G: "<<gCount<<", ";
          cout<<"T: "<<tCount<<", ";
          double gc = 1.0*(cCount + gCount)/(aCount + cCount + gCount +tCount);
          cout<<"gc: "<<gc<<endl;

          double gcskew = 0, taskew=0, pospos=0;
          double gcskews[4] = {0,0,0,0};
          double taskews[4] = {0,0,0,0};
          int acounts[4]={}, ccounts[4]={}, gcounts[4]={}, tcounts[4]={};
          int ngcount=0; // positions not in a gene

          int stopTAG=0, stopTAA=0, stopTGA=0, stopXXX=0;
          int startATG=0, startGTG=0, startTTG=0, startXXX=0;
          int protgenecount=0;
      
          GeneAnnotation last;
          last.gene=false;
          unsigned int gccount=0;

          /*
            Alternate start codons: 83% ATG (3542/4284), 14% (612) GTG, 3% (103) TTG
          */

          string lastName;
          int dnaApos=-1;
          int dnaAsense = -1;
          for(uint64_t s = 0 ; s < chr.size(); ++s) {
            unsigned genepos = s+1;
            // this discovers the edges of genes
            if(!(last.gene && last.startPos <= genepos && genepos < last.stopPos)) {
              last.gene=false;
              auto ans= gar.lookup(c.first, genepos);
              if(!ans.empty()) {
                /*
                cout<<"Got "<<ans.size()<<" annotations for "<<genepos<<" on "<<c.first<<": "<<endl;
                for(const auto& a: ans) {              
                  cout<<" " << a.startPos << " - " <<a.stopPos<<": name '"<<a.name<<"', type '"<<a.type<<"', biotype '"<<a.gene_biotype<<"'"<<endl;
                */
              }
              for(const auto& a: ans) {
                if(a.gene==true && a.gene_biotype == "protein_coding") {
                  //              cout<<"At position "<< genepos <<" of "<<c.first<<" we are in gene '"<<a.name<<"' which goes from " << a.startPos <<" to "<< a.stopPos <<", strand: "<<a.strand<<". Offset: "<< genepos - a.startPos;
                  if(a.name == lastName)
                    continue;
                  lastName=a.name;
                  if(a.name=="dnaA") {
                    dnaApos = a.strand ? a.startPos : a.stopPos;
                    dnaAsense = a.strand;
                  }
              
                  protgenecount++;
                  if(a.strand) { // positive sense
                    auto start = chr.getRange(a.startPos-1, 3);
                    //                cout<<" Start codon: "<<start;
                    if(auto as = start.toASCII();   as != "ATG" && as != "GTG" && as != "TTG") { // common start codons
                      //                  cout<<" !!!";
                      cout<<"Gene "<<a.name<<" of "<<c.first <<" had a weird start codon: "<<start<<endl;
                      startXXX++;
                    }
                    else {
                      if(as == "ATG") startATG++;
                      else if(as == "GTG") startGTG++;
                      else if(as == "TTG") startTTG++;
                    }
                    auto stop = chr.getRange(a.stopPos -3 , 3);
                    //                cout<<" STOP "<<stop;
                    if(auto as = stop.toASCII();   as != "TAG" && as!= "TAA" && as != "TGA") {
                      stopXXX++;
                      cout<<"Gene "<<a.name<<" of "<<c.first<<" had a weird stop codon: "<<stop<<endl;
                      //                  cout<<" ???";
                    }
                    else {
                      if(as == "TAG") stopTAG++;
                      else if(as == "TAA") stopTAA++;
                      else if(as == "TGA") stopTGA++;
                    }
                    //                cout<<"\n";
                  }
                  else { // negative sense
                    auto start = chr.getRange(a.stopPos-3, 3).getRC();
                    //                cout<<" Start codon: "<<start;
                    if(auto as = start.toASCII();   as != "ATG" && as != "GTG" && as != "TTG") { // common start codons
                      //                  cout<<" !!!";
                      cout<<"Gene "<<a.name<<" of "<< c.first<<" had a weird start codon: "<<start<<endl;
                      startXXX++;
                    }
                    else {
                      if(as == "ATG") startATG++;
                      else if(as == "GTG") startGTG++;
                      else if(as == "TTG") startTTG++;
                    }

                    auto stop = chr.getRange(a.startPos-1, 3).getRC();
                    //                cout<<" STOP "<< stop;
                    if(auto as = stop.toASCII();   as != "TAG" && as != "TAA" && as != "TGA") { 
                      //                  cout<<" ???";
                      cout<<"Gene "<<a.name<<" of "<<c.first<<" had a weird stop codon: "<<stop<<endl;
                      stopXXX++;
                    }
                    else {
                      if(as == "TAG") stopTAG++;
                      else if(as == "TAA") stopTAA++;
                      else if(as == "TGA") stopTGA++;
                  
                    }

                    //                cout<<"\n";
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
            else
              ngcount++;
	    
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
              ostringstream str;

              str<<c.first<<","<<1.0*s/chr.size()<<","<<s<<","<<gcskew<<","<<taskew<<","<<gcskews[0]<<","<<gcskews[1]<<","<<gcskews[2]<<"," <<gcskews[3]<<","<<taskews[0]<<","<<taskews[1]<<","<<taskews[2]<<"," <<taskews[3]<<","<<pospos<<","<<gccount<<","<<ngcount;
              for(int cpos = 0 ; cpos < 4 ; ++cpos) {
                str << "," << acounts[cpos]<<","<<ccounts[cpos]<<","<<gcounts[cpos]<<","<<tcounts[cpos];
              }
              str<<"\n";

              std::lock_guard<std::mutex> m(iolock);
              skplot<<str.str();
            }
          }

          vector<string> taxo= tr.get(gar.d_taxonID);
          taxo.resize(6); // so the printing is safe
          bool plasmid = (c.second.fullname.find("plasmid") != string::npos);
          {
            std::lock_guard<std::mutex> m(iolock);
            genomescsv<<c.first<<";"<<boost::replace_all_copy(c.second.fullname, ";", ",")<<";"<<aCount<<";"<<cCount<<";"<<gCount<<";"<<tCount<<";"<<plasmid<<";";
            genomescsv<<gar.d_taxonID<<";"<<taxo[1]<<";"<<taxo[2]<<";"<<taxo[3]<<";"<<taxo[4]<<";"<<taxo[5]<<";"<<protgenecount<<";"<<stopTAG<<";"<<stopTAA<<";"<<stopTGA<<";"<<stopXXX<<";"<<startATG<<";"<<startGTG<<";"<<startTTG<<";"<<startXXX<<";"<<dnaApos<<";"<<dnaAsense<<endl;
          }

        }
      }
      catch(std::exception& e) {
        cerr<<"Exception: "<<e.what()<<endl;
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

