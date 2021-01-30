#include "misc.hh"
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <set>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <IntervalTree/IntervalTree.h>

/* chromo  source  type    start   stop    ?       strand  ?       key "val" ; key "val";
   1       havana  gene    11869   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene";
*/

using namespace std;

std::string_view getNextWord(const std::string& line, string::size_type* pos)
{
  
  while(line.at(*pos)==' ' || line.at(*pos)=='\t')
    (*pos)++;
  auto startpos = *pos;
  while(line.at(*pos)!=' ' && line.at(*pos)!='\t')
    (*pos)++;

  return std::string_view(&line[startpos], *pos-startpos);
}


std::string_view getNextQuotedWord(const std::string& line, string::size_type* pos)
{
  
  while(line.at(*pos)==' ' || line.at(*pos)=='\t')
    (*pos)++;

  if(line.at(*pos)!='"')
    throw std::runtime_error("Quoted word was not quoted: " + line.substr(*pos));
  (*pos)++;
  auto startpos = *pos;
  while(line.at(*pos)!='"')
    (*pos)++;

  (*pos)++;
  return std::string_view(&line[startpos], *pos-startpos-1);
}


struct GeneData
{
  std::string chromo;
  uint32_t start, stop;
  bool strand;
  // named & versioned:
  std::map<string, std::map<int, vector<pair<uint32_t,uint32_t>>>> exons;
};

struct RepeatAnnotation
{
  std::string name;
  bool strand;
};

enum class LocStatus { gene = 1, exon = 2, cds=4, utr=8, decayedexon=16, repeat=32, alu=64, l1=128, rna=256 };

map<std::string, IntervalTree<unsigned int, RepeatAnnotation> > g_ras;

struct ChromoStatus
{
  void setBit(uint32_t pos, LocStatus ls)
  {
    if(pos >= d_stat.size())
      d_stat.resize(pos+1);
    d_stat[pos] |= (uint16_t) ls;
  }

  void setBits(uint32_t start, uint32_t stop, LocStatus ls)
  {
    if(stop >= d_stat.size())
      d_stat.resize(stop+1);
    for(auto pos = start; pos != stop; ++pos)
      d_stat[pos] |= (uint16_t)ls;
  }
  
  bool bitSet(uint32_t pos, LocStatus ls)
  {
    return d_stat.at(pos) & (uint16_t)ls;
  }
  size_t size() const
  {
    return d_stat.size();
  }
  vector<uint16_t> d_stat;
};

std::map<string, ChromoStatus> g_chromostatus;
void readRepeats(string_view fname)
{
  map<std::string, vector<Interval<unsigned int, RepeatAnnotation> > > ras;
  
  shared_ptr<FILE> fp(fopen(&fname[0], "r"), fclose);
  std::string line;
  unsigned int lineno=0;

  // #bin    swScore milliDiv        milliDel        milliIns        genoName        genoStart       genoEnd genoLeft        strand  repName repClass        repFamily       repStart        repEnd  repLeft id

  while(stringfgets(fp.get(), &line)) {
    ++lineno;
    if(lineno==1)
      continue;

    string::size_type pos = 0;
    getNextWord(line, & pos); // bin
    getNextWord(line, &pos); //swScore
    getNextWord(line, &pos); // milliDiv
    getNextWord(line, &pos); // milliDel
    getNextWord(line, &pos); // milliIns
    auto chromoName =     getNextWord(line, &pos);
    auto start = (unsigned int)atoi(&getNextWord(line, &pos)[0]);
    auto stop = (unsigned int)atoi(&getNextWord(line, &pos)[0]);
    getNextWord(line, &pos); // genoLeft
    bool strand =     getNextWord(line, &pos)=="+";
    auto repName =     getNextWord(line, &pos);
    auto repClass [[maybe_unused]] =     getNextWord(line, &pos);
    auto repFamily [[maybe_unused]] = getNextWord(line, &pos);

    RepeatAnnotation ra;
    ra.name = repName;
    ra.strand = strand;
    string shorter = boost::replace_first_copy((string)chromoName, "chr", "");
    
    ras[shorter].push_back(Interval<unsigned int, RepeatAnnotation>(start, stop, ra));

    auto& genestatus = g_chromostatus[shorter];

    bool alu = (repName.find("Alu")==0);
    bool l1 = (repName.find("L1")==0);
    for(auto pos = start; pos != stop; ++pos) {
      genestatus.setBit(pos, LocStatus::repeat);
      if(alu)
        genestatus.setBit(pos, LocStatus::alu);
      if(l1)
        genestatus.setBit(pos, LocStatus::l1);
    }

    
  }
  for(auto& ra : ras) {
    g_ras[ra.first] = IntervalTree<unsigned int, RepeatAnnotation>(std::move(ra.second));
  }
  cout<<"Read "<<lineno<<" repeats from "<<g_ras.size()<<" chromosomes"<<endl;

}


uint32_t binCount(const auto& store, std::initializer_list<LocStatus> mustlist={}, std::initializer_list<LocStatus> mustNotlist={})
{
  uint32_t ret=0;
  uint16_t must=0, mustNot=0;
  for(const auto& m : mustlist)
    must |= (uint16_t)m;
  for(const auto& m : mustNotlist)
    mustNot |= (uint16_t)m;
  
  for(unsigned int i =0 ; i < store.size(); ++i) {
    if((i & must) != must)
      continue;
    if(i & mustNot)
      continue;
    ret += store[i];
  }
  return ret;
}

int main(int arvg, char **argv)
{
  readRepeats(argv[1]);
  std::string line;
  ofstream csv("genes.csv");
  ofstream exonscsv("exons.csv");

  std::map<string, std::map<string, GeneData>> chromogenes;

  csv << "chromo name numexons exonsize size start stop strand\n";
  exonscsv << "chromo name start stop strand\n";
  while(stringfgets(stdin, &line)) {
    if(auto pos = line.find('#'); pos != string::npos)
      line.resize(pos);
    if(line.empty())
      continue;

    string::size_type pos = 0;
    auto chromo=getNextWord(line, &pos);
    getNextWord(line, &pos); // source
    auto type=getNextWord(line, &pos);
    auto start = (unsigned int) atoi(&getNextWord(line, & pos)[0]);
    auto stop = (unsigned int) atoi(&getNextWord(line, & pos)[0]);
    getNextWord(line, &pos); // unknown1
    bool sense = getNextWord(line, &pos)=="+";
    getNextWord(line, &pos); // unknown2 

    string_view rest{&line.at(pos), line.size()-pos-1};

    string_view gene_name, gene_biotype, transcript_name, transcript_biotype;
    int transcript_version=-1;
    while(pos +1  < line.size()) {
      auto key = getNextWord(line, &pos);
      auto val = getNextQuotedWord(line, &pos);
      if(key=="gene_name")
        gene_name = val;
      else if(key=="transcript_name")
        transcript_name = val;
      else if(key=="transcript_biotype")
        transcript_biotype = val;

      else if(key=="transcript_version")
        transcript_version = atoi(&val[0]);
      else if(key=="gene_biotype")
        gene_biotype = val;
      pos+=1;
    }

    if(type=="gene" && gene_biotype.find("RNA")!=string::npos) {
      auto& genestatus = g_chromostatus[(string)chromo];
      genestatus.setBits(start, stop, LocStatus::rna);
    }
    
    if(gene_biotype != "protein_coding")
      continue;
                 
    
    if(type=="gene") {
      auto& gene = chromogenes[(string)chromo][(string)gene_name];
      gene.start = start;
      gene.stop = stop;
      gene.strand = sense;
      gene.chromo = chromo;

      auto& genestatus = g_chromostatus[(string)chromo];
      genestatus.setBits(start, stop, LocStatus::gene);
      
    }
    if(type=="CDS") {
      auto& genestatus = g_chromostatus[(string)chromo];
      genestatus.setBits(start, stop, LocStatus::cds);
      
    }
    if(type.find("utr") != string::npos) {
      auto& genestatus = g_chromostatus[(string)chromo];
      genestatus.setBits(start, stop, LocStatus::utr);
    }
    if(type=="exon") {
      //cout<<chromo<<"\t"<<sense<<"\t" <<start<<"\t"<<stop<<"\t"<<rest<<endl;

    // gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"

      auto& genestatus = g_chromostatus[(string)chromo];

      auto& gene = chromogenes[(string)chromo][(string)gene_name];
      LocStatus val;
      if(transcript_biotype != "nonsense_mediated_decay") {
        val = LocStatus::exon;
        gene.exons[(string)transcript_name][transcript_version].push_back({start, stop});
      }
      else {
        val = LocStatus::decayedexon;;
        gene.exons[boost::algorithm::to_lower_copy((string)transcript_name)][transcript_version].push_back({start, stop});
      }
      genestatus.setBits(start, stop, val);

    }
  }

  uint32_t totintergenicct=0, totgenect=0, totexonct=0, totrepeatct=0;
  uint32_t totexorepeatct=0, totintronrepeatct=0, totinterrepeatct=0, totdecayexonct=0;
  uint32_t totintronaluct=0, totintergenaluct=0;
  uint32_t totintronl1ct=0, totintergenl1ct=0;
  uint32_t totct = 0, totcdsct = 0, totutrct = 0;
  uint32_t totrnact=0, totrnarepeatct=0;

  std::set<string> mustdo({"1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
        "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
        "21", "22", "X", "Y", "MT"});


  vector<uint32_t> counts(256*256);
  for(const auto& chromo : g_chromostatus) {
    cout<<"Chromosome "<<chromo.first<<" goes up to "<<chromo.second.size()<<endl;
    if(!mustdo.count(chromo.first)) {
      continue;
    }
    
    int intergenicct=0, genect=0, exonct=0, repeatct=0;
    int exorepeatct=0, intronrepeatct=0, interrepeatct=0, decayexonct=0;
    int intronaluct=0, intergenaluct=0;
    int intronl1ct=0, intergenl1ct=0;
    int rnact=0, rnarepeatct=0, cdsct = 0, utrct = 0;
    for(const auto& s : chromo.second.d_stat) {
      ++counts[s];
      
      if(s & (uint16_t) LocStatus::gene)
        ++genect;
      else if(s & (uint16_t) LocStatus::rna)
        ++rnact;
      else {
        ++intergenicct;
        if(s & (uint16_t) LocStatus::repeat) {
          ++interrepeatct;
          
          if(s & (uint16_t) LocStatus::alu)
            ++intergenaluct;
          if(s & (uint16_t) LocStatus::l1)
            ++intergenl1ct;

        }

      }

      if(s & (uint16_t) LocStatus::exon) { // exon
        ++exonct;
        if(s&4)
          ++exorepeatct;
      }

      if(s & (uint16_t) LocStatus::cds) {
        ++cdsct;
      }
      if(s & (uint16_t) LocStatus::utr) {
        ++utrct;
      }
      if(s & (uint16_t) LocStatus::repeat) {
        ++repeatct;
      }
      if(s & (uint16_t) LocStatus::decayedexon)
        ++decayexonct;

      // "a decayed exon is an intron"
      if((s & (uint16_t) LocStatus::gene) && !(s & (uint16_t) LocStatus::exon) &&  (s & (uint16_t) LocStatus::repeat)) {
        ++intronrepeatct;
        if(s & (uint16_t) LocStatus::alu)
          ++intronaluct;
        if(s & (uint16_t) LocStatus::l1)
          ++intronl1ct;
      }

      if((s & (uint16_t) LocStatus::rna) && (s & (uint16_t) LocStatus::repeat))
        ++rnarepeatct;
      
    }
    cout<<chromo.second.size()<<" ntds\n";
    cout<<intergenicct<<" intergenic ntds ("<< (intergenicct*100.0/chromo.second.size()) << "%)\n";
    cout<<"  "<<interrepeatct<<" intergenic repeat ntds ("<< (interrepeatct*100.0/intergenicct) << "%)\n";
    cout<<genect<< " gene ntds ("<< (genect*100.0/chromo.second.size()) << "%)\n";
    cout<<"  "<<exonct <<" exon ntds (" << (exonct*100.0/genect) <<"%, "<<(exonct*100.0/chromo.second.size())<<"%)\n";
    cout<<"    " << exorepeatct << " exonic repeats ntds (" << (exorepeatct*100.0/exonct)<< "%)\n";
    cout<<"  "<<decayexonct<<" decay mediated exon ntds (" << (decayexonct*100.0/exonct) << "%)\n";
    cout<<"  "<<(genect-decayexonct-exonct) <<" intron ntds (" << ( (genect-exonct)*100.0/genect) <<"%)\n";
    cout<<"    " << intronrepeatct << " intronic repeats ntds (" << (intronrepeatct*100.0/(genect-exonct))<< "%)\n";
    /*
198234321
100264834 intergenic, 
 57464310 intergeneic repeats, 
 97969487 part of gene, 
  5882277 part of exon, 
   810903 exonic repeats, 
 43342843 intronic repeats, 
    */

    if(mustdo.count(chromo.first)) {
      cout<<chromo.first<<" "<<chromo.second.size()<<" total, "<<intergenicct<<" intergenic"<<", "<< genect <<" part of gene, "<<exonct<<" part of exon, "<< repeatct <<" part of repeat, "<< exorepeatct<<" exonic repeats, "<<intronrepeatct<<" intronic repeats, "<< interrepeatct<<" intergeneic repeats, "<<decayexonct<<" part of decayed exon"<<endl;
      
      totct += chromo.second.size();
      
      totintergenicct += intergenicct;
      totgenect += genect; totexonct += exonct;
      totrepeatct += repeatct;
      
      totexorepeatct += exorepeatct; totintronrepeatct += intronrepeatct;
      totinterrepeatct += interrepeatct;
      totdecayexonct += decayexonct;
      totintronaluct += intronaluct;
      totintergenaluct += intergenaluct;
      totintronl1ct += intronl1ct;
      totintergenl1ct += intergenl1ct;
      totrnact += rnact;
      totrnarepeatct += rnarepeatct;
      totcdsct += cdsct;
      totutrct += utrct;
    }
  }
  boost::format fmt(R"([
    ('.', %d, [
        ('intergenic', %d, [
            ('Repeats', %d, [
               ('Alu', %d, []),
               ('L1', %d, []) 
        ]),
        ]),
        ('ncRNA', %d, [
           ('repeats', %d, [
            ('Alu', %d, []),
            ('L1', %d, [])
         ]) ]),
        ('genetic', %d, [
            ('exonic', %d, [
                ('cds', %d, []),
                ('utr', %d, []),
            ]),
            ('intronic', %d, [
                ('repeats', %d, [
                  ('Alu', %d, []),
                  ('L1', %d, []) 
               ])
            ]),
            ('decayed', %d, []),
        ]),
    ]),
])");


  auto bc = [&counts](const std::initializer_list<LocStatus>& a={}, const std::initializer_list<LocStatus> &b={}) {
    return binCount(counts, a, b);
  };
  cout << ( fmt % bc() % bc({}, {LocStatus::gene, LocStatus::rna}) %
            bc({LocStatus::repeat}, {LocStatus::gene, LocStatus::rna}) %
            bc({LocStatus::alu}, {LocStatus::gene, LocStatus::rna}) %
            bc({LocStatus::l1}, {LocStatus::gene, LocStatus::rna}) %
            bc({LocStatus::rna}) % bc({LocStatus::rna, LocStatus::repeat}) %
            bc({LocStatus::rna, LocStatus::alu}) % bc({LocStatus::rna, LocStatus::l1}) %
            bc({LocStatus::gene}) %
            bc({LocStatus::gene, LocStatus::exon}, {LocStatus::decayedexon}) % // exon
            bc({LocStatus::cds}, {LocStatus::utr}) % // cds
            bc({LocStatus::utr}, {LocStatus::cds}) % // utr
            bc({LocStatus::gene}, {LocStatus::exon, LocStatus::decayedexon}) % // "intron"
            bc({LocStatus::gene, LocStatus::repeat}, {LocStatus::exon}) %
            bc({LocStatus::gene, LocStatus::alu}, {LocStatus::exon}) %
            bc({LocStatus::gene, LocStatus::l1}, {LocStatus::exon}) %
            bc({LocStatus::gene, LocStatus::decayedexon}, {LocStatus::exon})) << endl;
            
            

  
  int printed=0;

  for(unsigned int i = 0 ; i < counts.size(); ++i) {
    if(!counts[i])
      continue;
    if(!(printed % 25))
      cout<<"\tgene\texon\tcds\tutr\tdecayed\trepeat\talu\tl1\trna\n";
    ++printed;
    
    cout << i <<"\t";
    for(int j = 0 ; j < 9; ++j) {
      if((1 << j) & i)
        cout << "*";
      cout<<"\t";
    }
    cout<<counts[i]<<endl;
  }
  
  return  0;
  cout<<"Have "<<chromogenes.size()<<" chromosomes"<<endl;
  for(auto& chromo : chromogenes) {
    cout<<"Chromosome "<<chromo.first<<" has "<<chromo.second.size()<<" gene names"<<endl;
    for(auto& gene : chromo.second) {
      unsigned int genesize = 1 + gene.second.stop - gene.second.start;
      cout<<gene.first<<" spans "<< genesize <<" ntds & has "<<gene.second.exons.size()<<" transcripts. Strand " << (gene.second.strand  ? '+' : '-') <<endl;
      int totexonsize=0;

      vector<vector<string>> exonic(genesize);
      
      if(!gene.second.exons.empty()) {
        for(auto& gexons : gene.second.exons) {
          cout <<"   transcript "<<gexons.first<<" has "<<gexons.second.size()<<" versions.\n";
          for(auto& version : gexons.second) {
            cout<< "      version "<<version.first<<" has "<<version.second.size()<<" exons (";
            totexonsize = 0;
            sort(version.second.begin(), version.second.end(), [](const auto& a, const auto& b)
                 {
                   return a.first < b.first;
                 });
            for(auto exon = version.second.begin(); exon != version.second.end(); ++exon) {
              if(exon->first < gene.second.start || exon->second > gene.second.stop)
                cerr<<"No way! "<<gene.first<<" "<<exon->first<<" " <<gene.second.start<<", "<<exon->second<<" " <<gene.second.stop<<endl;
              else {
                for(auto s = exon->first - gene.second.start; s != exon->second - gene.second.start; ++s)
                  exonic.at(s).push_back(gexons.first);
              }
              
              cout<< (exon->second - exon->first) <<"";

              auto results = g_ras[chromo.first].findOverlapping(exon->first, exon->second);
              cout<<"[";
              for(const auto& res : results) {
                cout<<res.value.name<<" ";
              }
              cout<<"] ";

              if(exon + 1 != version.second.end()) {
                results = g_ras[chromo.first].findOverlapping(exon->second, std::next(exon)->first);
                cout<<" <- " << std::next(exon)->first - exon->second<<" ";
                for(const auto& res : results) {
                  cout << (res.value.strand ? '+' : '-');
                  cout<<res.value.name<<" ";
                }
                cout<<"-> ";
                
              }
              
              totexonsize += 1 + exon->second - exon->first;
            }
            cout<<"), "<<totexonsize<<" exonic ntds\n";
          }
        }
        unsigned int touched=0;

        ofstream genefs("./genes/"+gene.first);
        set<string> transcripts;
        for(auto& b : exonic) {
          if(!b.empty()) 
            ++touched;
                    
          std::sort(b.begin(), b.end());
          for(const auto& l : b) 
            transcripts.insert(l);
        }

        auto results = g_ras[chromo.first].findOverlapping(gene.second.start, gene.second.stop);
        genefs<<"Repeats: ";
        for(const auto& res : results) {
          genefs<<res.value.name<<" ";
        }
        genefs<<"\n";
        /*
        for(const auto& transcript: transcripts) {
          genefs << transcript<<": ";
          for(const auto& b : exonic) {
            if(binary_search(b.begin(), b.end(), transcript))
              genefs<<'*';
            else
              genefs<<' ';
          }
          genefs<<"\n";

        }
        genefs<<endl;
        */

        
        csv << gene.second.chromo<<" "<<gene.first << " "<< 0 << " " <<touched <<" " << (gene.second.stop - gene.second.start) <<" "<<gene.second.start <<" " <<gene.second.stop <<" "<<gene.second.strand<<"\n";

      }
    }
  }
}
