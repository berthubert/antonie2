#include <iostream>
#include "misc.hh"
#include "fastq.hh"
using namespace std;

bool tryMerge(const FastQRead& one, const FastQRead& two, FastQRead* together)
{
  FastQRead inv;
  const FastQRead* a = &one, *b = &two;

  for(int tries = 0; tries < 2; ++tries) {
    if(tries) {
      a = &two;
      b = &one;
    }

    inv.d_nucleotides = b->d_nucleotides;
    inv.reverse();

    if(inv.d_nucleotides.find(a->d_nucleotides.substr(a->d_nucleotides.length()-20)) == string::npos)
      continue;
  
    for(int overlap = a->d_nucleotides.length() ; overlap > 20; --overlap) {
      if(a->d_nucleotides.substr(a->d_nucleotides.length()-overlap) == inv.d_nucleotides.substr(0, overlap)) {
	//      cerr<<"Got overlap of "<<overlap<<endl;
	//      cerr<<a->d_nucleotides<<endl;
	//      cerr<<string(a->d_nucleotides.length()-overlap,' ')<<inv.d_nucleotides<<endl;
	together->d_nucleotides = a->d_nucleotides;
	together->d_nucleotides += inv.d_nucleotides.substr(overlap);
	//      cerr<<together->d_nucleotides<<endl;
	return true;
      }
    }
  }
  return false;
}

bool findInRead(const FastQRead& fqr, const string& search, const string& rsearch)
{
  auto pos = fqr.d_nucleotides.find(search);
  if(pos != string::npos) {
    cout<<"F: "<<fqr.d_nucleotides.substr(pos)<<endl;
    return true;
  }
  else if((pos=fqr.d_nucleotides.find(rsearch)) != string::npos) {
    FastQRead inv = fqr;
    inv.reverse();
    pos = inv.d_nucleotides.find(search);
    cout<<"R: "<<inv.d_nucleotides.substr(pos)<<endl;
    return true;
  }
  return false;
}

int main(int argc, char**argv)
{
  string search(argv[1]);
  string rsearch(search);
  reverseNucleotides(&rsearch);

  StereoFASTQReader fqreader(argv[2], argv[3], 33);
  FastQRead fqr1, fqr2, merged;

  int mergedCount=0, total=0;
  while(fqreader.getReadPair(&fqr1, &fqr2)) {
    total++;
    if(tryMerge(fqr1, fqr2, &merged)) {
      mergedCount++;
      if(findInRead(merged, search, rsearch))
	cout<<"^ merged"<<endl;
    }
    else {
      findInRead(fqr1, search, rsearch);
      findInRead(fqr2, search, rsearch);
    }
  }
  cout<<"Got "<<total<<" pairs of which "<< mergedCount*100.0/total <<"% could be merged"<<endl;
}


#if 0
#endif
