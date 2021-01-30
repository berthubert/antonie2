#include "taxoreader.hh"
#include "misc.hh"
#include <regex>
#include <iostream>
#include <boost/algorithm/string.hpp>

using namespace std;

vector<std::string_view> vstringtok (const string_view in,
           const char * const delimiters = " \t\n")
{
  const string::size_type len = in.length();
  string::size_type i = 0;
  vector<string_view> container;
  while (i<len) {
    // eat leading whitespace
    i = in.find_first_not_of (delimiters, i);
    if (i == string::npos)
      return container;   // nothing left but white space
    
    // find the end of the token
    string::size_type j = in.find_first_of (delimiters, i);
    
    // push token
    if (j == string::npos) {
      container.push_back (string_view(&in[0]+i, len-i));
      return container;
    } else
      container.push_back (string_view(&in[0]+i, j-i));
    
    // set up for next loop
    i = j + 1;
  }
  return container;
}


TaxoReader::TaxoReader(const std::string& fname)
{
  FILE* fp=fopen(fname.c_str(), "rb");
  if(!fp)
    throw runtime_error("Unable to open '"+fname+"' for taxonomy reading");

  std::string line;
  int id;
  int counter=0;
  while(stringfgets(fp, &line)) {
    vector<std::string_view> parts = vstringtok(line, "|");
    if(parts.size() < 3)
      continue;
    id=atoi(&(parts[0][0]));
    //    cout<<"\n\tTaxonomy: ";
    if(parts.size() > 2) {
      vector<std::string_view> parts2 = vstringtok(parts[2], ";");
      vector<string> store;
      store.reserve(parts2.size());
      for(const auto& p2: parts2) {
	//	cout<<"'"<<p2<<"', ";
	string r = boost::trim_copy((string)p2);
	if(!r.empty())
	  store.push_back(r);
      }
      d_store[id]=store;
    }
  }
}

// 562	|	Escherichia coli	|	cellular organisms; Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Escherichia; 	|

vector<string> TaxoReader::get(int id) const
{
  vector<string> ret;

  if(auto iter = d_store.find(id); iter != d_store.end())
    return iter->second;
  return ret;
}
