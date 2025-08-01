#include "misc.hh"
#include <string.h>
#include <stdexcept>
#include <algorithm>
using namespace std;
#include <sys/types.h>
#include <sys/stat.h>
#include <boost/lexical_cast.hpp>
#include <fstream>

//! read a line of text from a FILE* to a std::string, returns false on 'no data'
bool stringfgets(FILE* fp, std::string* line)
{
  char buffer[1024];   
  line->clear();
  
  do {
    if(!fgets(buffer, sizeof(buffer), fp))
      return !line->empty();
    
    line->append(buffer);
  } while(!strchr(buffer, '\n'));
  return true;
}

//! read a line of text from a gzfile to a std::string, returns false on 'no data'
bool stringfgets(gzFile fp, std::string* line)
{
  char buffer[1024];   
  line->clear();
  
  do {
    if(!gzgets(fp, buffer, sizeof(buffer)))
      return !line->empty();
    
    line->append(buffer);
  } while(!strchr(buffer, '\n'));
  return true;
}


uint64_t filesize(const char* name)
{
  struct stat buf;
  if(!stat(name, &buf)) {
    return buf.st_size;
  }
  return 0;
}



char* sfgets(char* p, int num, FILE* fp)
{
  char *ret = fgets(p, num, fp);
  if(!ret) 
    throw std::runtime_error("Unexpected EOF");
  return ret;
}

void chomp(char* line)
{
  char *p;
  p = strchr(line, '\r');
  if(p)*p=0;
  p = strchr(line, '\n');
  if(p)*p=0;
}
#if 0
// thanks jeff sipek
static void rev_and_comp_tbl_small(const char *in, char *out, size_t len)
{
        const char tbl[8] = {
                [1] = 'T',
                [4] = 'A',
                [3] = 'G',
                [7] = 'C',
        };

        if (!in || !out || !len)
                return;

        out[len] = '\0';

        while (len) {
                *out = tbl[(int) in[len - 1] & 0x7];

                len--;
                out++;
        }
}
#endif

void reverseNucleotides(std::string* nucleotides)
{
  std::reverse(nucleotides->begin(), nucleotides->end());
  for(string::iterator iter = nucleotides->begin(); iter != nucleotides->end(); ++iter) {
    if(*iter == 'C')
      *iter = 'G';
    else if(*iter == 'G')
      *iter = 'C';
    else if(*iter == 'A')
      *iter = 'T';
    else if(*iter == 'T')
      *iter = 'A';
  }
}

string compilerVersion()
{
#if defined(__clang__)
  return string("clang " __clang_version__);
#elif defined(__GNUC__)
  return string("gcc " __VERSION__);
#elif defined(_MSC_FULL_VER)
  return string("Microsoft Visual Studio " + boost::lexical_cast<string>(_MSC_FULL_VER));
#else  // add other compilers here 
  return string("Unknown compiler");
#endif
}

// blah 1.fna @file-with-more-fnas
// expands this with the contents of 'file-with-more-fnas', one per line
vector<string> expandArguments(int argc, char** argv)
{
  vector<string> todo;
  for(int n=0; n < argc; ++n) {
    if(argv[n][0]=='@') {
      ifstream ifs(argv[n]+1);
      string fname;
      while(getline(ifs, fname)) {
	if(!fname.empty() && fname[fname.size()-1]=='\n')
	  fname.resize(fname.size()-1);
	todo.push_back(fname);
      }
    }
    else todo.push_back(argv[n]);
  }

  return todo;
}


void visitAllNgrams(std::function<void(const std::string&)> exec, unsigned int chars, std::string start)
{
  if(!chars) {
    exec(start);
    return;
  }
  --chars;
  visitAllNgrams(exec, chars, start+"A");
  visitAllNgrams(exec, chars, start+"C");
  visitAllNgrams(exec, chars, start+"G");
  visitAllNgrams(exec, chars, start+"T");
}
