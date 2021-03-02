#pragma once
#include <string>
#include <vector>
#include <boost/utility/string_ref.hpp>
#include <string.h>
#include <iostream>
extern "C" {
#include "hash.h"
}

class NucleotideStore
{
public:
  explicit NucleotideStore(const boost::string_ref& in)
  {
    append(in);
  }
  NucleotideStore() {}
  void append(char c);
  void append(const boost::string_ref& line);

  int getNum(size_t pos) const
  {
    uint8_t byte;
    if(pos/4 < d_storage.size())
      byte=d_storage.at(pos/4);
    else
      byte=d_curval;
    
    byte >>= ((pos%4)*2);
    return (byte & 0x3);
  }

  char get(size_t pos) const
  {
    return "ACGT"[getNum(pos)];
  }
  
  
  char operator[](size_t pos) const
  {
    return get(pos);
  }
  void set(size_t pos, char c);
  NucleotideStore getRange(size_t pos, size_t len) const;
  NucleotideStore getRC() const;
  size_t size() const
  {
    return 4*d_storage.size() + bitpos/2;
  }

  struct Delta
  {
    uint32_t pos;
    char o;
    enum class Action {Replace, Delete, Insert} a;
    bool operator==(const Delta& rhs) const
    {
      return pos==rhs.pos && o==rhs.o && a==rhs.a;
    }
  };

  std::vector<Delta> getDelta(const NucleotideStore& b, double mispen=1, double gappen=2, double skwpen=0) const;
  void applyDelta(std::vector<Delta>& delta);
  size_t hash() const
  {
    /*
    if(d_storage.size()==4) {
      uint32_t ret;
      memcpy((char*)&ret, d_storage.c_str(), 4);
      return ret;
    }
    */
    return qhash(d_storage.c_str(), d_storage.size(), bitpos ? d_curval : 0);
  }

  size_t overlap(const NucleotideStore& rhs) const;
  size_t fuzOverlap(const NucleotideStore& rhs, int ratio) const;

  bool isCanonical() const
  {
    return (*this < getRC());
  }
  bool operator==(const NucleotideStore& rhs) const
  {
    return d_storage == rhs.d_storage && bitpos == rhs.bitpos && d_curval == rhs.d_curval;
  }

  bool operator<(const NucleotideStore& rhs) const
  {
    if(d_storage < rhs.d_storage)
      return true;
    if(d_storage > rhs.d_storage)
      return false;

    // we have to think about it

    auto ourpos=d_storage.size()*4;
    auto rhspos=rhs.d_storage.size()*4;

    //    std::cerr<<"Ok, thinking about it, "<<ourpos<<", "<<rhspos<<std::endl;
    if(size() == rhs.size()) // optimization
      return d_curval < rhs.d_curval;
    
    for(; ourpos < size() && rhspos < rhs.size(); ++ourpos, ++rhspos) {
      //      std::cerr<<"Is "<<get(ourpos)<<" < " << rhs.get(rhspos)<<"?"<<std::endl;
      if(get(ourpos) < rhs.get(rhspos))
	return true;
      if(get(ourpos) > rhs.get(rhspos))
	return false;
    }
    //    std::cerr<<"Out of things to compare, shortest should now win: "<<ourpos<<", "<<rhspos<<std::endl;
    if(ourpos == size() && rhspos != rhs.size())
      return true;
    return false;
  }
  
  static char getVal(char c)
  {
    switch(c) {
    case 0:
    case 'A':
    case 'a':
      return 0;
      
      
    case 1:
    case 'C':
    case 'c':
      return 1;
      
    case 2:
    case 'G':
    case 'g':
      return 2;
      
    case 3:
    case 'T':
    case 't':
      return 3;
    default:
      return 0;
    }
    throw std::runtime_error("Impossible nucleotide: "+std::string(1, c));
  }

  std::string getString() const { return d_storage; }
  uint32_t getUInt32() const {
    if(size() != 16)
      throw std::runtime_error("GetUInt32 only works on 16 nucleotides");
    uint32_t ret;
    memcpy(&ret, d_storage.c_str(), 4);
    return ret;
  }
  void setString(const std::string& str) { d_storage = str; }
  std::string toASCII() const;
private:
  uint8_t d_curval{0};
  uint8_t bitpos{0};
  std::string d_storage;
};

std::ostream& operator<<(std::ostream& os, const NucleotideStore& ns);
std::ostream& operator<<(std::ostream& os, const NucleotideStore::Delta& delta);
