#pragma once
#include <map>
#include <string>
#include <vector>
class TaxoReader
{
public:
  explicit TaxoReader(const std::string& fname);
  std::vector<std::string> get(int id) const;
  size_t size() const
  {
    return d_store.size();
  }
private:
  std::map<int, std::vector<std::string> > d_store;
};
