#include <boost/test/unit_test.hpp>
#include "dnamisc.hh"
#include "misc.hh"
#include <set>
using namespace std;

BOOST_AUTO_TEST_SUITE(misc_hh)

BOOST_AUTO_TEST_CASE(test_kmerMapper) {
  BOOST_CHECK_EQUAL(kmerMapper("AAAA", 0, 4), 0U);
  BOOST_CHECK_EQUAL(kmerMapper("AAAAAAAA", 0, 8), 0U);
  BOOST_CHECK_EQUAL(kmerMapper("AAAAAAAAAAAA", 0, 12), 0U);
  BOOST_CHECK_EQUAL(kmerMapper("AAAAAAAAAAAAAAAA", 0, 16), 0U);

  BOOST_CHECK_EQUAL(kmerMapper("CCCC", 0, 4), 85U);
  BOOST_CHECK_EQUAL(kmerMapper("CCCCCCCC", 0, 8), 21845U);
  BOOST_CHECK_EQUAL(kmerMapper("CCCCCCCCCCCC", 0, 12), 5592405U);
  BOOST_CHECK_EQUAL(kmerMapper("CCCCCCCCCCCCCCCC", 0, 16), 1431655765U);

  BOOST_CHECK_EQUAL(DNAToAminoAcid("GCC"), 'A');
  BOOST_CHECK_EQUAL(AminoAcidName('A'), "Alanine");
}

BOOST_AUTO_TEST_CASE(test_visitAllNgrams) {
  set<string> all;
  auto func = [&all](const std::string& ngram) {
    all.insert(ngram);
  };
  visitAllNgrams(func, 3);
  BOOST_CHECK_EQUAL(all.size(), 64);
  all.clear();
  visitAllNgrams(func, 1);
  BOOST_CHECK_EQUAL(all.size(), 4);
  BOOST_CHECK_EQUAL(all.count("A"), 1);
  BOOST_CHECK_EQUAL(all.count("C"), 1);
  BOOST_CHECK_EQUAL(all.count("G"), 1);
  BOOST_CHECK_EQUAL(all.count("T"), 1);

  all.clear();
  visitAllNgrams(func, 6);
  BOOST_CHECK_EQUAL(all.size(), 4096);
  BOOST_CHECK_EQUAL(all.count("AAACCC"), 1);
  BOOST_CHECK_EQUAL(all.count("TTTTTT"), 1);

  
  
}

BOOST_AUTO_TEST_SUITE_END()
