#include <boost/test/unit_test.hpp>
#include "dnamisc.hh"
#include "nucstore.hh"
#include <iostream>

BOOST_AUTO_TEST_SUITE(nucstore_hh)

BOOST_AUTO_TEST_CASE(test_nucstore_basic) {
  NucleotideStore ns;
  ns.append('A');
  BOOST_CHECK_EQUAL(ns.get(0),'A');
  ns.append('C');
  BOOST_CHECK_EQUAL(ns.get(0),'A');
  BOOST_CHECK_EQUAL(ns.get(1),'C');
  ns.append('G');
  BOOST_CHECK_EQUAL(ns.get(2),'G');
  ns.append('T');

  BOOST_CHECK_EQUAL(ns.get(3),'T');

  BOOST_CHECK_EQUAL(ns.get(0),'A');
  ns.append('C');
  BOOST_CHECK_EQUAL(ns.get(0),'A');
  BOOST_CHECK_EQUAL(ns.get(1),'C');
  ns.append('G');
  BOOST_CHECK_EQUAL(ns.get(2),'G');
  ns.append('T');

  BOOST_CHECK_EQUAL(ns.get(3),'T');


  BOOST_CHECK_EQUAL(ns.size(), 7);
  
  // AACG TAACG

  NucleotideStore sep;
  sep.append("ACGTCGT");
  BOOST_CHECK_EQUAL(ns, sep);


  sep.set(0, 'C');
  BOOST_CHECK_EQUAL(sep.get(0), 'C');


  sep.set(4, 'C');
  BOOST_CHECK_EQUAL(sep.get(0), 'C');
}

BOOST_AUTO_TEST_CASE(test_nucstore_val) {
  NucleotideStore a("ACGTTC");
  BOOST_CHECK_EQUAL(a.get(0), 'A');
  BOOST_CHECK_EQUAL(a.getNum(0), 0);
  BOOST_CHECK_EQUAL(a.getNum(3), 3);
  BOOST_CHECK_EQUAL(a.getNum(4), 3);
  BOOST_CHECK_EQUAL(a.getNum(5), 1);
  
}

BOOST_AUTO_TEST_CASE(test_nucstore_comp) {
  using namespace std;
  NucleotideStore a("AAA"), b("AA");
  BOOST_CHECK_LT(b, a);

  BOOST_CHECK_LT(NucleotideStore("A"), NucleotideStore("C"));
  BOOST_CHECK_LT(NucleotideStore("C"), NucleotideStore("G"));
  BOOST_CHECK_LT(NucleotideStore("G"), NucleotideStore("T"));
  
  BOOST_CHECK_LT(NucleotideStore("A"), NucleotideStore("C"));
  BOOST_CHECK_LT(NucleotideStore("AA"), NucleotideStore("CC"));
  BOOST_CHECK_LT(NucleotideStore("AAA"), NucleotideStore("CCC"));
  BOOST_CHECK_LT(NucleotideStore("AAAA"), NucleotideStore("CCCC"));
  BOOST_CHECK_LT(NucleotideStore("AAAAA"), NucleotideStore("CCCCC"));
  BOOST_CHECK_LT(NucleotideStore("AAAAC"),
		 NucleotideStore("AAAAG"));

  BOOST_CHECK_LT(NucleotideStore("AAAAC"),
		 NucleotideStore("AAAACC"));

  BOOST_CHECK_LT(NucleotideStore("AAAAC"),
		 NucleotideStore("AAAACG"));

  BOOST_CHECK_LT(NucleotideStore("AAAAC"),
		 NucleotideStore("AAAACCG"));

  BOOST_CHECK_LT(NucleotideStore("AAAACC"),
		 NucleotideStore("AAAACCG"));

  BOOST_CHECK_LT(NucleotideStore("AAAACCC"),
		 NucleotideStore("AAAACCG"));

  
  BOOST_CHECK_LT(NucleotideStore("AAAAA"),
		 NucleotideStore("TAAAA"));

  BOOST_CHECK_LT(NucleotideStore("AAAACCCCG"),
		 NucleotideStore("AAAACCCCT"));

  BOOST_CHECK_LT(NucleotideStore("AACC"),
		 NucleotideStore("AAGG"));
}

BOOST_AUTO_TEST_CASE(test_canonicalpalindrome) {
  NucleotideStore a("ACT");
  BOOST_CHECK(a.isCanonical());
  NucleotideStore a1("AC"), b1("GT");
  BOOST_CHECK(!a1.isDNAPalindrome());
  BOOST_CHECK(!a1.isDNAPalindrome());
  BOOST_CHECK(a1.getRC() == b1);

  BOOST_CHECK(!a.isDNAPalindrome());
  NucleotideStore p1("CCCGGG"), p2("CG");
  BOOST_CHECK(p1.isDNAPalindrome());
  BOOST_CHECK(p2.isDNAPalindrome());

  NucleotideStore e;
  BOOST_CHECK(!e.isDNAPalindrome());
  
}

BOOST_AUTO_TEST_CASE(test_delta) {
  using namespace std;
  NucleotideStore a("ACGTTGCA"), b("ACGTTTCA"), c;
  auto ds = a.getDelta(b);
  cout<<a<<endl<<b<<endl;
  for(const auto& d : ds) {
    cout<<d<<endl;
  }
  cout<<"---"<<endl;

  vector<NucleotideStore::Delta> expected({{(uint32_t)5, 'T', NucleotideStore::Delta::Action::Replace}});
  BOOST_CHECK(ds==expected);
  
  auto ds2= a.getDelta(c);
  for(const auto& d : ds2) {
    cout<<d<<endl;
  }
  
  NucleotideStore d("AGCCTTTCCGGA"), e("AGCCTTTCCCGGA");
  auto ds3 = d.getDelta(e);
  cout<<"---"<<endl;

  for(const auto& d : ds3) {
    cout<<d<<endl;
  }
  
  NucleotideStore f("AGCCTTTCCCGGGA");
  cout<<"---"<<endl;
  for(const auto& de : d.getDelta(f))
    cout<<de<<endl;

  
}
BOOST_AUTO_TEST_SUITE_END()
