#include <string>
#include <vector>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

using namespace std;

string greet(const std::string& in)
{
  return in +", " + in;
}

vector<string> more()
{
  return vector<string>{"boeh", "bah", "beh"};
}

boost::python::list moredoub()
{
  boost::python::list ret;
  for(double n = 0 ; n < 1000000; ++n)
    ret.append(n);
  
  return ret;
}

BOOST_PYTHON_MODULE(libbridge)
{
    using namespace boost::python;
    def("greet", greet);
    
    class_<std::vector<string> >("stl_vector_string")
      .def(vector_indexing_suite<std::vector<string> >())
      ;

    class_<std::vector<double> >("stl_vector_double")
      .def(vector_indexing_suite<std::vector<double> >())
      ;

    
    def("more", more);
    def("moredoub", moredoub);


}
