#include "csv-parser/csv.hpp"
#include <iostream>
#include <map>
#include <fstream>
#include <optional>
#define _USESTDVECTOR_
#define _USENRERRORCLASS_
#include <fmt/core.h>
#include <fmt/os.h>
#include "nr3.h"
#include "amoeba.h"

using namespace std;

/* this reads the output of gcstats ('skfit.csv') and attempts to fit
   models on all the skews reported therein.
   
   This is an all or nothing thing, it will do all fits at once, and produce loads of csv output

   results.csv contains the summary, the _fit files are ready to plot to see how well it worked  
*/

/*
Nice fits:
  NC_000117.1: pretty ok, could be better
  NC_006347.1: pretty ok, could be better

Biology not cooperating:
  NC_002663.1: You could fit a more advanced model to this

Bad fits:
  NC_013132.1: ends with a GC excess
  NC_002967.9: ends with GC excess


Bad fits that should just work better:
  NZ_CP007747.1
  NC_014119.1

  NZ_LR890466.1

REVERSE GC/TA:
  NC_015949.1
*/

void doDump(const auto& g)
{
  auto out = fmt::output_file(g.first+"_fit.csv");

  vector<string> hdrs({"pos", "gcskew", "predgcskew",
			    "taskew", "predtaskew",
			    "sbskew", "predsbskew",
			    "gc0skew", "predgc0skew",
			    "gc1skew", "predgc1skew",
			    "gc2skew", "predgc2skew",
			    "ta0skew", "predta0skew",
			    "ta1skew", "predta1skew",
			    "ta2skew", "predta2skew",
                            
			    "gcNGskew", "predgcNGskew",
                            "taNGskew", "predtaNGskew",
                            "predleading"
    });
  bool firsthdr=true;
  for(const auto& h : hdrs) {
    if(!firsthdr)
      out.print(",");
    else
      firsthdr=false;
    out.print("{}", h);
  }
  out.print("\n");
  for(auto& skp : g.second) {
    out.print("{}", skp.pos);
    auto vals= make_tuple(
			 skp.gc.skew, skp.gc.predskew, 
			 skp.ta.skew, skp.ta.predskew,
			 skp.sb.skew, skp.sb.predskew,
			 skp.gc0.skew, skp.gc0.predskew,
			 skp.gc1.skew, skp.gc1.predskew,
			 skp.gc2.skew, skp.gc2.predskew,
			 skp.ta0.skew, skp.ta0.predskew,
			 skp.ta1.skew, skp.ta1.predskew,
			 skp.ta2.skew, skp.ta2.predskew,
 
			 skp.gcNG.skew, skp.gcNG.predskew,
                         skp.taNG.skew, skp.taNG.predskew,
                         skp.predleading
                          );
    std::apply([&out](auto&&... args) {((out.print(",{}", args)), ...);}, vals);
    out.print("\n");
               
  }
}

struct BiasStats
{
  int minpos;
  int maxpos;
  int shift;
  double div;
  double alpha1{0};
  double alpha2{0};
  double origRMS{0};
  double rms{0};
};

struct SKPos
{
  int pos;
  int gccount, ngcount;
  int acounts2, ccounts2, gcounts2, tcounts2;

  struct SkewDeets
  {
    double skew;
    double skewdiff{0};
    double predskew{0};
    double predskewdiff{0};
  } gc, ta, sb,
    gc0, gc1, gc2,
    ta0, ta1, ta2,
    gcNG, taNG;
  int predleading;
};

/* 
   Two fitting models, one that fits two steepness parameters and two inflection points.

   Another fitting model is simpler: a chromosome shift parameter and a single steepness parameter

   To evaluate the second model:

    Shift parameter S
    Steepness parameter alpha

    The relative position L is locus/len
    The prediction is alpha*L up to L=0.5 and from there on it is 
    0.5*alpha - (L-0.5)*alpha =alpha*(0.5 -L + 0.5) = alpha*(1-L)

    To compare, iterate over L and take the chromosome value at L-shift.
    In addition, shift down the chromosome by cromo[shift]
    If L-shift is negative, add the cromosome length
*/



BiasStats doAnalysis(std::function<SKPos::SkewDeets&(SKPos&)> getter, vector<SKPos>& chromo, std::optional<double> leshift = std::optional<double>(), std::optional<double> lediv = std::optional<double>())
{
    auto getAbsAvg = [&getter](vector<SKPos>& v) {
		       double tot=0.0;
		       unsigned int count=1;
		       for(auto& sk : v) {
			 tot += fabs(getter(sk).skew);
			 count++;
		       }
		       return tot/count;
		     };


    // this is called by the optimizer to evaluate a fit based on the parameters
    auto func2 = [&chromo, &getter](const vector<double>& params) {
		  // alpha, shift
		  auto& alpha1 = params[0];
                  auto& alpha2 = params[1];
		  auto shift = params[2];
                  auto div = params[3];
                  //                  cout<<"func2 alpha1 "<<alpha1<<" alpha2 "<<alpha2<<" shift "<<shift<<" div "<<div<<endl;
		  static int counter=0;
		  double missRMS=0;
		  double chromolen = chromo.rbegin()->pos;
                  double cumul=0;
                  int prevpos = 0;
		  for(auto& skp : chromo) {
		    int rpos = skp.pos;
                    if((rpos < shift || rpos > shift + div*chromolen) && (shift > 0 || rpos < chromolen + shift)  ) { // going down
                      getter(skp).predskewdiff = -alpha1;
                      cumul -= (skp.pos - prevpos)*alpha1;
                      getter(skp).predskew=cumul;
                      skp.predleading=false;
		    }
		    else {
                      getter(skp).predskewdiff = alpha2;
                      cumul += (skp.pos - prevpos)*alpha2;
                      getter(skp).predskew=cumul;                      
                      skp.predleading=true;
		    }
                    prevpos = skp.pos;
		    //		    fits<<counter<<";"<<shift<<";"<<alpha<<";"<<skp.pos<<";"<<rpos<<";"<<getter(skp).skew<<";"<<getter(skp).predskew<<endl;
		    
                    //		    missRMS += pow(getter(skp).skewdiff - getter(skp).predskewdiff, 2);
                    missRMS += pow(getter(skp).skew - getter(skp).predskew, 2);
                  }
		  ++counter;
                  //                  cout<<"Returning missRMS "<<missRMS<<endl;
		  return missRMS;
		};

    
    
    Amoeba amo(0.00001); // was 0.0001
    //    vector<double> params({alpha, alpha, 1.0*maxpos, 1.0*minpos});


    // create derivative
    for(auto iter = chromo.begin() ; iter != chromo.end() ; ++iter) {
      if(iter==chromo.begin()) {
        getter(*iter).skewdiff=0;
        continue;
      }

      getter(*iter).skewdiff = getter(*iter).skew -  getter(*std::prev(iter)).skew;
    }

    // initial guess at alphas
    double totalpha=0;
    auto citer = chromo.begin();
    for(; citer < chromo.begin() + chromo.size()/2; ++citer) {
      totalpha += getter(*citer).skewdiff;
    }
    for(; citer != chromo.end(); ++citer) {
      totalpha -= getter(*citer).skewdiff;
    }

    
    //    cout<<"minval: "<<minval<<", maxval: "<<maxval<<", alpha ";
    double alpha2 = totalpha/chromo.rbegin()->pos;
    cout<<"Estimated alpha "<<alpha2<<endl;
    
    vector<double> params({alpha2, alpha2});
    vector<double> deltas({10,10});
    
    if(leshift) {
      params.push_back(*leshift);
      deltas.push_back(0); // the shift is fixed
      params.push_back(*lediv); // excess, for uneven lagging and leading strands
      deltas.push_back(0); // fixed

    }
    else {
      params.push_back(0);
      deltas.push_back(0.1*chromo.rbegin()->pos);

      params.push_back(0.5); // the divider
      deltas.push_back(0.01); // 
    }
    //    cout<<"Initial: alpha "<<alpha<<" maxpos "<< maxpos<<" minpos "<<minpos<<" FuncRMS: "<<sqrt(func2(params))/chromo.size()/getAbsAvg(chromo)<<endl;

    // first attempt, based on an educated guess
    double newRMS;

    cout<<"Optimizing: ";
    auto np = amo.minimize(params, deltas, func2);
    
    newRMS = func2(np);
    //    auto nv = sqrt(func2(np));
    //    cout << "alpha: " << np[0]<< ", shift: "<<np[1]<<", newRMS: "<<newRMS<<", sqrt(rms): ";
    //    cout << nv <<", chromo.size(): "<<chromo.size()<<", absavg: "<<getAbsAvg(chromo)<<endl;

    params=np;

    // seeing if we can improve on this with some brute force
    map<double, vector<double>> scores;
    scores[newRMS]=params; // make sure the original fit is still in there

    for(int x=1; x < 5 ; ++x) {
      cout<<"Optimizing round "<<x<<" ";
      auto np = amo.minimize(params, deltas, func2);
      newRMS = func2(np);
      cout << "alpha1: " << np[0]<< ", alpha2: "<<np[1] <<", newRMS: "<<newRMS<<endl;
      params=np;
      scores[newRMS]=params;
    }

    auto& best = *scores.begin();
    cout<<"Best result: "<<best.first<<" ";
    func2(best.second);  // this fills out the struct for the graphs
    np = best.second;
    //    cout << "alpha: " << np[0]<< ", shift: "<<np[1]<<", new RMS: "<<best.first<<endl;
    newRMS = best.first;
    params=np;

    BiasStats ret;
    ret.alpha1= params[0];
    ret.alpha2= params[1];
    ret.shift = params[2];
    ret.div = params[3];

    // the chromo.size() is CORRECT here - it represents the number of samples, not the size of the chromosome!
    // the newRMS is built up out of chromo.size() measurements!
    ret.rms = sqrt(newRMS/chromo.size())  / getAbsAvg(chromo);
    cout<<"alpha1 "<<ret.alpha1<<" alpha2 "<<ret.alpha2 <<". Done, adjusted RMS: " << ret.rms <<", raw "<<newRMS  << endl;
    // XXX
    return ret;
}


int main(int argc, char **argv)
{

  string fname("skplot.csv");
  if(argc > 1)
    fname=argv[1];
  cerr<<"Reading 'skplot.csv' file from "<<fname<<endl;
  csv::CSVReader reader(fname);

  set<string> filter;
  for(int n = 2; n < argc; ++n) {
    cerr<<"Filtering on "<<argv[n]<<endl;
    filter.insert(argv[n]);
  }
  
  /*
  for(const auto& s: reader.get_col_names()) {
    cout<<s<<": "<<reader.index_of(s)<<endl;
  }
  */

  auto safeindex = [&reader](const char* name) {
		     int ret = reader.index_of(name);
		     if(ret < 0)
		       throw runtime_error("Could not find index for "+string(name));
		     return ret;
		   };
    
  
  int name = safeindex("name");
  int abspos = safeindex("abspos"); //
  int gcskew = safeindex("gcskew");
  int taskew = safeindex("taskew");
  int sbskew = safeindex("pospos"); // 
  int gc0 = safeindex("gcskew0");
  int gc1 = safeindex("gcskew1");
  int gc2 = safeindex("gcskew2");
  int ta0 = safeindex("taskew0");
  int ta1 = safeindex("taskew1");
  int ta2 = safeindex("taskew2");

  int acounts2 = safeindex("acounts2");
  int ccounts2 = safeindex("ccounts2");
  int gcounts2 = safeindex("gcounts2");
  int tcounts2 = safeindex("tcounts2"); 
  
  int gcNG = safeindex("gcskewNG");
  int taNG = safeindex("taskewNG");
  int gccount = safeindex("gccount");
  int ngcount = safeindex("ngcount");

  struct {
    BiasStats gc, ta, sb,
      gc0, gc1, gc2, gcng,
      ta0, ta1, ta2, tang;
  } biases;

  map<string, vector<SKPos>> chromosomes;
  cout<<"Reading "<<fname<<endl;
  for (csv::CSVRow& row: reader) { // Input iterator
    if(!filter.empty() && !filter.count(row[name].get<string>()))
      continue;

    chromosomes[row[name].get<string>()].push_back({row[abspos].get<int>(),
					      row[gccount].get<int>(),
					      row[ngcount].get<int>(),                                                    
					      row[acounts2].get<int>(),
					      row[ccounts2].get<int>(),
					      row[gcounts2].get<int>(),
					      row[tcounts2].get<int>(),
					      {row[gcskew].get<double>()},
					      {row[taskew].get<double>()},
					      {row[sbskew].get<double>()},
					      {row[gc0].get<double>()},
					      {row[gc1].get<double>()},
					      {row[gc2].get<double>()},

    					      {row[ta0].get<double>()},
					      {row[ta1].get<double>()},
					      {row[ta2].get<double>()},

					      {row[gcNG].get<double>()},
 					      {row[taNG].get<double>()}
			


      });
  }

  cout<<"Got "<<chromosomes.size()<<" chromosomes"<<endl;
  auto resos = fmt::output_file("results.csv");

  vector<string> headers({"name", "flipped", "siz", "gccount", "ngcount", "acounts2", "ccounts2", "gcounts2", "tcounts2", "maxpos",
			       "alpha1gc", "minpos", "alpha2gc", "shift", "div",
			       "alpha1ta", "alpha2ta",
			       "alpha1sb", "alpha2sb",
			       
			       "alpha1gc0", "alpha2gc0",
			       "alpha1gc1", "alpha2gc1",
			       "alpha1gc2", "alpha2gc2",
			       
			       "alpha1ta0", "alpha2ta0",
			       "alpha1ta1", "alpha2ta1",
			       "alpha1ta2", "alpha2ta2",
			       
			       "alpha1gcNG", "alpha2gcNG",
			       "alpha1taNG", "alpha2taNG",
			       
			       "rmsGC", "rmsTA", "rmsSB",
			       "rmsGC0", "rmsGC1", "rmsGC2",
			       "rmsTA0", "rmsTA1", "rmsTA2",
                          "rmsGCNG", "rmsTANG"});
  bool firsthdr=true;
  for(const auto& h: headers) {
    if(!firsthdr)
      resos.print(",{}", h);
    else {
      resos.print("{}", h);
      firsthdr=false;
    }
  }
  resos.print("\n");


  int counter=0;
  for(auto& g : chromosomes) {
    ++counter;
    if(!filter.empty() && !filter.count(g.first))
      continue;
    cout<<"Start with chromosome "<<g.first<<", "<<counter<<"/"<<chromosomes.size()<<endl;
    //    double minskew = 1e10, maxskew = -1e10;
    int minpos = -1, maxpos = -1;


    // de-trend gc0, gc1, gc2, gcNG
    double detrend0 = g.second.rbegin()->gc0.skew / g.second.size();
    double detrend1 = g.second.rbegin()->gc1.skew / g.second.size();
    double detrend2 = g.second.rbegin()->gc2.skew / g.second.size();
    double detrendNG = g.second.rbegin()->gcNG.skew / g.second.size();
    int pos = 0;
    for(auto& skp : g.second) {
      skp.gc0.skew -= detrend0 * pos;
      skp.gc1.skew -= detrend1 * pos;
      skp.gc2.skew -= detrend2 * pos;
      skp.gcNG.skew -= detrendNG * pos;
      pos++;
    }

    // de-trend ta0, ta1, ta2, taNG
    detrend0 = g.second.rbegin()->ta0.skew / g.second.size();
    detrend1 = g.second.rbegin()->ta1.skew / g.second.size();
    detrend2 = g.second.rbegin()->ta2.skew / g.second.size();
    detrendNG = g.second.rbegin()->taNG.skew / g.second.size();
    pos = 0;
    for(auto& skp : g.second) {
      skp.ta0.skew -= detrend0 * pos;
      skp.ta1.skew -= detrend1 * pos;
      skp.ta2.skew -= detrend2 * pos;
      skp.taNG.skew -= detrendNG * pos;
      pos++;
    }

    
    bool flipped{false};

    // this is an adaptor for the evaluator below
    auto gcgetter = [](auto& skp) -> SKPos::SkewDeets& {
		      return skp.gc;
		  };

    //    double alpha = (maxskew - minskew) / (minpos - maxpos);

    cout << "GC skew"<<endl;
    biases.gc = doAnalysis(gcgetter, g.second);
    cout<<"Shift: "<<biases.gc.shift<<", gc.alpha1: "<<biases.gc.alpha1<<", div: "<<biases.gc.div << endl;

    if(biases.gc.alpha1 < 0) {
      cout<<"Flipping "<<g.first<<endl;
      flipped=true;
      for(auto& skp : g.second) {
	skp.gc.skew   *= -1;
	skp.gc0.skew  *= -1;
	skp.gc1.skew  *= -1;
	skp.gc2.skew  *= -1;
	skp.gcNG.skew *= -1;

	skp.ta.skew   *= -1;
	skp.ta0.skew  *= -1;
	skp.ta1.skew  *= -1;
	skp.ta2.skew  *= -1;

	skp.taNG.skew *= -1;
        skp.sb.skew *= -1; // perhaps not?
      }

      biases.gc = doAnalysis(gcgetter, g.second);
      cout<<"Shift: "<<biases.gc.shift<<endl;
    }
    
    /// now TA 
    auto tagetter = [](SKPos& skp) -> SKPos::SkewDeets& {
	       return skp.ta;
	     };
    cout << "TA skew"<<endl;
    biases.ta = doAnalysis(tagetter, g.second, biases.gc.shift, biases.gc.div);

    // SB
    auto sbgetter = [](SKPos& skp) -> SKPos::SkewDeets& {
	       return skp.sb;
	     };

    cout << "SB skew"<<endl;
    biases.sb = doAnalysis(sbgetter, g.second, biases.gc.shift, biases.gc.div);

    // GC0
    auto gc0getter = [](SKPos& skp) -> SKPos::SkewDeets& {
	       return skp.gc0;
	     };
    cout << "gc0 skew"<<endl;
    biases.gc0 = doAnalysis(gc0getter, g.second, biases.gc.shift, biases.gc.div);

    // GC1
    auto gc1getter = [](SKPos& skp) -> SKPos::SkewDeets& {
	       return skp.gc1;
	     };
    //    cout << "gc1 skew"<<endl;
    biases.gc1 = doAnalysis(gc1getter, g.second, biases.gc.shift, biases.gc.div);

    // GC2
    auto gc2getter = [](SKPos& skp) -> SKPos::SkewDeets& {
	       return skp.gc2;
	     };

    //    cout << "gc2 skew"<<endl;
    biases.gc2 = doAnalysis(gc2getter, g.second, biases.gc.shift, biases.gc.div);

    // gcng
    auto gcnggetter = [](SKPos& skp) -> SKPos::SkewDeets& {
	       return skp.gcNG;
	     };
    //    cout << "ng skew"<<endl;
    biases.gcng = doAnalysis(gcnggetter, g.second, biases.gc.shift, biases.gc.div);

    // TA0
    auto ta0getter = [](SKPos& skp) -> SKPos::SkewDeets& {
	       return skp.ta0;
	     };
    //    cout << "ta0 skew"<<endl;
    biases.ta0 = doAnalysis(ta0getter, g.second, biases.gc.shift, biases.gc.div);

    // TA1
    auto ta1getter = [](SKPos& skp) -> SKPos::SkewDeets& {
	       return skp.ta1;
	     };
    //    cout << "ta1 skew"<<endl;
    biases.ta1 = doAnalysis(ta1getter, g.second, biases.gc.shift, biases.gc.div);

    // TA2
    auto ta2getter = [](SKPos& skp) -> SKPos::SkewDeets& {
	       return skp.ta2;
	     };

    //    cout << "ta2 skew"<<endl;
    biases.ta2 = doAnalysis(ta2getter, g.second, biases.gc.shift, biases.gc.div);

    // tang
    auto nggetter = [](SKPos& skp) -> SKPos::SkewDeets& {
	       return skp.taNG;
	     };
    //    cout << "ng skew"<<endl;
    biases.tang = doAnalysis(nggetter, g.second, biases.gc.shift, biases.gc.div);
    
    doDump(g);
    cout<<"Writing out alpha1ta: "<<biases.ta.alpha1<<endl;
    resos.print("{}", g.first);
    auto vals=make_tuple((int)flipped, g.second.rbegin()->pos, g.second.rbegin()->gccount, g.second.rbegin()->ngcount,
			    g.second.rbegin()->acounts2,
			    g.second.rbegin()->ccounts2,
			    g.second.rbegin()->gcounts2,
			    g.second.rbegin()->tcounts2, 
			    maxpos, biases.gc.alpha1, minpos, biases.gc.alpha2, biases.gc.shift, biases.gc.div,
			    biases.ta.alpha1, biases.ta.alpha2,
			    biases.sb.alpha1, biases.sb.alpha2,
			    
			    biases.gc0.alpha1, biases.gc0.alpha2,
			    biases.gc1.alpha1, biases.gc1.alpha2,
			    biases.gc2.alpha1, biases.gc2.alpha2,

			    biases.ta0.alpha1, biases.ta0.alpha2,
			    biases.ta1.alpha1, biases.ta1.alpha2,
			    biases.ta2.alpha1, biases.ta2.alpha2,

			    biases.gcng.alpha1, biases.gcng.alpha2,
			    biases.tang.alpha1, biases.tang.alpha2,
			    
			    biases.gc.rms, biases.ta.rms, biases.sb.rms,
			    biases.gc0.rms, biases.gc1.rms, biases.gc2.rms,
			    biases.ta0.rms, biases.ta1.rms, biases.ta2.rms,
			    biases.gcng.rms, biases.tang.rms);
    std::apply([&resos](auto&&... args) {((resos.print(",{}", args)), ...);}, vals);
    resos.print("\n");
  }
}





#if 0    
    // this is called by the optimizer to evaluate a fit based on the parameters
    auto func = [&chromo, &getter](const vector<double>& params) {
		  // alpha1, alpha2, maxpos, minpos
		  auto& alpha1 = params[0];
		  auto& alpha2 = params[1];
		  auto& maxpos = params[2];
		  auto& minpos = params[3];
		  double missRMS=0;
		  double maxskew = alpha1 * maxpos;
		  double minskew = maxskew - (minpos - maxpos) * alpha2;
			  
		  for(auto& skp : chromo) {
		    if(skp.pos < maxpos) { // straight up until maxpos
		      getter(skp).predskew = alpha1 * skp.pos;
		    }
		    else if(skp.pos < minpos) {  // then we go down again
		      getter(skp).predskew = maxskew - alpha2 * ( skp.pos - maxpos);
		    }
		    else { // skp.pos > minpos, and up again
		      getter(skp).predskew = minskew + alpha1 * (skp.pos - minpos);
		    }
		    missRMS += pow(getter(skp).skew - getter(skp).predskew, 2);
		  }
		  return missRMS;
		};
#endif



    /*
    // pick some reasonable defaults for the optimizer to start its work with
    for(const auto& skp : g.second) {
      if(skp.gc.skew > maxskew) {
	maxskew = skp.gc.skew;
	maxpos = skp.pos;
      }
    }

    for(const auto& skp : g.second) {
      if(skp.pos > maxpos && skp.gc.skew < minskew) {
	minskew = skp.gc.skew;
	minpos = skp.pos;
      }
    }

    if(minpos < 0) { // pick a global mininum then
      for(const auto& skp : g.second) {
	if(skp.gc.skew < minskew) {
	  minskew = skp.gc.skew;
	  minpos = skp.pos;
	}
      }
    }
*/
