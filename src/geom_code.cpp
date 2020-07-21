// [[Rcpp::depends(BH)]]
#include <iostream>
#include <Rcpp.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>

using namespace Rcpp;
using namespace boost::geometry;
namespace bg = boost::geometry;
using namespace std;


//##################################################################
// Function to project points on lines
//##################################################################

typedef boost::geometry::model::d2::point_xy<double> point_type;

// note : adapted from http://www.sunshine2k.de/coding/java/PointOnLine/PointOnLine.html#step5
bg::model::point<double, 2, bg::cs::cartesian> getProjectedPointOnLine(bg::model::point<double, 2, bg::cs::cartesian> p, bg::model::point<double, 2, bg::cs::cartesian> v1, bg::model::point<double, 2, bg::cs::cartesian> v2){

  // get dot product of e1, e2
  bg::model::point<double, 2, bg::cs::cartesian> e1(bg::get<0>(v2)-bg::get<0>(v1), bg::get<1>(v2)-bg::get<1>(v1));
  bg::model::point<double, 2, bg::cs::cartesian> e2(bg::get<0>(p)-bg::get<0>(v1), bg::get<1>(p) - bg::get<1>(v1));
  double valDp = bg::get<0>(e1) * bg::get<0>(e2) + bg::get<1>(e1) * bg::get<1>(e2); //dotProduct(e1, e2);
  // get squared length of e1
  double len2 = bg::get<0>(e1) * bg::get<0>(e1) + bg::get<1>(e1) * bg::get<1>(e1);
  bg::model::point<double, 2, bg::cs::cartesian> p2(bg::get<0>(v1) + (valDp * bg::get<0>(e1)) / len2,
                                                    bg::get<1>(v1) + (valDp * bg::get<1>(e1)) / len2
                                                    );
  return p2;
}


bg::model::point<double, 2, bg::cs::cartesian> ProjectedPointOnLine(bg::model::point<double, 2, bg::cs::cartesian> p, bg::model::linestring<bg::model::point<double, 2, bg::cs::cartesian>> line){
  int cnt_pt = line.size();
  bg::model::point<double, 2, bg::cs::cartesian> v1;
  bg::model::point<double, 2, bg::cs::cartesian> v2;
  double refdist = -1;
  //first : iterating over all the segments of the line
  for(int i=0; i < cnt_pt+1; ++i){
    bg::model::point<double, 2, bg::cs::cartesian> s1 = line[i];
    bg::model::point<double, 2, bg::cs::cartesian> s2 = line[i+1];
    double d1 = sqrt(pow(bg::get<0>(s1) - bg::get<0>(p),2) + pow(bg::get<1>(s1) - bg::get<1>(p),2));
    double d2 = sqrt(pow(bg::get<0>(s2) - bg::get<0>(p),2) + pow(bg::get<1>(s2) - bg::get<1>(p),2));
    double dist = d1+d2;
    if(dist < refdist || refdist==-1){
      refdist = dist;
      v1 = s1;
      v2 = s2;
    }
  }
  bg::model::point<double, 2, bg::cs::cartesian> new_point = getProjectedPointOnLine(p, v1, v2);
  return new_point;
}

// [[Rcpp::export]]
void test_func(std::string WKT, NumericVector coords){

  typedef boost::geometry::model::d2::point_xy<double> point_type;
  using segment_type = model::segment<point_type>;
  using linestring_type = model::linestring<point_type>;
  point_type point;
  boost::geometry::model::linestring<bg::model::point<double, 2, bg::cs::cartesian>> line;

  bg::model::point<double, 2, bg::cs::cartesian> p(coords[0],coords[1]);
  boost::geometry::read_wkt(WKT, line);
  //ProjectedPointOnLine(p, line);
}


//##################################################################
// Function to interpolate points on lines
//##################################################################

// [[Rcpp::export]]
DataFrame interpolate_pts(StringVector WKT, NumericVector dists){

  NumericVector Xs;
  NumericVector Ys;
  IntegerVector oids;

  typedef boost::geometry::model::d2::point_xy<double> point_type;
  using segment_type = model::segment<point_type>;
  using linestring_type = model::linestring<point_type>;
  point_type point;

  boost::geometry::model::linestring<point_type> line;

  int cnt = WKT.length();
  for(int i=0; i < cnt; ++i){
    std::string txt = Rcpp::as< std::string >(WKT(i));
    boost::geometry::read_wkt(txt, line);
    float dist = dists[i];
    line_interpolate(line, dist, point);
    double x = boost::geometry::get<0>(point);
    double y = boost::geometry::get<1>(point);
    Xs.push_back(x);
    Ys.push_back(y);
    oids.push_back(i);
  }
  DataFrame df = DataFrame::create(Named("oid") = oids ,Named("x") = Xs ,Named("y") = Ys);
  return df;
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
interpolate_pts(c("LINESTRING(0 0,2 2,3 1)",
                  "LINESTRING(0 0,2 2,3 1)",
                  "LINESTRING(0 0,2 2,3 1)"),
                c(0,2.5,3)
                )
*/
