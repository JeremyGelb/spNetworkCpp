#include <Rcpp.h>
#include <iostream>
#include<sstream>
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;
using namespace std;
// [[Rcpp::depends(RcppProgress)]]

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// #### General helper functions ####
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


void rcpp_rprintf(NumericVector v){
  // printing values of all the elements of Rcpp vector
  for(int i=0; i<v.length(); ++i){
    Rprintf("the value of v[%i] : %f \n", i, v[i]);
  }
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// #### Graph helper functions ####
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IntegerVector get_edges(int v1, IntegerVector v2s, List edge_list){
  IntegerVector Ly;
  int n = v2s.length();
  stringstream ss;
  ss << v1;
  string n1 = ss.str();
  for(int j=0; j < n; ++j){
    stringstream ss;
    int v2 = v2s[j];
    ss << v2;
    string n2 = ss.str();
    string edge_name = n1+"_"+n2;
    int edge_id = edge_list[edge_name];
    Ly.push_back(edge_id,edge_name);
  }
  return Ly;
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// #### NKDE discontinuous functions ####
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' The recursive function to calculate continuous NKDE
//'
//' @param samples_k a vector of current KDE values at sample points
//' @param graph The igraph graph object representing the graph
//' @param v the starting node
//' @param v1 the ending node of l1
//' @param l1 the considered edge during this recursion
//' @param d the current distance value (updated during recursion)
//' @param alpha the current alpha value (updated during recursion)
//' @param bw the kernel bandwidth
//' @param kernel_func the function to calculate the kernel from distance
//' @param samples a DataFrame of the samples (with spatial coordinates and belonging edge)
//' @param nodes a DataFrame representing the nodes of the graph (with spatial coordinates)
//' @param line_list a DataFrame representing the lines of the graph
//' @param depth the actual recursion depth (updated during recursion)
//' @param max_depth the maximum recursion depth (after which recursion is stopped)
//' @return a numeric vector of the kernel values for each sample
//'
NumericVector esc_kernel_rcpp(NumericVector samples_k, List edge_list, List neighbour_list, int v, int v1, int l1, double d,double alpha, double bw, Function kernel_func,  NumericVector line_weights, IntegerVector samples_edgeid, NumericVector samples_x, NumericVector samples_y, IntegerVector samples_oid, NumericVector nodes_x, NumericVector nodes_y , int depth, int max_depth){

  //step1 find the index of the right samples
  LogicalVector test = samples_edgeid==l1;
  NumericVector sampling_x = samples_x[test];
  NumericVector sampling_y = samples_y[test];
  IntegerVector sampling_oid = samples_oid[test];
  sampling_oid = sampling_oid-1;
  //extracting the X and Y coordinates of the starting node
  int v_m1 = v-1;
  double node_x = nodes_x[v_m1];
  double node_y = nodes_y[v_m1];
  //step2 calculating the distances for each sample
  NumericVector part1 = (sampling_x - node_x)*(sampling_x - node_x);
  NumericVector part2 = (sampling_y - node_y)*(sampling_y - node_y);
  NumericVector x_dists = sqrt(part1 + part2) + d;
  //step3 calculating the values of the new kernel
  NumericVector new_k = kernel_func(x_dists,bw);
  //Rcout << "    here are the new_k2 values : " << "\n";
  NumericVector new_k2 = new_k*alpha;
  //rcpp_rprintf(new_k2);
  NumericVector old_values = samples_k[sampling_oid];
  NumericVector K_final = new_k2 + old_values;
  samples_k[sampling_oid] = K_final;

  //mettre a jour d
  double d2 = line_weights[l1-1] + d;
  int new_depth = depth+1;
  if((bw>d2) & (new_depth < max_depth)){
    //on veut trouver toutes les lignes emannant de v (Lv)
    IntegerVector v_neighbours = neighbour_list[v1-1];
    IntegerVector Lv = get_edges(v1,v_neighbours,edge_list);
    int n = v_neighbours.length();
    if(n>1){
      for(int j=0; j < n; ++j){
        int li = Lv[j];
        int vi = v_neighbours[j];
        if(li==l1){
          double p2 = (n-2.0)/n;
          double n_alpha = -1 * alpha * p2;
          samples_k = esc_kernel_rcpp(samples_k, edge_list, neighbour_list, v1,vi, li, d2, n_alpha, bw,kernel_func, line_weights, samples_edgeid, samples_x, samples_y, samples_oid, nodes_x, nodes_y, new_depth, max_depth);
        }else{
          double n_alpha = alpha * (2.0/n);
          samples_k = esc_kernel_rcpp(samples_k, edge_list, neighbour_list, v1,vi, li, d2, n_alpha, bw,kernel_func, line_weights, samples_edgeid, samples_x, samples_y, samples_oid, nodes_x, nodes_y, new_depth, max_depth);
        };
      };
    };
  };
  return samples_k;
}


//' @title The main function to calculate continuous NKDE
//'
//' @param edge_list A list of edges, accessible by their names
//' @param neighbour_list a list of the neighbours of each node
//' @param events a numeric vector of the node id of each event
//' @param weights a numeric vector of the weight of each event
//' @param samples a DataFrame of the samples (with spatial coordinates and belonging edge)
//' @param bw the kernel bandwidth
//' @param kernel_func the function to calculate the kernel from distance
//' @param nodes a DataFrame representing the nodes of the graph (with spatial coordinates)
//' @param line_list a DataFrame representing the lines of the graph
//' @param max_depth the maximum recursion depth (after which recursion is stopped)
//' @return a DataFrame with two columns : the kernel values (sum_k) and the number of event reaching each samples (n)
//' @export
//'
// [[Rcpp::export]]
DataFrame continuous_nkde_cpp(List edge_list, List neighbour_list, NumericVector events, NumericVector weights, DataFrame samples, double bw, Function kernel_func, DataFrame nodes, DataFrame line_list, int max_depth, bool verbose){

  //step0 extract the columns of the dataframe
  NumericVector line_weights = line_list["weight"];
  IntegerVector samples_edgeid = samples["edge_id"];
  NumericVector samples_x = samples["X_coords"];
  NumericVector samples_y = samples["Y_coords"];
  IntegerVector samples_oid = samples["oid"];
  NumericVector nodes_x = nodes["X_coords"];
  NumericVector nodes_y = nodes["Y_coords"];

  //step 1 : mettre toutes les valeurs a 0
  NumericVector base_k = rep(0.0,samples.nrow());
  NumericVector base_count = rep(0.0,samples.nrow());

  //step2 : iterer sur chaque event
  int cnt_e = events.length()-1;
  Progress p(cnt_e, verbose);
  for(int i=0; i <= cnt_e; ++i){
    p.increment(); // update progress
    //preparer les differentes valeurs de departs pour l'event y
    int y = events[i];
    double w = weights[i];
    //on veut trouver toutes les voisins emannant de y
    IntegerVector y_neighbours = neighbour_list[y-1];
    IntegerVector Ly = get_edges(y,y_neighbours,edge_list);
    int cnt_y = Ly.length()-1;
    for(int j=0; j <= cnt_y; ++j){
      //preparing all values for this loop
      int li = Ly[j];
      int vi = y_neighbours[j];
      int d = 0 ;
      int depth = 0 ;
      double alpha = 1 ;
      //preparing the 0 values
      NumericVector samples_k = rep(0.0,samples.nrow());
      // launching recursion
      NumericVector k = esc_kernel_rcpp(samples_k, edge_list, neighbour_list ,y,vi,li,d,alpha,bw,kernel_func, line_weights, samples_edgeid, samples_x, samples_y, samples_oid, nodes_x, nodes_y, depth,max_depth);
      // getting back the sample_k values
      //NumericVector kw = (k*w);
      // getting the actual base_k values (summed at each iteration)
      NumericVector new_base_k = base_k + k*w;
      base_k = new_base_k;
      // calculating the new value
      NumericVector count = rep(0.0,k.length());
      //LogicalVector test = k>0;
      count[k>0] = 1;
      //NumericVector new_base_count = base_count + (count*w);
      base_count = base_count + (count*w);
    };

  };
  DataFrame df =  DataFrame::create( Named("sum_k") = base_k ,Named("n") = base_count);
  return df;
}


