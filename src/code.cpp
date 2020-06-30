#include <Rcpp.h>
#include <iostream>
#include<sstream>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <unordered_map>
#include <queue>
using namespace Rcpp;
using namespace std;
// [[Rcpp::depends(RcppProgress)]]
//[[Rcpp::plugins(cpp11)]]

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// #### General helper functions ####
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


void rcpp_rprintf(NumericVector v){
  // printing values of all the elements of Rcpp vector
  for(int i=0; i<v.length(); ++i){
    Rprintf("the value of v[%i] : %f \n", i, v[i]);
  }
}

void rcpp_rprintf_int(IntegerVector v){
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
// #### NKDE continuous functions ####
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
  NumericVector new_k2 = new_k*alpha;
  NumericVector old_values = samples_k[sampling_oid];
  NumericVector K_final = new_k2 + old_values;
  samples_k[sampling_oid] = K_final;

  //mettre a jour d
  double d2 = line_weights[l1-1] + d;
  if((bw>d2) & (depth < max_depth)){
    //on veut trouver toutes les lignes emannant de v (Lv)
    IntegerVector v_neighbours = neighbour_list[v1-1];
    IntegerVector Lv = get_edges(v1,v_neighbours,edge_list);
    int n = v_neighbours.length();
    int new_depth;
    //updating depth
    if(n>2){
      new_depth = depth+1;
    }else{
      new_depth = depth+0;
    }
    if(n>1){
      for(int j=0; j < n; ++j){
        int li = Lv[j];
        int vi = v_neighbours[j];
        if(li==l1){
          double p2 = (n-2.0)/n;
          double n_alpha = -1.0 * alpha * p2;
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
//' @param verbose a boolean indicating if the function must print its progress
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


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// #### NKDE discontinuous functions ####
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NumericVector titi_func(NumericVector d, double bw){
  NumericVector u = d/bw ;
  NumericVector k = ((15/16)*pow((1-pow(u,2)),2)) / bw;
  LogicalVector test = d>bw;
  k[test] = 0;
  return k;
}


// arma::vec titi_func(arma::vec d, double bw){
//   arma::vec u = d/bw ;
//   arma::vec k = ((15/16)*arma::pow((1-arma::pow(u,2)),2)) / bw;
//   arma::uvec ids = find(d>bw);
//   k.elem(ids).fill(0.0);
//   return k;
// }



//I have to rewrite it in an efficient way !

//this dict work as expected !
unordered_map<int,unordered_map<int,int>> make_map(DataFrame df){
  unordered_map<int,unordered_map<int,int>> mymap;
  IntegerVector starts = df["start_oid"];
  IntegerVector ends = df["end_oid"];
  IntegerVector edge_id = df["graph_id"];
  int cnt = starts.length();
  for(int i=0; i < cnt; ++i){
    mymap[starts[i]][ends[i]]=edge_id[i];
    mymap[ends[i]][starts[i]]=edge_id[i];
  }
  return mymap;
}


void rcpp_rcout_num(NumericVector v){
  // printing value of vector
  Rcout << "" << v << "\n";
}

void rcpp_rcout_int(IntegerVector v){
  // printing value of vector
  Rcout << "The value of v : " << v << "\n";
}


// this is the recursive function
NumericVector esd_kernel_rcpp(NumericVector samples_k, unordered_map<int,unordered_map<int,int>> edge_dict,
                              List neighbour_list ,int v,int prev_node,double d, double alpha, double bw,
                              Function kernel_func, NumericVector line_weights, IntegerVector samples_edgeid,
                              NumericVector samples_x, NumericVector samples_y, IntegerVector samples_oid,
                              NumericVector nodes_x, NumericVector nodes_y, int depth, int max_depth){
  //step1 : find all the neighbours
  IntegerVector neighbours = neighbour_list[v-1];

  //step2 : iterate over the neighbours
  int cnt_n = neighbours.length();
  int new_depth;
  if(cnt_n>2){
    new_depth = depth+1;
  }else{
    new_depth = depth;
  }

  double new_alpha = alpha * (1.0/(cnt_n-1.0));
  //if we have only one neighbour, we must stop
  if(cnt_n>1 or prev_node<=0){
    for(int i=0; i < cnt_n; ++i){
      int v2 = neighbours[i];
      //on ne veut pas revenir en arriere !
      if(v2!=prev_node){
        //find the edge between the two nodes
        int edge_id = edge_dict[v][v2];
        //find the samples on that edge
        LogicalVector test = samples_edgeid == edge_id;
        NumericVector sampling_x = samples_x[test];
        NumericVector sampling_y = samples_y[test];

        //extracting the X and Y coordinates of the starting node
        int v_m1 = v-1;
        double node_x = nodes_x[v_m1];
        double node_y = nodes_y[v_m1];

        //calculating the distances
        NumericVector part1 = (sampling_x - node_x)*(sampling_x - node_x) + (sampling_y - node_y)*(sampling_y - node_y);
        NumericVector x_dists = sqrt(part1) + d;

        //step3 calculating the values of the new kernel
        NumericVector new_k = kernel_func(x_dists,bw);
        new_k = new_k*new_alpha;
        NumericVector old_values = samples_k[test];
        new_k = new_k + old_values;
        samples_k[test] = new_k;

        //evaluating for the next move
        double d2 = line_weights[edge_id-1] + d;

        if (d2<bw and new_depth<max_depth){
          //si les conditions sont remplies, on repart pour une iteration
          samples_k = esd_kernel_rcpp(samples_k, edge_dict, neighbour_list ,v2,v,d2,new_alpha,bw,kernel_func, line_weights, samples_edgeid, samples_x, samples_y, samples_oid, nodes_x, nodes_y, new_depth,max_depth);
        }
      }
    }
  }

  return samples_k;

}

// this is the inline function
NumericVector esd_kernel_rcpp2(NumericVector samples_k, unordered_map<int,unordered_map<int,int>> edge_dict,
                              List neighbour_list ,int v, double bw,
                              Function kernel_func, NumericVector line_weights, IntegerVector samples_edgeid,
                              NumericVector samples_x, NumericVector samples_y, IntegerVector samples_oid,
                              NumericVector nodes_x, NumericVector nodes_y, int depth, int max_depth){
  //step0 : generate the queue
  queue <List> data_holder;
  //step1 : generate the first case
  List cas1 = List::create(Named("d")=0.0,
                           Named("alpha")=1.0,
                           Named("v") = v,
                           Named("prev_node") = -999,
                           Named("depth") = 0
                           );


  data_holder.push(cas1);

  //lancement des iterations
  while(data_holder.empty()==FALSE){
    //unpacking
    List cas = data_holder.front();
    data_holder.pop();
    int v = cas["v"];
    int depth = cas["depth"];
    double d = cas["d"];
    double alpha = cas["alpha"];
    int prev_node = cas["prev_node"];

    //step1 : find all the neighbours
    IntegerVector neighbours = neighbour_list[v-1];

    //step2 : iterate over the neighbours
    int cnt_n = neighbours.length();
    int new_depth;
    if(cnt_n>2){
      new_depth = depth+1;
    }else{
      new_depth = depth;
    }

    double new_alpha = alpha * (1.0/(cnt_n-1.0));

    //if we have only one neighbour, we must stop
    if(cnt_n>1 or prev_node<=0){
      for(int i=0; i < cnt_n; ++i){
        int v2 = neighbours[i];
        //on ne veut pas revenir en arriere !
        if(v2!=prev_node){
          //find the edge between the two nodes
          int edge_id = edge_dict[v][v2];
          //find the samples on that edge
          LogicalVector test = samples_edgeid == edge_id;
          NumericVector sampling_x = samples_x[test];
          NumericVector sampling_y = samples_y[test];

          //extracting the X and Y coordinates of the starting node
          int v_m1 = v-1;
          double node_x = nodes_x[v_m1];
          double node_y = nodes_y[v_m1];

          //calculating the distances
          NumericVector part1 = (sampling_x - node_x)*(sampling_x - node_x) + (sampling_y - node_y)*(sampling_y - node_y);
          NumericVector x_dists = sqrt(part1) + d;

          //step3 calculating the values of the new kernel
          NumericVector new_k = kernel_func(x_dists,bw);
          new_k = new_k*new_alpha;
          NumericVector old_values = samples_k[test];
          new_k = new_k + old_values;
          samples_k[test] = new_k;

          //evaluating for the next move
          double d2 = line_weights[edge_id-1] + d;

          if (d2<bw and new_depth<max_depth){
            List new_cas = List::create(Named("d")=d2,
                                        Named("alpha")=new_alpha,
                                        Named("v") = v2,
                                        Named("prev_node") = v,
                                        Named("depth") = new_depth
            );
            data_holder.push(new_cas);
          }
        }
      }
    }
  }

  return samples_k;

}




//' @title The main function to calculate discontinuous NKDE
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
//' @param verbose a boolean indicating if the function must print its progress
//' @return a DataFrame with two columns : the kernel values (sum_k) and the number of event reaching each samples (n)
//' @export
//'
// [[Rcpp::export]]
DataFrame discontinuous_nkde_cpp(List edge_list, List neighbour_list, NumericVector events, NumericVector weights, DataFrame samples, double bw, Function kernel_func, DataFrame nodes, DataFrame line_list, int max_depth, bool verbose){

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

  //calculer le dictionnaire des lignes
  unordered_map<int,unordered_map<int,int>> edge_dict = make_map(line_list);

  //step2 : iterer sur chaque event
  int cnt_e = events.length()-1;
  Progress p(cnt_e, verbose);
  for(int i=0; i <= cnt_e; ++i){
    p.increment(); // update progress
    //preparer les differentes valeurs de departs pour l'event y
    int y = events[i];
    double w = weights[i];
    //preparing all values for this loop
    int prev_node = -999;
    int d = 0 ;
    int depth = 0 ;
    double alpha = 1 ;
    NumericVector samples_k = rep(0.0,samples.nrow());
    // launching recursion
    NumericVector k = esd_kernel_rcpp(samples_k, edge_dict, neighbour_list ,y,prev_node,d,alpha,bw,kernel_func, line_weights, samples_edgeid, samples_x, samples_y, samples_oid, nodes_x, nodes_y, depth,max_depth);
    //NumericVector k = esd_kernel_rcpp2(samples_k, edge_dict, neighbour_list ,y,bw,kernel_func, line_weights, samples_edgeid, samples_x, samples_y, samples_oid, nodes_x, nodes_y, depth,max_depth);
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
  DataFrame df =  DataFrame::create( Named("sum_k") = base_k ,Named("n") = base_count);
  return df;
}


