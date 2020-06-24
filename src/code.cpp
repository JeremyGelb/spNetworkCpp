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


NumericVector esd_kernel_rcpp(int y, List edge_list, List neighbour_list, double bw, Function kernel_func,  NumericVector line_weights, IntegerVector samples_edgeid, NumericVector samples_x, NumericVector samples_y, IntegerVector samples_oid, NumericVector nodes_x, NumericVector nodes_y){
  //definir les premiere valeurs a 0
  NumericVector samples_k = rep(0.0,samples_x.length());
  //definir la premiere serie de parametres
  List all_parameters = List::create(
      List::create(Named("v")=y,
          Named("d")=0.0,
          Named("alpha") = 1.0,
          Named("prev_node")=-999)
  );

  //lancement des iterations
  while(all_parameters.length()>0){
    //on cree une petite liste vide
    List new_parameters = List::create();
    int cnt_list = all_parameters.length();
    //on itere sur les cas en cours
    for(int i=0; i < cnt_list; ++i){
      List params = all_parameters[i];

      //step1 : unpacking the values
      int v = params["v"];
      double alpha = params["alpha"];
      int prev_node = params["prev_node"];
      double d = params["d"];
      //step2 : on trouve les voisins de v (et on enleve le precedent node)
      IntegerVector v_neighbours_all = neighbour_list[(v-1)];
      LogicalVector test = v_neighbours_all != prev_node;
      IntegerVector v_neighbours = v_neighbours_all[test];
      //avec ces voisins, on peut setter le new_alpha
      double new_alpha;
      if(prev_node>=0){
        int n_nei = v_neighbours.length();
        new_alpha = (1.0/((double)n_nei)) * alpha;
      }else{
        new_alpha = 1.0;
      }
      if(v_neighbours.length()>0){
        //step3 on trouve les edges entre v et ses voisins
        IntegerVector edges = get_edges(v,v_neighbours,edge_list);

        double node_x = nodes_x[(v-1)];
        double node_y = nodes_y[(v-1)];

        //step4 :  on va iterer sur chacune de ces lignes
        int cnt_j = edges.length();
        for(int j=0; j < cnt_j; ++j){
          int li = edges[j];

          int vi = v_neighbours[j];

          //il faut trouver les echantillons concernes
          LogicalVector test = samples_edgeid == li;
          NumericVector sub_samples_x = samples_x[test];
          NumericVector sub_samples_y = samples_y[test];

          // il faut maintenant calculer les distances entre ces observations
          NumericVector part1 = pow((node_x - sub_samples_x),2);
          NumericVector part2 = pow((node_y - sub_samples_y),2);
          NumericVector d1 = sqrt(part1+part2);
          NumericVector d2 = d1 + d;

          //on calcule maintenant la valeur kernel
          NumericVector k = kernel_func(d2,bw);
          NumericVector k1 = k * new_alpha;
          //et on l'ajoute a la valeur precedente
          NumericVector old_k = samples_k[test];
          NumericVector new_k = old_k + k1;
          samples_k[test] = new_k;

          //il ne reste plus que a voir si on peut continuer sur le prochain noeud

          double d3 = d + line_weights[(li-1)];

          if(d3<bw){

            List new_params = List::create(
                 Named("v") = vi,
                 Named("prev_node") = v,
                 Named("d") = d3,
                 Named("alpha") = new_alpha
                 );

            new_parameters.push_back(new_params);

          }

        }
      }

    }
    //on reset les nouveaux parameters
    for(int j=0; j < all_parameters.length(); ++j){
      all_parameters.erase(j);
    }
    for(int j=0; j < new_parameters.length(); ++j){
      all_parameters.push_back(new_parameters[j]);
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
//' @param verbose a boolean indicating if the function must print its progress
//' @return a DataFrame with two columns : the kernel values (sum_k) and the number of event reaching each samples (n)
//' @export
//'
// [[Rcpp::export]]
DataFrame discontinuous_nkde_cpp(List edge_list, List neighbour_list, NumericVector events, NumericVector weights, DataFrame samples, double bw, Function kernel_func, DataFrame nodes, DataFrame line_list, bool verbose){
  //step 1 : mettre toutes les valeurs a 0
  NumericVector base_k = rep(0.0,samples.nrow());
  NumericVector base_count = rep(0.0,samples.nrow());

  //step0 : extraire tous les vecteurs necessaires
  NumericVector line_weights = line_list["weight"];
  IntegerVector samples_edgeid = samples["edge_id"];
  NumericVector samples_x = samples["X_coords"];
  NumericVector samples_y = samples["Y_coords"];
  IntegerVector samples_oid = samples["oid"];
  NumericVector nodes_x = nodes["X_coords"];
  NumericVector nodes_y = nodes["Y_coords"];

  //step2 : iterer sur chaque event
  int cnt_e = events.length();
  Progress p(cnt_e, verbose);
  for(int i=0; i < cnt_e; ++i){
    p.increment(); // update progress
    //#preparer les differentes valeurs de departs pour l'event y
    int y = events[i];
    double w = weights[i];
    NumericVector samples_k = esd_kernel_rcpp(y, edge_list,
                                              neighbour_list, bw,
                                              kernel_func, line_weights,
                                              samples_edgeid, samples_x,
                                              samples_y, samples_oid,
                                              nodes_x, nodes_y);
    base_k = base_k + samples_k;
    NumericVector count = rep(0.0,samples.nrow());
    LogicalVector test = samples_k>0;
    count[test] = 1.0;
    base_count = base_count + (count*w);

  };
  DataFrame df =  DataFrame::create( Named("sum_k") = base_k ,Named("n") = base_count);
  return df;
}

