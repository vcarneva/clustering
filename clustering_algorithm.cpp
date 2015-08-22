#include <vector>
#include <list>
#include "clustering_algorithm.h"
#include <math.h>

// Two clustering methods are implemented : jarvis-patrick and DBSCAN (Density-Based Spatial Clustering 
// of Applications with Noise) . 
//
// jarvis patrick clustering algorithm : two elements of the dataset belong to the same cluster
// if each one is in the neighbor list of the other (m closest elements) AND at least p neighbors
// are shared (present in both lists).
//Jarvis RA, Patrick EA, Clustering Using a Similarity Measure Based on Shared Near Neighbors. 
//Computers, IEEE Transactions 22(11) 1025-1034, 1973.
//
// Density-Based Spatial Clustering of Applications with Noise : a point q is 
// in the same cluster as point p if q is density reachable from p, i.e. if their distance is 
// less than epsilon to p AND there are at least minpts points closer than
// epsilon around p.
// Ester M, Kriegel HP, Sander J, Xu Xi. A Density-Based Algorithm for Discovering Clusters in Large Spatial Databases with Noise.
// In Evangelos Simoudis, Jiawei Han, Usama M. Fayyad. Proceedings of the Second International 
// Conference on Knowledge Discovery and Data Mining (KDD-96). AAAI Press. pp. 226–231, 1996.



void clustering_algorithm::list_members ()
{

// this is just to build a list of members for each label

  std::list<int> cl_mem;
  
  for (int i=0; i<map_labels.size(); i++) {

    for (int j=0; j<cluster_indicators.size(); j++) {
      
      if ( cluster_indicators[j] == map_labels[i] ) {
        cl_mem.push_back (j);
      }

    } 
   
    cluster_members.push_back (cl_mem);
    cl_mem.clear(); 

  }

}



void clustering_algorithm::sort_clusters ()
{
 
// we want to sort the lists of members according to the size
// we will build clusters_ranking a sorted list of cluster indexes (not labels)
// first of all get the size of each cluster
 
  for (int i=0; i<cluster_members.size(); i++) {

    cluster_size.push_back (cluster_members[i].size());

  }

  std::vector<int> index_done;
  index_done.assign (cluster_size.size(), 0);

// now sort : double loop, not efficient but the number of clusters should be small enough

  int index;
  int maxsize;

  for (int i=0; i<cluster_size.size(); i++) {

    maxsize=0;

    for (int j=0; j<cluster_size.size(); j++) {

      if (cluster_size[j] > maxsize && index_done[j] != 1) {
      
        maxsize=cluster_size[j];
        index=j;
      
      }

    }

// every time we discover the maximum size we keep the index by pushing it into clusters_ranking 
// and we pull the index out of the pool by switching the corresponding value of index_done

    clusters_ranking.push_back (index);
    index_done[index]=1;

  }

// clusters_ranking is what we use for the output,
// e.g. to list the elements of the most populated cluster get the first element stored in clusters_ranking
// then use this index, let's say i, to access the elements of cluster_members[i] 

}


jarvis_patrick::jarvis_patrick (int m, int p, std::vector< std::vector<double> > &  data_set) {

  size_of_nneigh_list = m;
  shared_neigh_thres = p;
  read_configurations(data_set);
  do_clustering(data_set);
  list_members();
  sort_clusters();

}


void jarvis_patrick::read_configurations ( std::vector< std::vector<double> > &  configs)
{

  std::list<int> iniitalize_nneigh;
  std::list<double> iniitalize_nneigh_dist;
  iniitalize_nneigh.assign (size_of_nneigh_list, -1);
  iniitalize_nneigh_dist.assign (size_of_nneigh_list, 1.0e+30);

  for (size_t i = 0; i < configs.size(); i++) {

    nearest_neighbors.push_back (iniitalize_nneigh);
    nearest_neighbors_distances.push_back (iniitalize_nneigh_dist);  

    for (size_t j = 0; j < configs.size(); j++) {
//      double dist = *(configs[i]).dist (*(configs[j]));

//   use this for the moment 
     
     double dist = sqrt(pow((configs[i][0]-configs[j][0]),2)+pow((configs[i][1]-configs[j][1]),2)+pow((configs[i][2]-configs[j][2]),2)+pow((configs[i][3]-configs[j][3]),2)+pow((configs[i][4]-configs[j][4]),2)+pow((configs[i][5]-configs[j][5]),2)); 

      std::list<int>::iterator nn = nearest_neighbors[i].begin();
      std::list<double>::iterator nnd = nearest_neighbors_distances[i].begin();

// for each point we build the list of its nearest neighbors up to size_of_nneigh_list
// distances are stored in nearest_neighbors_distances in decreasing order

      for ( ; nn != nearest_neighbors[i].end(); nn++, nnd++) {
        if (dist > *nnd)
          break;
      }
      if ( (nnd != nearest_neighbors_distances[i].begin()) && (j != i)   )  {
        nearest_neighbors[i].insert (nn, (int)j);
        nearest_neighbors_distances[i].insert (nnd, dist);
        nearest_neighbors[i].pop_front();
        nearest_neighbors_distances[i].pop_front();
      }
    }
  }



}

void jarvis_patrick::build_sorted_k_distance_plot ( std::vector< std::vector<double> > &  configs ) {

// This is to create sorted-k-distance-plots, useful for analyzing the dataset and to guess the 
// paramters for the clustering, but not strictly required in any algorithm. 
// k-distance of element i is the distance of its kth nearest neighbor, the sorted list of all the 
// k-distances across the dataset is the k-distance plot 
// the lists will be accessed in reversed order because we need them in descending order 

  int position;
  std::list<double> list_of_distances;
  std::list<double>::reverse_iterator rit;
  for (int k = 0; k < size_of_nneigh_list; k++) {
    for (size_t i = 0; i < configs.size(); i++) {
      position = 0;
      for ( rit=nearest_neighbors_distances[i].rbegin() ; rit != nearest_neighbors_distances[i].rend(); ++rit ) {
        if (k == position) {
          list_of_distances.push_back (*rit);
          break;
        }
        position++;
      }
    }
    list_of_distances.sort();
    sorted_k_distances_lists.push_back (list_of_distances);
    list_of_distances.clear();
  }
}

void jarvis_patrick::do_clustering ( std::vector< std::vector<double> > &  configs )
{

// cluster_indicators will be used to store the cluster-id of each point  
// each id is defined as the lowest index in the cluster 
// initially each point is a cluster with 1 element

  for (int i = 0; i < (int)configs.size(); i++) {
    cluster_indicators.push_back (i);
  }

  for (size_t i = 0; i < configs.size(); i++) {
    std::list<int>::iterator nni;

    for (size_t j = i+1; j < configs.size(); j++) {
      std::list<int>::iterator nnj;
      bool i_neigh_j = false;

// first check if j is a neghbor of i : find j into nearest_neighbors[i]

      for (nni = nearest_neighbors[i].begin() ; nni != nearest_neighbors[i].end(); nni++) {
        if ( (int)j == *nni ) {
          i_neigh_j = true;
          break;
        } 
      }


      if (i_neigh_j) { 

        bool j_neigh_i =false;

// now check if i is a neghbor of j : remember that the relationship is not symmetric

        for (nnj = nearest_neighbors[j].begin() ; nnj != nearest_neighbors[j].end(); nnj++) {
          if ( (int)i == *nnj ) {
            j_neigh_i = true;
            break;
          } 
        }

        if ( i_neigh_j && j_neigh_i ) {
          int nshared_neighbors = 0;

// now let's count the number of nearest neighbors shared by i and j

          for (nni = nearest_neighbors[i].begin() ; nni != nearest_neighbors[i].end(); nni++) {
            for (nnj = nearest_neighbors[j].begin() ; nnj != nearest_neighbors[j].end(); nnj++) {
              if (*nni == *nnj && *nnj != i && *nni != j) {
                nshared_neighbors++;
              } 
            }
          }

// clustering now is only about changing the clid of the two points but we change them only if needed : 
// it is possible that the two points are already in the same cluster 
        
          if (nshared_neighbors >= shared_neigh_thres && cluster_indicators[i]!=cluster_indicators[j]) {  
  

// we keep the lowest index : remember that we have to change also all the instances of the old index
// in cluster_indicators

            if (cluster_indicators[j]>cluster_indicators[i]) {

              for (int k = 0; k < (int)configs.size() ; k++) {
                if (cluster_indicators[k] == cluster_indicators[j]) {
                  cluster_indicators[k]=cluster_indicators[i];
                }
              }

              cluster_indicators[j]=cluster_indicators[i];
            }
            else {

              for (int k = 0; k < (int)configs.size() ; k++) {
                if (cluster_indicators[k] == cluster_indicators[i]) {
                  cluster_indicators[k]=cluster_indicators[j];
                }
              }

              cluster_indicators[i]=cluster_indicators[j];
            }
 
          }

        }

      }
    
    } 

  }

// with this definition of clid each cluster possesses one and only one element 
// such that the index is equal to the clid : we use this property to count 
// the clusters and to build the list of clids (map_labels)

  num_clusters = 0;
  for (int i = 0; i < (int)configs.size(); i++) {

    if ( cluster_indicators[i]== i ){
      map_labels.push_back (cluster_indicators[i]);
      num_clusters++;
    }

  }

}


dbscan::dbscan (double eps, int p, std::vector< std::vector<double> > &  data_set) {

  epsilon = eps;
  minpts = p;
  read_configurations(data_set);
  do_clustering(data_set);
  list_members();
  sort_clusters();

}


void dbscan::read_configurations ( std::vector< std::vector<double> > &  configs)
{

  std::list<int> eps_neigh;
// for each point we build the list of points closer than epsilon (epsilon neighborhood) 
  for (size_t i = 0; i < configs.size(); i++) {

    for (size_t j = 0; j < configs.size(); j++) {
//      double dist = *(configs[i]).dist (*(configs[j]));
//   use this for the moment 
     
     double dist = sqrt(pow((configs[i][0]-configs[j][0]),2)+pow((configs[i][1]-configs[j][1]),2)+pow((configs[i][2]-configs[j][2]),2)); 

      if ( dist < epsilon )  {
     
        eps_neigh.push_back (j); 

      }

    }

  epsilon_neighborhood.push_back (eps_neigh);
  eps_neigh.clear();

  }

}



void dbscan::do_clustering ( std::vector< std::vector<double> > &  configs)
{

// cluster_indicators will be used to store the cluster-id of each point  
// each id is defined as the lowest index in the cluster 
// initially each point is a cluster with 1 element

  for (int i = 0; i < (int)configs.size(); i++) {
    cluster_indicators.push_back (i);
  }

  for (size_t i = 0; i < configs.size(); i++) {
    std::list<int>::iterator nni;

// first check if i is a core point : i possesses an epsilon neighborhood of size at least minpts

    bool core_point_i = false;
    if ( epsilon_neighborhood[i].size() >= minpts ) {
      core_point_i = true; 
    }

// if i is a core point then i + its epsilon neighborhood is a cluster and we change the labels of the
// points by keeping the lowest index

    if ( core_point_i ) {
      for (nni = epsilon_neighborhood[i].begin(); nni != epsilon_neighborhood[i].end(); nni++) {
// we change the index only if needed : it is possible that the two points are already in the same cluster 
        if (cluster_indicators[i]!=cluster_indicators[*nni]) {
          if (cluster_indicators[*nni]>cluster_indicators[i]) {
            for (size_t k = 0; k < configs.size() ; k++) {
              if (cluster_indicators[k] == cluster_indicators[*nni]) {
                cluster_indicators[k]=cluster_indicators[i];
              }
            }
            cluster_indicators[*nni]=cluster_indicators[i];
          }
          else {
            for (size_t k = 0; k < configs.size() ; k++) {
              if (cluster_indicators[k] == cluster_indicators[i]) {
                cluster_indicators[k]=cluster_indicators[*nni];
              }
            }
            cluster_indicators[i]=cluster_indicators[*nni];
          }
        }
      }
    }
 
    for (size_t j = i+1; j < configs.size(); j++) {
      std::list<int>::iterator nnj;

// now check if j is a core point 
      bool core_point_j = false;
      if ( epsilon_neighborhood[j].size() >= minpts ) {
        core_point_j = true; 
      }

// now check if i and j are neighbors, i.e. they are closer than epsilon (symmetric relationship) 
      bool neighbors = false;
      for ( nni = epsilon_neighborhood[i].begin(); nni != epsilon_neighborhood[i].end(); nni++) {
        if ( j == *nni ) {
          neighbors = true;
          break;
        } 
      }

      bool density_reachable = false;
      if ( neighbors && ( core_point_i || core_point_j ) ) {
        density_reachable = true;
      }

// the two points belong to the same cluster if they are neighbors and at least one of the 
// two is a core point (one of the points is density reachable from the other)
 
      if (density_reachable) { 
// if j is not a core point then i is a core point : add only j to the cluster containing i
        if ( !core_point_j ) {  
          if (cluster_indicators[i]!=cluster_indicators[j]) {
            if (cluster_indicators[j]>cluster_indicators[i]) {
              for (size_t k = 0; k < configs.size() ; k++) {
                if (cluster_indicators[k] == cluster_indicators[j]) {
                  cluster_indicators[k]=cluster_indicators[i];
                }
              }
              cluster_indicators[j]=cluster_indicators[i];
            }
            else {
              for (size_t k = 0; k < configs.size() ; k++) {
                if (cluster_indicators[k] == cluster_indicators[i]) {
                  cluster_indicators[k]=cluster_indicators[j];
                }
              }
              cluster_indicators[i]=cluster_indicators[j];
            }
          }
        }
// if j is a core point then we add j + its epsilon neighborhood to the cluster containing i
        else {
          for (nnj = epsilon_neighborhood[j].begin(); nnj != epsilon_neighborhood[j].end(); nnj++) {
            if (cluster_indicators[i]!=cluster_indicators[*nnj]) {
              if (cluster_indicators[*nnj]>cluster_indicators[i]) {
                for (size_t k = 0; k < configs.size() ; k++) {
                  if (cluster_indicators[k] == cluster_indicators[*nnj]) {
                    cluster_indicators[k]=cluster_indicators[i];
                  }
                }
                cluster_indicators[*nnj]=cluster_indicators[i];
              }
              else {
                for (size_t k = 0; k < configs.size() ; k++) {
                  if (cluster_indicators[k] == cluster_indicators[i]) {
                    cluster_indicators[k]=cluster_indicators[*nnj];
                  }
                }
                cluster_indicators[i]=cluster_indicators[*nnj];
              }
            }
          }
        }
      }
    }
  }

// with this definition of clid each cluster possesses one and only one element 
// such that the index is equal to the clid : we use this property to count 
// the clusters and to build the list of clids (map_labels)

  num_clusters = 0;
  for (size_t i = 0; i < configs.size(); i++) {
    if (cluster_indicators[i]==i){
      map_labels.push_back (cluster_indicators[i]);
      num_clusters++;
    }
  }

}


