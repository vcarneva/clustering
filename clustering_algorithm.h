// members of the class clustering algorithm share the definition of 
// the functions list_members and sort_clusters, these two operate on
// two vectors cluster_indicators and map_labels that are supposed to
// have always the same structure and meaning
#include <iostream>

class clustering_algorithm {

protected:
  std::vector<int> cluster_indicators;
  std::vector<int> map_labels;

  void list_members ();
  void sort_clusters ();

public:
  int num_clusters;
  std::vector<std::list<int> > cluster_members;
  std::vector<int> cluster_size;
  std::vector<int> clusters_ranking;


};


class jarvis_patrick : public clustering_algorithm {

protected:
  std::vector<std::list<int> > nearest_neighbors;
  std::vector<std::list<double> > nearest_neighbors_distances;
  int size_of_nneigh_list;
  int shared_neigh_thres;

  void read_configurations ( std::vector< std::vector<double> > &  configs );
  void do_clustering ( std::vector< std::vector<double> > &  configs );
  void build_sorted_k_distance_plot ( std::vector< std::vector<double> > &  configs );

public:
  std::vector<std::list<double> > sorted_k_distances_lists;

  jarvis_patrick (int m, int p, std::vector< std::vector<double> > &  data_set);

};


class dbscan : public clustering_algorithm {

protected:
  std::vector<std::list<int> > epsilon_neighborhood;
  double epsilon;
  int minpts;

  void read_configurations ( std::vector< std::vector<double> > &  configs);
  void do_clustering ( std::vector< std::vector<double> > &  configs);

public:

  dbscan (double eps, int p, std::vector< std::vector<double> > &  data_set);

};
