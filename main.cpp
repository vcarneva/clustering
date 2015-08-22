#include <stdio.h>
#include <vector>
#include <list>
#include "clustering_algorithm.h"

int main ()
{
  FILE * coor;
  std::vector< std::vector<double> > configuration;

  std::vector<double> cartesian (6, 0.0);
  coor = fopen ("coordinates.txt","r");
  while (fscanf (coor, "%lf %lf %lf %lf %lf %lf", &cartesian[0], &cartesian[1], &cartesian[2], &cartesian[3], &cartesian[4], &cartesian[5])!=EOF ) {
    configuration.push_back(cartesian);
  }
  fclose (coor);

// pass to the constractor the two parameters (by value) and the configuration object (by reference) 
  jarvis_patrick jp (8, 3, configuration);
  dbscan db (3.0, 5, configuration);  

// this is only for testing   
  printf ("*********************************************************\n\n");
  printf ("JARVIS PATRICK\n\n");
  printf ("There are %d clusters\n\n", jp.num_clusters);

  for (size_t j = 0 ; j < jp.clusters_ranking.size() ; j++){
    printf ("There are %d elements belonging to cluster %d:\n", (int)jp.cluster_members[jp.clusters_ranking[j]].size(),(int)j);
    printf ("begin cluster\n");

    for (std::list<int>::iterator i = jp.cluster_members[jp.clusters_ranking[j]].begin(); i != jp.cluster_members[jp.clusters_ranking[j]].end(); i++ ) {
      printf (" %d ", *i);
    }

    printf ("\nend cluster\n\n");
  }

  printf ("*********************************************************\n\n");
  printf ("DBSCAN\n\n");
  printf ("There are %d clusters\n\n", db.num_clusters);

  for (size_t j = 0 ; j < db.clusters_ranking.size() ; j++){
    printf ("There are %d elements belonging to cluster %d:\n", (int)db.cluster_members[db.clusters_ranking[j]].size(),(int)j);
    printf ("begin cluster\n");

    for (std::list<int>::iterator i = db.cluster_members[db.clusters_ranking[j]].begin(); i != db.cluster_members[db.clusters_ranking[j]].end(); i++ ) {
      printf (" %d ", *i);
    }

    printf ("\nend cluster\n\n");
  }
  printf ("*********************************************************\n\n");
//end testing

}
