#include<stdio.h>
#include <stdlib.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<map>
#include "graph.h"
#include "token.h"

using namespace std;

Graph DAG_tree;
int main(int argc, char ** argv) { 
  //Takes in a .dag.out file from unfoldsim.f90, outputs a dot file that contains the max 
  //path highlighted
    map<string, double> Nodes;
  double Tcutoff=0.01;
  if (argc < 2) {
    cerr << "Usage: ./maxTraffic <INPUT.dag.out> [traffic_cutoff=0.01 dagfile_url ]" << endl;
    cerr << "Converts a traffic file to a DOT file" << endl;
    cerr << "Optional cutoff prunes paths with low traffic" << endl;
    cerr << "maxTraffic.cpp v.  Wed Jun 18 14:57:20 EDT 2014 C.Bystroff" << endl;
    exit(1);
  }
  cout << "#CREATED_BY maxTraffic.cpp v.  Thu Aug  1 15:20:02 EDT 2013 " << endl;
  cout << "#Input DAG file  (argv[1]) =" <<  argv[1] << endl;

  ifstream INFILE(argv[1]);
  char* dagfile=argv[1];
  if (!INFILE) { cerr << "maxTraffix:: DAG file not found : " << dagfile << endl; exit(1) ; }
  dagfile=token(argv[1]);

  if (argc > 2) {
    Tcutoff = atof(argv[2]);
    cout << "#Minimum traffic to show (argv[2]) =" <<  Tcutoff << endl;
  } else {
    cout << "#Minimum traffic to show (default) =" << Tcutoff << endl;
  }

  if (argc > 3) {
    dagfile=argv[3];
    cout << "#Output URL file (argv[3]) =" << dagfile << endl;
  } else {
    cout << "#No output URL file " << endl;
  }

  string name;
  while ( !INFILE.eof() ) { 
    INFILE >> name;
    if (name == "TSTATE") { //now add the transition states
      string node;
      string from;
      string to1;
      string to2;
      string type;
      string discard;
      string discard2;
      double traffic;
      bool stop = false;
      // INFILE >> node >> from >> to1 >> to2 >> discard >> type >> discard2 >> traffic;
      INFILE >> node >> from >> to1 >> to2 >> discard >> type >> discard2 >> traffic;
      if (traffic > Tcutoff) { //pruning step 1
	    if ( node == "" ) //check if last line
	      break;
	    if (type == "p" || type == "b" || type == "h" || type == "s" ) { 
	      node = "t"+type+node;
	      from = "n"+from;
	      to1 = "n"+to1;
	      to2 = "n"+to2;
	      DAG_tree.addNode(node, traffic);
	      DAG_tree.addNode(from, 0); //itermediate
	      DAG_tree.addNode(to1, 0); //intermediate
	      if (type != "s" ) { DAG_tree.addNode(to2, 0); }  //intermediate
	      DAG_tree.addEdge(from, node); //from to tstate
	      DAG_tree.addEdge(node, to1); //from tstate to one child
	      if (type != "s" ) { DAG_tree.addEdge(node, to2); } //from tstate to other child
	    }
	    //cout << "the captured line is: transition " << node;
	    //cout << " goes from "<< from << " to nodes " << to1;
	    //cout << " and " << to2 << " with " << "traffic value of: ";
	    //cout << traffic << endl;
	    //}
      }
    }
  }
  cout << "#starting maxTraffic algorithm..." << endl;
  DAG_tree.getNode("n1")->highlight = true;
  DAG_tree.getNode("n1")->max_highlight = true;
  DAG_tree.maxTraffic( DAG_tree.getNode("n1") );
  INFILE.close();
  cout << "#number of nodes in maxTraffic: " << DAG_tree.max_val << endl;
  cout << "#Traffic cutoff used: " << Tcutoff << endl;
  cout << "digraph FILENAME" << endl;
  cout << "{" << endl;
  cout << "size=\"15, 25\"";
  DAG_tree.print_tree(DAG_tree.getNode("n1"), dagfile );
  cout << "}" << endl;
  return 0;
}
