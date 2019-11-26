#ifndef __parametermap_H__
#define __parametermap_H__
#include <map>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include "randy.h"

using namespace std;

//---------------------------------------------------------------
//This code helps read and parse parameters, and is based on an STL map
//
//Example:
// CparameterMap mymap;
// mymap.ReadParsFromFile("parameters.txt");
// arraysize = mymap.getI("ARRAYSIZE",100);   // 100 would be the default value
//
// Here parameters.txt might have a line, which would over-ride the defaul
// ARRAYSIZE 200
//
//---------------------------------------------------------------

class CparameterMap : public map<string,string> {
public:
  bool   getB(string ,bool);
  int    getI(string ,int);
  string getS(string ,string);
  double getD(string ,double);
  vector< double > getV(string, string);
  vector< string > getVS(string, string);
  vector< vector< double > > getM(string, double);
  void set(string, double);
  void set(string, int);
  void set(string, bool);
  void set(string, string);
  void set(string, char*);
  void set(string, vector< double >);
  void set(string, vector< string >);
  void set(string, vector< vector< double > >);
  void ReadParsFromFile(const char *filename);
  void ReadParsFromFile(string filename);
  void PrintPars();
};

#endif
