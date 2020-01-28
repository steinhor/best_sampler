#ifndef parameterMap_CC
#define parameterMap_CC

#include "pratt_sampler/parametermap.h"
#include <iostream>

//Returns an integer from the map.
int CparameterMap::getI(string key,int def)
{
  int param;
  map<string,string>::iterator itr;
  itr = this->find(key);
  if(itr!=this->end()){
    stringstream ss(itr->second);
    ss>>param;
  }else{
    param = def;
  }
  return param;
}

//Returns an bool from the map.
bool CparameterMap::getB(string key,bool def)
{
  bool param=false;
  string pstring,ppstring;
  stringstream ss;
  map<string,string>::iterator itr;
  itr = this->find(key);
  if(itr!=this->end()){
    pstring=itr->second;
    //printf("pstring=%s\n",pstring.c_str());
    char tf[10];
    strncpy(tf,pstring.c_str(),1);
    //printf("tf=%s\n",tf);
    if(tf[0]=='t' || tf[0]=='1')
			param=true;
    else if(tf[0]=='f' || tf[0]=='0')
			param=false;
    else{
      printf("parameterMap::getB(), boolean parameter with key %s read in with value %s, set to false\n",key.c_str(),pstring.c_str());
    }
  }
  else{
    //printf("string %s not defined\n",key.c_str());
    param = def;
  }
  return param;
}

//Returns a string from the map.
string CparameterMap::getS(string key,string def)
{
  string param;
  map<string,string>::iterator itr;
  itr = this->find(key);  //find the value associated with string "key" in the parameter map
  if(itr!=this->end()){
    param = itr->second;
  }else{
    param = def;  //default to second string if not found
  }
  return param;
}

//Returns a double from the map.
double CparameterMap::getD(string key,double def)
{
  double param;
  map<string,string>::iterator itr;
  itr = this->find(key);
  if(itr!=this->end()){
    stringstream ss(itr->second);
    ss>>param;
  }else{
    param = def;
  }
  return param;
}

//Adds a double to the map.
void CparameterMap::set(string key,double val)
{
  string sval;
  stringstream ss;
  ss<<val;ss>>sval;
  set(key,sval);
}
//Adds an int to the map.
void CparameterMap::set(string key,int val)
{
  string sval;
  stringstream ss;
  ss<<val;ss>>sval;
  set(key,sval);
}
//Adds a bool to the map.
void CparameterMap::set(string key,bool val)
{
  string sval;
  stringstream ss;
  ss<<val;ss>>sval;
  set(key,sval);
}
//Adds a char* to the map.
void CparameterMap::set(string key,char* val)
{
  string sval(val);
  set(key,sval);
}
//Adds a string to the map.
void CparameterMap::set(string key,string val){
	this->insert(make_pair(key,val));
}

//Adds a vector to the map.
void CparameterMap::set( string key, vector< double > val){
  stringstream ss;
  for (int i=0;i<static_cast<int>(val.size());++i){
		ss<<"    "<<val[i]<<"\n";
	}
  set(key,ss.str());
}

//Adds a vector to the map.
void CparameterMap::set( string key, vector< string > val){
  stringstream ss;
  for (int i=0;i<static_cast<int>(val.size());++i) {ss<<"    "<<val[i]<<"\n";}
  set(key,ss.str());
}

//Adds a matrix to the map.
void CparameterMap::set( string key, vector< vector< double > > val){
  stringstream ss;
  for (int i=0;i<static_cast<int>(val.size());++i) {
    ss << "    ";
    for (int j=0;j<static_cast<int>(val[i].size());++j) {ss<<val[i][j]<<" ";}
    ss << "\n";
  }
  set(key,ss.str());
}

// Read parameters from file in format "type  key value", e.g.,
// int nqmax 50
// A # sign denotes a comment line
void CparameterMap::ReadParsFromFile(const char *filename){
  ifstream parsfile;
	//char line[150];
  string type,key,line,value;
	stringstream ss;
  parsfile.open(filename);
  if(! parsfile){
    printf("attempting to read non-existent parameter file %s\n",filename);
    exit(1);
  }

  while(!parsfile.eof()){
		//uses a string instead of a c string, removes line length limitation
    getline(parsfile, line, '\n');
    if(line.c_str()[0]!='#' && strlen(line.c_str())>4){
			ss.str("");
			ss.clear();
			ss << line;
			ss >> type;
			if(type=="double" || type=="int" || type=="bool" || type=="string")
				ss >> key;
			else
				key=type;

			//these lines allow for vector data to be read in, in the form of a set of delimited values.
			//this line sets value to be a string that is the rest of the line.
      value = ss.str().substr(ss.tellg());

			//these lines strip out any leading or trailing whitespace from the strings.
      size_t beginning = value.find_first_not_of(" \t");
      if(beginning == string::npos){
        cout << "error: blank parameter value for parameter " << key << "." << endl;
        exit(1);
      }
      size_t end = value.find_last_not_of(" \t");
      size_t range = end-beginning + 1;
      value = value.substr(beginning, range);

			// cout << "Storing:" << endl;
			// cout << "Key:" << key << endl;
			// cout << "Value:" << value << endl;
      set(key,value);
      ss.flush();
    }
  }
	// while(parsfile.getline(line,150)){
	// 		if(line[0]!='#' && strlen(line)>4){
	// 			ss << line;
	// 			ss >> type >> key >> value;
	// 			set(key,value);
	// 			ss.clear();
	// 		}
	// 	}
  parsfile.close();
}

void CparameterMap::ReadParsFromFile(string filename){
  ReadParsFromFile(filename.c_str());
}

void CparameterMap::PrintPars(){
  map<string,string>::iterator itr;
  for(itr=this->begin(); itr !=this->end(); ++itr){
    cout << itr->first << " = " << itr->second << endl;
  }
}

//Returns a matrix from the map
vector< vector< double > > CparameterMap::getM( string key, double def){
  vector< vector< double > > mtx;
  double tmp;
  map<string,string>::iterator itr;
  itr = this->find(key);
  if(itr!=this->end()){
    stringstream ss(itr->second);
    string line("");
    while (ss.good()){
      vector< double > vec;
      while ((line=="")&&ss.good()) {
        getline(ss,line);
        stringstream buf(line);
        while (buf>>tmp) {vec.push_back(tmp);}
      }
      line = "";
      if (vec.size()!=0) mtx.push_back(vec);
    }
  }else{
    vector< double > vec(1,0.);
    mtx.push_back(vec);
  }
  return mtx;
}

#endif
