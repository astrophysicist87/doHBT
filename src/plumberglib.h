#ifndef PLUMBERGLIB_H
#define PLUMBERGLIB_H

#include<string>
#include<fstream>

using namespace std;

string truestring = "true";
string falsestring = "false";

bool fexists(const char *filename)
{
  ifstream ifile(filename);
  return ifile;
}

#endif
