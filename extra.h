#ifndef _EXTRA_H_INCLUDED_
#define _EXTRA_H_INCLUDED_

#include <fstream>
#include <Eigen/Dense>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <string>

using namespace std;
using namespace Eigen;

void matrix_output(MatrixXcd Mc, MatrixXcd Mcx, MatrixXcd Mcy, MatrixXcd Mcz)
{
  cout << setprecision(3) <<  Mc << endl << endl << Mcx << endl << endl << Mcy << endl << endl << Mcz << endl << endl;
}

void generate_scriptname(string& scriptname) {
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer,sizeof(buffer),"%S-%M-%I-%Y-%m-%d",timeinfo);
    string str(buffer);
    string base= "wolframscripts/script";
    string extension=".nb";
    scriptname= base+str+extension;
}

void eigenvalues_Mathematica(MatrixXcf Mc, ofstream& fout, string scriptname)
{
  int size = Mc.rows()/2;
  fout.open(scriptname);
  //fout << "#!/usr/local/bin/WolframScript -linewise -script" << endl;
  fout << "Print[Eigenvalues[N[{";
  for(int i=0; i<2*size; i++)
  {
    fout << "{ ";
    for(int j=0; j<2*size; j++)
      {
        fout << Mc.real()(i,j) << "+ I (" << Mc.imag()(i,j);
         if(j==2*size-1) fout << ") ";
         else fout << "),";
      }
   if(i==2*size-1) fout << "} ";
   else fout << "},";
  }
  fout << "}]]];" << endl;
  fout.close();
}

void make_exec(string str)
{
   const char *scriptname = str.c_str();
   chmod(scriptname, S_IRWXU);
   cout << "To generate the eigenvalues using Mathematica for verification, run "<< str << endl;
}

void generate_wlscript(MatrixXcf Mc, ofstream& fout)
{
  int size = Mc.rows()/2;
  string scriptname;
  generate_scriptname(scriptname);
  eigenvalues_Mathematica(Mc, fout, scriptname);
  make_exec(scriptname);
}

void progress_percent_desc(double initial_temp, double final_temp, double temperature)
{
   int progress_percent;
   progress_percent=int((final_temp-temperature)*100/(final_temp-initial_temp));
   cout.flush();
   cout << "\r [ "<< progress_percent+1 << "% ] ";
   for(int i=0; i<progress_percent; i++)
   {
     cout << "#";
   }
}

void progress_percent_asc(double initial_temp, double final_temp, double temperature)
{
   int progress_percent;
   progress_percent=int((temperature-initial_temp)*100/(final_temp-initial_temp));
   cout.flush();
   cout << "\r [ "<< progress_percent+1 << "% ] ";
   for(int i=0; i<progress_percent; i++)
   {
     cout << "#";
   }
}

void spinarrangement_Mathematica_output(MatrixXd M)
{
  string notebookname; generate_scriptname(notebookname);
  ofstream fout(notebookname);
  double a;
  cout << "Enter lattice separation (a): ";
  cin >> a;

  fout << "Show[Graphics3D[{" << endl;
  for(int i=0; i< M.rows(); i++)
  {
    fout << "Arrow[{{" << a*i << ", 0, 0}, {" << M(i,0)+a*i << ","  << M(i,1) << "," << M(i,2) << "}}]";
    if(i!=M.rows()-1) fout << ",\n";
  }
  fout <<"}] ]";
}

void show_time(milliseconds begin_ms, milliseconds end_ms, string s)
{
   long int t = (end_ms.count()-begin_ms.count())/1000;
    if (t<=60)
    { cout <<  s << " took " << t << " seconds." << endl; }
    else if (t>60 && t<3600)
    {
      int minutes=int(t/60); int seconds= t%minutes;
      cout << s << " took " << minutes << " minute and " << seconds << " seconds." << endl;
    }
    else if(t>=3600)
    {
      int hrs= int(t/3600); int minutes=int((t-3600*hrs)/60); int seconds= int(t-3600*hrs-60*minutes);
      cout << s << " took " << hrs << " hour, " << minutes << " minutes and " << seconds << " seconds. ";
    }
    else
    {cout << s << " took " << t << "time. Wrong t received.\n"; }
}

#endif
