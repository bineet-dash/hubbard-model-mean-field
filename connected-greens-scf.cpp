#include "tca_library.hpp"
#include "extra.hpp"

int Lc;

void show_eigenvalues(MatrixXcd H)
{
  std::vector<double> eigenvalues;
  diagonalize(H, eigenvalues);
  for(auto it=eigenvalues.begin(); it!= eigenvalues.end(); it++) cout << *it << " ";
  cout << endl << endl;
}

// void greens_sigma_generate(MatrixXd& suggested_randsigma, int lattice_index, long & idum)
// {
//   if(ran0(&idum)<=0.5) suggested_randsigma(lattice_index,2) *= -1;
// }

// MatrixXcd g_d(MatrixXcd H, double omega)
// {
//   MatrixXcd Z = cd(omega,eta)*MatrixXcd::Identity(H.rows(),H.rows());
//   MatrixXcd G = invert(Z-H);
//   return G;
// }

MatrixXcd g_c_L(MatrixXcd Hl, double omega, MatrixXcd Gr)
{
  MatrixXcd Z = cd(omega,eta)*MatrixXcd::Identity(Hl.rows(),Hl.rows());

  // MatrixXcd tau = MatrixXcd::Zero(Gr.rows(),1);
  // tau(0) = tau(Lc) = -t;
  // cd c_part = (tau.adjoint()*(Gr*tau)).sum(); //1x1 matrix, so sum() converts to cd
  // MatrixXcd corr = MatrixXcd::Zero(Hl.rows(),Hl.cols());
  // corr(Lc-1,Lc-1) = corr(2*Lc-1,2*Lc-1) = c_part;
  // cout << corr << endl << endl;

  MatrixXcd tau = MatrixXcd::Zero(Gr.rows(),Gr.rows());
  // tau(tau.rows()-1,0) = -t;
  tau(Lc-1,0) = tau(2*Lc-1, Lc) = -t;
  // cout << tau << endl << endl;
  // cout << Z << endl << endl << Hl << endl << endl <<  tau*Gr*tau.adjoint() << endl;
  // exit(1);
  MatrixXcd G = invert(Z- Hl- tau*Gr*tau.adjoint());
  return G;
}

MatrixXcd g_c_R(MatrixXcd Hr, double omega, MatrixXcd Gl)
{
  MatrixXcd Z = cd(omega,eta)*MatrixXcd::Identity(Hr.rows(),Hr.rows());
  MatrixXcd G_d_r = invert(Z-Hr); //will be used again and again.

  MatrixXcd tau = MatrixXcd::Zero(Gl.rows(),Hr.cols());
  // tau(tau.rows()-1,0) = -t;
  tau(Lc-1,0) = tau(2*Lc-1, Lc) = -t;

  MatrixXcd G = G_d_r + G_d_r*tau.adjoint()*Gl*tau*G_d_r;
  return G;
}


double scf_spectral_weight(MatrixXcd G, double k, double omega)
{
  // MatrixXcd G = g_d(H,omega);
  double weight = G.trace().imag();
  return -1/M_PI*weight;
}

double get_dos(MatrixXcd G, double omega)
{
  double DoS = 0;
  for(int n=0; n<size; n++)
  {
    double k = n*M_PI/(size*a);// -M_PI/a ;
    DoS += scf_spectral_weight(G, k, omega);
  }
  return DoS/size;
}

// int write_dos_count_1 = 0;
// int write_dos_count_2 = 0;
//
// void write_dos(MatrixXcd H, std::vector<MatrixXcd>& v, int half)
// {
//   int& num = (half==1)?write_dos_count_1:write_dos_count_2;
//   // if(half==1) num=write_dos_count_1;
//   // else if(half==2) num=write_dos_count_2;
//
//   string filename = "scf/scf_greens_test_"+ to_string(half)+ "_" + to_string(num)+ ".txt";
//   ofstream fout(filename);
//   int count = 0;
//
//   for(double omega=omega_L; omega<=omega_U; omega+=0.1)
//   {
//     MatrixXcd Gfunc = g_c(H, omega, v.at(count), half);
//     fout << omega << " " << get_dos(Gfunc, omega) << endl;
//     v.at(count) = Gfunc;
//     count++;
//   }
//   fout.close();
//   num++;
// }

int num_L = 0;
int num_R = 0;

void write_L_dos(MatrixXcd H, std::vector<MatrixXcd> G_r_list, std::vector<MatrixXcd>& G_l_list)
{
  string filename = "scf/scf_greens_test_L_" + to_string(num_L)+ ".txt";
  ofstream fout(filename);
  int count = 0;

  for(double omega=omega_L; omega<=omega_U; omega+=0.1)
  {
    MatrixXcd Gfunc = g_c_L(H, omega, G_r_list.at(count));
    fout << omega << " " << get_dos(Gfunc, omega) << endl;
    G_l_list.at(count) = Gfunc;
    count++;
  }
  fout.close();
  num_L++;
}

void write_R_dos(MatrixXcd H, std::vector<MatrixXcd> G_l_list, std::vector<MatrixXcd>& G_r_list)
{
  string filename = "scf/scf_greens_test_R_" + to_string(num_R)+ ".txt";
  ofstream fout(filename);
  int count = 0;

  for(double omega=omega_L; omega<=omega_U; omega+=0.1)
  {
    MatrixXcd Gfunc = g_c_R(H, omega, G_l_list.at(count));
    fout << omega << " " << get_dos(Gfunc, omega) << endl;
    G_r_list.at(count) = Gfunc;
    count++;
  }
  fout.close();
  num_R++;
}

int main(int argc, char* argv[])
{

  cout << "Enter the size, U: ";
  cin >> size >> U;
  Lc = size/2;

  std::vector<MatrixXcd> Gr_list;
  std::vector<MatrixXcd> Gl_list;
  for(double omega=omega_L; omega<=omega_U; omega+=0.1)
  {
    Gr_list.push_back(MatrixXcd::Zero(2*Lc,2*Lc));
    Gl_list.push_back(MatrixXcd::Zero(2*Lc,2*Lc));
  }

  MatrixXd randsigma = MatrixXd::Zero(size, 3);
  long idum = -1; //time(NULL);
  for(int i=0; i<size; i++)  randsigma(i,2) = 1;
  // for(int i=0; i<size; i++)  greens_sigma_generate(randsigma, i, idum);
  // cout << "spins\n\n" << randsigma << endl << endl;

  MatrixXd selected_randsigma = MatrixXd::Zero(Lc,3);

  int lattice_index = 0;
  for(int i=0; i<Lc; i++) selected_randsigma.row(i) = randsigma.row((lattice_index+i)%size);
  MatrixXcd H_cl_1 = cluster_h0(Lc,0)-U/2*cluster_sigmaz(Lc,selected_randsigma);

  // cout << randsigma.transpose() << endl << endl << selected_randsigma.transpose() << endl << endl;
  // cout << H_cl_1 << endl << endl; exit(11);

  lattice_index += Lc;
  for(int i=0; i<Lc; i++) selected_randsigma.row(i) = randsigma.row((lattice_index+i)%size);
  MatrixXcd H_cl_2 = cluster_h0(Lc)-U/2*cluster_sigmaz(Lc,selected_randsigma);



  for(int master_loop=0; master_loop<50; master_loop++)
  {
    std::cerr << "master_loop=" << master_loop << '\n';
    write_L_dos(H_cl_1, Gr_list, Gl_list);
    std::cerr << "H_cl_1 complete" << '\n';
    write_R_dos(H_cl_2, Gl_list, Gr_list);
    std::cerr << "H_cl_2 complete" << '\n';
  }

  exit(1943);

  // int no_sweeps = 20;
  // int initial_exp = -2;
  // int final_exp = -1;
  //
  // double final_temp = 10*pow(10,final_exp);
  //
  // milliseconds begin_ms, end_ms;
  // begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
  //
  // MatrixXd randsigma = MatrixXd::Zero(size, 3);
  // long idum = time(NULL);
  // for(int i=0; i<size; i++)  randsigma(i,2) = 1;
  // for(int i=0; i<size; i++)  greens_sigma_generate(randsigma, i, idum);
  //
  //
  // string filename, latticedata;
  // latticedata = "_U="+to_string(int(U))+"_size="+to_string(size)+"_sweeps="+to_string(no_sweeps)+"_Lc="+to_string(Lc);
  // filename="wolframscripts/tca_spin_arrangement"+current_time_str()+latticedata+".nb"; ofstream outfile_spinarr(filename);
  // spinarrangement_Mathematica_output(randsigma,outfile_spinarr);
  // filename="data/tca_m_length_vs_temp"+ current_time_str()+latticedata+".txt"; ofstream outfile_mlength(filename);
  // // filename="data/free_energy_vs_temp"+current_time_str()+latticedata+".txt"; ofstream outfile_freeenergy(filename);
  // // filename="data/mcdetails"+current_time_str()+latticedata+".txt"; ofstream outfile_mcdetails(filename);
  // cout << "==============================\n"<< "filename is: " << filename << "\n========================\n";
  //
  // for(int j=final_exp; j>=initial_exp; j--)
  // {
  //   for(double i=10; i>=2; i-=1)
  //   {
  //     double temperature = i*pow(10,j);
  //     for(int sweep=0; sweep<0.75*no_sweeps; sweep++)
  //     {
  //       for(int lattice_index=0; lattice_index<size; lattice_index++)
  //       {
  //         MatrixXd selected_randsigma = MatrixXd::Zero(Lc,3);
  //         for(int i=0; i<Lc; i++) selected_randsigma.row(i) = randsigma.row((lattice_index+i)%size);
  //         MatrixXcd H_cl = cluster_h0(Lc)-U/2*cluster_sigmaz(Lc,selected_randsigma);
  //         double free_energy = //find_canonical_free_energy(H_cl, temperature, randsigma);
  //                               filled_E(H_cl,Lc);
  //
  //         MatrixXd suggested_randsigma_selected = selected_randsigma;
  //         greens_sigma_generate(suggested_randsigma_selected,0, idum); //present lattice index is 0 in selected part
  //         MatrixXcd suggested_H_cl = cluster_h0(Lc)-U/2*cluster_sigmaz(Lc,suggested_randsigma_selected);
  //         double suggested_free_energy = //find_canonical_free_energy(suggested_H_cl, temperature, randsigma);
  //                                       filled_E(suggested_H_cl,Lc);
  //
  //         double uniform_rv = ran0(&idum); double move_prob = exp((free_energy - suggested_free_energy)/temperature);
  //         if(uniform_rv <= move_prob)
  //         {
  //           free_energy = suggested_free_energy;
  //           randsigma.row(lattice_index) = suggested_randsigma_selected.row(0);
  //         }
  //         else
  //         {
  //           suggested_randsigma_selected=selected_randsigma;
  //         }
  //       }
  //       cout << "\r sweep = " << sweep << " done."; cout.flush();
  //     }
  //
  //     for(int sweep= int(0.75*no_sweeps); sweep<no_sweeps; sweep++)
  //     {
  //       for(int lattice_index=0; lattice_index<size; lattice_index++)
  //       {
  //         MatrixXd selected_randsigma = MatrixXd::Zero(Lc,3);
  //         for(int i=0; i<Lc; i++) selected_randsigma.row(i) = randsigma.row((lattice_index+i)%size);
  //         MatrixXcd H_cl = cluster_h0(Lc)-U/2*cluster_sigmaz(Lc,selected_randsigma);
  //         double free_energy = //find_canonical_free_energy(H_cl, temperature, randsigma);
  //                              filled_E(H_cl,Lc);
  //
  //         MatrixXd suggested_randsigma_selected = selected_randsigma;
  //         greens_sigma_generate(suggested_randsigma_selected,0, idum); //present lattice index is 0 in selected part
  //         MatrixXcd suggested_H_cl = cluster_h0(Lc)-U/2*cluster_sigmaz(Lc,suggested_randsigma_selected);
  //         double suggested_free_energy =// find_canonical_free_energy(suggested_H_cl, temperature, randsigma);
  //                                       filled_E(suggested_H_cl,Lc);
  //
  //         double uniform_rv = ran0(&idum); double move_prob = exp((free_energy - suggested_free_energy)/temperature);
  //         if(uniform_rv <= move_prob)
  //         {
  //           free_energy = suggested_free_energy;
  //           randsigma.row(lattice_index) = suggested_randsigma_selected.row(0);
  //         }
  //         else
  //         {
  //           suggested_randsigma_selected=selected_randsigma;
  //         }
  //       }
  //       cout << "\r sweep = " << sweep << " done."; cout.flush();
  //     }
  //
  //     outfile_mlength << setw(5) << temperature << " ";
  //     for(int j=0; j<size; j++) outfile_mlength << " " << setw(5) << randsigma(j,2)  << " ";
  //     outfile_mlength << endl;
  //     // outfile_mlength << " \t" << m_length_avg/double(size) <<  endl;
  //     // outfile_freeenergy << temperature << " " << final_free_energy/double(count_free_energy) << endl;
  //
  //     cout << "\rtemperature = " << temperature << " done."; cout.flush();
  //   }
  // }
  //
  // cout << endl;
  // end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
  // show_time(begin_ms, end_ms,"TCA calculation");
  // spinarrangement_Mathematica_output(randsigma,outfile_spinarr);
  //
  // outfile_mlength.close();
  // // outfile_mcdetails.close();
  // // outfile_freeenergy.close();
  // outfile_spinarr.close();
  return 0;
}
