#ifndef GREEN_LIBRARY_HPP_INCLUDED
#define GREEN_LIBRARY_HPP_INCLUDED

#include "spa_library.hpp"
#include "extra.hpp"
#include "tca_library.hpp"
#include <iomanip>

milliseconds begin_ms, end_ms;

double omega_L = -8;
double omega_U = 8;

double eta = 0.05;
double a = 1;

MatrixXcd cluster_h0(int Lc, int pbc_pref = 1)
{
  MatrixXcd Mc = MatrixXcd::Zero(2*Lc,2*Lc);
  for(int row=0; row<2*Lc-1; row++) Mc(row,row+1)=Mc(row+1,row)=-t;
  Mc(Lc,Lc-1)= Mc(Lc-1,Lc)=0; //up and down spins don't mix

  if (pbc_pref==1) {
    Mc(2*Lc-1,Lc)=Mc(Lc,2*Lc-1)= -t; //PBC
    Mc(Lc-1,0)=Mc(0,Lc-1)=-t; //PBC
  }
  return Mc;
}

MatrixXcd cluster_sigmaz(int Lc, MatrixXd rs)
{
  MatrixXcd Mcz = MatrixXcd::Zero(2*Lc,2*Lc);
  for(int row=0; row<Lc; row++)  Mcz(row,row)= cd(rs(row,2),0);
  for(int row=Lc; row<2*Lc; row++) Mcz(row, row)=cd(-rs(row-Lc,2),0);
  return Mcz;
}

void greens_sigma_generate(MatrixXd& suggested_randsigma, int lattice_index, long & idum)
{
  if(ran0(&idum)<=0.5) suggested_randsigma(lattice_index,2) *= -1;
}

MatrixXcd invert(MatrixXcd A)
{
  int N = A.rows();
  int *IPIV = new int[N+1];
  int LWORK = N*N;
  __complex__ double* WORK= new __complex__ double [LWORK];
  int INFO;

  zgetrf_(&N,&N, reinterpret_cast <__complex__ double*> (A.data()), &N,IPIV,&INFO);
  zgetri_(&N, reinterpret_cast <__complex__ double*> (A.data()), &N,IPIV,WORK,&LWORK,&INFO);

  delete IPIV;
  delete WORK;
  if(INFO==0) return A;
  else {cout << "Inversion failed. Exiting..."; exit(28);}
}

// MatrixXd invert(MatrixXd A)
// {
//   int N = A.rows();
//   int *IPIV = new int[N+1];
//   int LWORK = N*N;
//   double* WORK= new double [LWORK];
//   int INFO;

//   dgetrf_(&N,&N, A.data(),&N,IPIV,&INFO);
//   dgetri_(&N, A.data(), &N,IPIV,WORK,&LWORK,&INFO);

//   delete IPIV;
//   delete WORK;
//   if(INFO==0) return A;
//   else {cout << "Inversion failed. Exiting..."; exit(28);}
// }

// double spectral_weight(MatrixXcd H, double k, double omega)
// {
//   MatrixXcd Z = cd(omega,eta)*MatrixXcd::Identity(H.rows(),H.rows());
//   MatrixXcd G = invert(Z-H);
//   double weight = G.trace().imag();
//   return -1/M_PI*weight;
// }

// double dos(MatrixXcd H, double omega)
// {
//   double DoS = 0;
//   for(int n=0; n<size; n++)
//   {
//     double k = n*M_PI/(size*a);// -M_PI/a ;
//     DoS += spectral_weight(H, k, omega);
//   }
//   return DoS/size;
// }

double dos(MatrixXcd H, double omega)
{
  MatrixXcd Z = cd(omega,eta)*MatrixXcd::Identity(H.rows(),H.rows());
  MatrixXcd G = invert(Z-H);
  double DoS = -1/M_PI*G.trace().imag();
  return DoS;
}

double mu(MatrixXcd H)
{
  double step = 0.02;
  double trapez_sum = 0.0; double omega = omega_L;
  while(trapez_sum < size)
  {
    trapez_sum += dos(H, omega)*step;
    omega += step;
  }
  return omega;
}

double check_e(MatrixXcd H)
{
  double ne=0.0;
  for(double omega=omega_L; omega<omega_U; omega += 0.01)
  {
    ne += dos(H,omega)*0.01;
  }
  cout << ne << endl;
  exit(1);
}

double filled_E(MatrixXcd H, int fill=size)
{
  double step = 0.04;
  double energy = 0.0;
  double no_of_electrons = 0.0;
  double omega;
  // begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());

  for(omega = omega_L; omega <= omega_U; omega += step)
  {
    if(no_of_electrons < fill)
    {
      energy += omega*dos(H, omega)*step;
      no_of_electrons += dos(H, omega)*step;

      // end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
      // cout << setw(5) << omega << " " << setw(5) << (end_ms.count()-begin_ms.count()) << endl;
    }
    else break;
  }
  return energy;
}

MatrixXd select_block(MatrixXd m, int li, int p)
{
  if(li+p < m.rows()){
    return m.block(li,0,p,3);
  }
  else{
    MatrixXd trun = MatrixXd::Zero(p,m.cols());
    trun.topRows(m.rows()-li) = m.bottomRows(m.rows()-li);
    trun.bottomRows(li+p-m.rows()) = m.topRows(li+p-m.rows());  
    return trun;
  }
}

double filled_E_ed_connected(MatrixXd rs, double temperature, int li, int Lc) //li ~ lattice_index
{
  // double step = 0.04;
  // double energy = 0.0;
  // double no_of_electrons = 0.0;
  // double omega;
  // int count_for_omega = 0;
  // begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());

  MatrixXd selected_rs_l = select_block(rs,li,Lc);
  MatrixXd selected_rs_r = select_block(rs, (li+Lc)%size, Lc);
  MatrixXcd Hl = cluster_h0(Lc)-U/2*cluster_sigmaz(Lc,selected_rs_l);
  MatrixXcd Hr = cluster_h0(Lc)-U/2*cluster_sigmaz(Lc,selected_rs_r);

  
  MatrixXcd tau = MatrixXcd::Zero(Hl.rows(),Hr.cols());
  tau(Lc-1,0) = tau(2*Lc-1, Lc) = -t;
  tau(0,Lc-1) = tau(Lc, 2*Lc-1) = -t; //outer PBC


  MatrixXcd H = MatrixXd::Zero(2*Hl.rows(),2*Hl.cols());
  H.block(0,0,Hl.rows(),Hl.cols()) = Hl;
  H.block(Hl.rows(),Hl.cols(),Hr.rows(),Hr.cols()) = Hl;
  H.block(0,Hl.cols(), Hl.rows(),Hr.cols()) = tau;
  H.block(Hl.rows(),0, Hr.rows(),Hl.cols()) = tau.adjoint();

  // cout << Hl << endl << endl << Hr << endl << endl << tau << endl << endl;
  // cout << H << endl << endl;
  // exit(1);

  double free_energy = find_canonical_free_energy(H,temperature,rs);
  

  // end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
  // cout << setw(5) << omega << " " << setw(5) << (end_ms.count()-begin_ms.count()) << endl;
  return free_energy;
}

double filled_E_ed_disconnected(MatrixXd rs, double temperature, int li, int Lc) //li ~ lattice_index
{

  MatrixXd selected_rs_l = select_block(rs,li,Lc);
  MatrixXcd Hl = cluster_h0(Lc)-U/2*cluster_sigmaz(Lc,selected_rs_l);

  double free_energy = find_canonical_free_energy(Hl,temperature,rs);
  return free_energy;
}


double get_mu_gf(vector <double> a_w, double dw, double temperature, int fill)
{
  double lb = omega_L;
  double ub = omega_U;
  double mu, ne; 
  // int loop=0;
  do
    {
      ne = 0.0;
      mu = (lb+ub)/2;

      // if(loop > 10) exit(1812);
      // cout << "lb=" << lb << " ub=" << ub  << " mu=" << mu;

      for(auto it=a_w.begin(); it!=a_w.end(); it++) ne += (*it)*dw/(exp((*it-mu)/temperature)+1);
      if(ne<fill) lb = mu;
      else if(ne > fill) ub = mu;

      // loop++;
      // cout << " ne=" << ne << endl;
    }
  while(abs(ne-fill)>0.001);
  return mu;
}

double get_mu_gf(MatrixXcd H, double temperature, int fill)
{
  vector <double> a_w;
  double omega_step = 0.02;
  for(double omega=omega_L; omega<omega_U; omega+= omega_step) a_w.push_back(dos(H,omega));
  return get_mu_gf(a_w, omega_step, temperature, fill);
}

double filled_E_disconnected(MatrixXd rs, int li, int Lc, double temperature) //li ~ lattice_index
{
  MatrixXd selected_rs_l = select_block(rs,li,Lc);
  MatrixXd selected_rs_r = select_block(rs, (li+Lc)%size, Lc);
  MatrixXcd Hl = cluster_h0(Lc)-U/2*cluster_sigmaz(Lc,selected_rs_l);
  MatrixXcd Hr = cluster_h0(Lc)-U/2*cluster_sigmaz(Lc,selected_rs_r);
  
  MatrixXcd tau = MatrixXcd::Zero(Hl.rows(),Hr.cols());
  tau(Lc-1,0) = tau(2*Lc-1, Lc) = -t;
  tau(0,Lc-1) = tau(Lc, 2*Lc-1) = -t; //outer PBC

  vector <double> a_w;
  vector <double> omega_arr;
  double step = 0.01;
  for(double omega = omega_L; omega <= omega_U; omega += step) omega_arr.push_back(omega);

  for(auto it=omega_arr.begin(); it!=omega_arr.end(); it++)
  {  
    MatrixXcd Z = cd(*it,eta)*MatrixXcd::Identity(Hl.rows(),Hl.rows());
    MatrixXcd Gl = invert(Z-Hl);
    a_w.push_back(-1/M_PI*Gl.trace().imag());
  }
  double mu = get_mu_gf(a_w,step, temperature,Lc);
  double energy = 0.0;
  for(int i=0; i<omega_arr.size(); i++) energy += omega_arr.at(i)*a_w.at(i)/(exp(omega_arr.at(i)-mu)/temperature+1) ;
  // cout << "temperature = " << temperature << ", mu = " << mu << endl;
  return energy;
}


double filled_E_connected(MatrixXd rs, int li, int Lc, double temperature) //li ~ lattice_index
{
  MatrixXd selected_rs_l = select_block(rs,li,Lc);
  MatrixXd selected_rs_r = select_block(rs, (li+Lc)%size, Lc);
  MatrixXcd Hl = cluster_h0(Lc)-U/2*cluster_sigmaz(Lc,selected_rs_l);
  MatrixXcd Hr = cluster_h0(Lc)-U/2*cluster_sigmaz(Lc,selected_rs_r);
  
  MatrixXcd tau = MatrixXcd::Zero(Hl.rows(),Hr.cols());
  tau(Lc-1,0) = tau(2*Lc-1, Lc) = -t;
  tau(0,Lc-1) = tau(Lc, 2*Lc-1) = -t; //outer PBC

  vector <double> a_w;
  vector <double> omega_arr;
  double step = 0.01;
  for(double omega = omega_L; omega <= omega_U; omega += step) omega_arr.push_back(omega);

  for(auto it=omega_arr.begin(); it!=omega_arr.end(); it++)
  {  
    MatrixXcd Z = cd(*it,eta)*MatrixXcd::Identity(Hl.rows(),Hl.rows());
    MatrixXcd Gl = invert(Z-Hl-tau*invert(Z-Hr)*tau.adjoint());
    a_w.push_back(-1/M_PI*Gl.trace().imag());
  }
  double mu = get_mu_gf(a_w,step, temperature,Lc);
  double energy = 0.0;
  for(int i=0; i<omega_arr.size(); i++) energy += omega_arr.at(i)*a_w.at(i)/(exp(omega_arr.at(i)-mu)/temperature+1) ;
  // cout << "temperature = " << temperature << ", mu = " << mu << endl;
  return energy;
}

double dos_and_filled_E_disconnected(MatrixXd rs, int li, int Lc, MatrixXd& dos_l, MatrixXd& dos_r) //li ~ lattice_index
{
  double step = 0.02;
  double energy = 0.0;
  double no_of_electrons = 0.0;
  double omega;
  dos_l = MatrixXd::Zero((omega_U-omega_L)/step+1, 2);
  dos_r = MatrixXd::Zero((omega_U-omega_L)/step+1, 2);
  int count_for_omega = 0;
  // begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());

  MatrixXd selected_rs_l = select_block(rs,li,Lc);
  MatrixXd selected_rs_r = select_block(rs, (li+Lc)%size, Lc);
  MatrixXcd Hl = cluster_h0(Lc)-U/2*cluster_sigmaz(Lc,selected_rs_l);
  MatrixXcd Hr = cluster_h0(Lc)-U/2*cluster_sigmaz(Lc,selected_rs_r);

  for(omega = omega_L; omega <= omega_U; omega += step)
  {
    double energy = 0.0;
    double no_of_electrons = 0.0;
  
    if(no_of_electrons < Lc)
    {
      MatrixXcd Z = cd(omega,eta)*MatrixXcd::Identity(Hl.rows(),Hl.rows());
      MatrixXcd tau = MatrixXcd::Zero(Hl.rows(),Hr.cols());
      tau(Lc-1,0) = tau(2*Lc-1, Lc) = -t;
      tau(0,Lc-1) = tau(Lc, 2*Lc-1) = -t; //outer PBC
      MatrixXcd Gl = invert(Z-Hl);
      MatrixXcd Gr = invert(Z-Hr); 

      double N_wl = -1/M_PI*Gl.trace().imag();
      energy += omega*N_wl*step;
      no_of_electrons += N_wl*step;
      dos_l(count_for_omega,0) = dos_r(count_for_omega,0) = omega;
      dos_l(count_for_omega,1) = N_wl;
      dos_r(count_for_omega,1) = -1/M_PI*Gr.trace().imag();
      count_for_omega++;
    }
    else 
      break;
  }
    // end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
    // cout << setw(5) << omega << " " << setw(5) << (end_ms.count()-begin_ms.count()) << endl;
    return energy;
}

double dos_and_filled_E_connected(MatrixXd rs, int li, int Lc, MatrixXd& dos_l, MatrixXd& dos_r) //li ~ lattice_index
{
  double step = 0.04;
  double energy = 0.0;
  double no_of_electrons = 0.0;
  double omega;
  dos_l = MatrixXd::Zero((omega_U-omega_L)/step+1, 2);
  dos_r = MatrixXd::Zero((omega_U-omega_L)/step+1, 2);
  int count_for_omega = 0;
  // begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());

  MatrixXd selected_rs_l = select_block(rs,li,Lc);
  MatrixXd selected_rs_r = select_block(rs, (li+Lc)%size, Lc);
  MatrixXcd Hl = cluster_h0(Lc)-U/2*cluster_sigmaz(Lc,selected_rs_l);
  MatrixXcd Hr = cluster_h0(Lc)-U/2*cluster_sigmaz(Lc,selected_rs_r);

  for(omega = omega_L; omega <= omega_U; omega += step)
  {
    double energy = 0.0;
    double no_of_electrons = 0.0;
  
    if(no_of_electrons < Lc)
    {
      MatrixXcd Z = cd(omega,eta)*MatrixXcd::Identity(Hl.rows(),Hl.rows());
      MatrixXcd tau = MatrixXcd::Zero(Hl.rows(),Hr.cols());
      tau(Lc-1,0) = tau(2*Lc-1, Lc) = -t;
      tau(0,Lc-1) = tau(Lc, 2*Lc-1) = -t; //outer PBC
      MatrixXcd Gl = invert(Z-Hl-tau*invert(Z-Hr)*tau.adjoint());
      MatrixXcd G_d_r = invert(Z-Hr); 
      MatrixXcd Gr = G_d_r + G_d_r*tau.adjoint()*Gl*tau*G_d_r;

      double N_wl = -1/M_PI*Gl.trace().imag();
      energy += omega*N_wl*step;
      no_of_electrons += N_wl*step;
      dos_l(count_for_omega,0) = dos_r(count_for_omega,0) = omega;
      dos_l(count_for_omega,1) = N_wl;
      dos_r(count_for_omega,1) = -1/M_PI*Gr.trace().imag();
      count_for_omega++;
    }
    else break;
  }


  // end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
  // cout << setw(5) << omega << " " << setw(5) << (end_ms.count()-begin_ms.count()) << endl;
  return energy;
}

#endif
