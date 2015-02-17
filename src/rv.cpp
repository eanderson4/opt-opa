#include "rv.h"

double ranvar::g(double fnormalized, double L, double p, double pc){
  double f=fnormalized;
  double prob=0;
  double Uc=L + pc*(1-L)/p;
  double a=-p*L/(1-L);
  double b=p/(1-L);

  if(f<=L) prob=0;
  else if (f<Uc) prob=a+b*f;
  else prob=1;

  return prob;
}
double ranvar::ginv(double eps, double L, double p, double pc){
  double LARGE=pow(10,4);
  double Uc=L + pc*(1-L)/p;
  double a=-p*L/(1-L);
  double b=p/(1-L);

  if(eps>pc) return LARGE;
  else if (eps==pc) return Uc;
  else if (eps>0 && eps <pc) return (eps-a)/b;
  else return 0;
}

double ranvar::simProb(double L, double p, double pc){
  double totalProb=0;
  double probInt=1/double(_N);

  double Uc=L + pc*(1-L)/p;
  double a=-p*L/(1-L);
  double b=p/(1-L);

  for(int i=0; i<getNum(); i++){
    double r=getValue(i);
    double g;

    if(r<=L)   g=0;
    else if (r<Uc){
      g= a + b*r;
    }
    else g=1;

    totalProb=totalProb+probInt*g;

  }
  return totalProb;
}

double ranvar::deriveMu(double L, double p, double pc, double mean, double stdv){
  double b=p/(1-L);
  double mu = mean;
  double sigma = stdv;
  double alpha_low=(L - mu)/sigma;
  double derive=b*(1-PHI(alpha_low));


  return derive;

}
double ranvar::d2Mu(double L, double p, double pc, double mean, double stdv){
  double a=-p*L/(1-L);
  double b=p/(1-L);

  double mu = mean;
  double sigma = stdv;

  double alpha_low=(L - mu)/sigma;
  
  double phi_low=phi(alpha_low);
  
  double derive;

  //Ignoring right side truncation
  derive = a*alpha_low*phi_low/sigma/sigma
    + b*mu*alpha_low*phi_low/sigma/sigma 
    + b*phi_low/sigma
    + b*alpha_low*alpha_low*phi_low/sigma;

  derive = b*phi_low/sigma;

  return derive;

}
double ranvar::d2Sigma(double L, double p, double pc, double mean, double stdv){
  double a=-p*L/(1-L);
  double b=p/(1-L);

  double mu = mean;
  double sigma = stdv;

  double alpha_low=(L - mu)/sigma;
  
  double phi_low=phi(alpha_low);
  
  double derive;

  //Ignoring right side truncation
  derive = (a+b*mu)*alpha_low*phi_low*(alpha_low*alpha_low - 2)/sigma/sigma
    + b*alpha_low*alpha_low*phi_low*(alpha_low*alpha_low-2)/sigma
    + b*alpha_low*alpha_low*phi_low/sigma;

  derive=b*phi_low/sigma*pow(alpha_low,2);
  
  return derive;

}
double ranvar::dSigmaMu(double L, double p, double pc, double mean, double stdv){

  double a=-p*L/(1-L);
  double b=p/(1-L);

  double mu = mean;
  double sigma = stdv;

  double alpha_low=(L - mu)/sigma;
  
  double phi_low=phi(alpha_low);
  
  double derive,dcheck;

  //Ignoring right side truncation
  dcheck = (a + b*mu)*phi_low*(alpha_low*alpha_low-1)/sigma/sigma
    + b*alpha_low*alpha_low*alpha_low*phi_low/sigma;

  derive=b*phi_low*alpha_low/sigma;
  cout<<"Error: "<<dcheck-derive<<endl;

  cout<<"alpha_low: "<<alpha_low<<endl;
  return derive;

}


double ranvar::deriveSigma(double L, double p, double pc, double mean, double stdv){
  double b=p/(1-L);
  double mu = mean;
  double sigma = stdv;
  double alpha_low=(L - mu)/sigma;

  double derive=b*phi(alpha_low);

  //  if (abs(deriveC - derive) >= _tol) throw analerr;

  return derive;

}

double ranvar::anaProb(double L, double p, double pc, double mean, double stdv){
  double Uc=L + pc*(1-L)/p;
  double a=-p*L/(1-L);
  double b=p/(1-L);

  double mu = mean;
  double sigma = stdv;

  if(sigma==0){
    double tp=a+b*mu;
    if(tp<0) tp=0;
    return tp;
  }

  double low = L;
  double high = Uc;
  
  double alpha_low=(low - mu)/sigma;
  double alpha_high=(high - mu)/sigma;
  double p_iL = PHI(alpha_low);
  double p_LU = PHI(alpha_high)-PHI(alpha_low);
  double ef_LU;
  if(p_LU==0) ef_LU = 0; 
  else ef_LU = mu+sigma*(( phi(alpha_low)- phi(alpha_high))/( PHI(alpha_high)- PHI(alpha_low)));

  double p_Ui = 1- PHI(alpha_high);
  double t_p = a*p_LU + b*ef_LU*p_LU + p_Ui;

  if (p_Ui >_tol){
    //    cout<<"Probability f>=Uc: "<<p_Ui<<endl;
    //throw analerr;
  }

  if(1==0){
  cout<<"\nAnalytic Prob"<<endl;
  cout<<"\nmu: "<<mu<<endl;
  cout<<"sigma: "<<sigma<<endl;
  cout<<"L: "<<L<<endl;
  cout<<"Uc: "<<Uc<<endl;
  cout<<"alpha_l: "<<alpha_low<<endl;
  cout<<"alpha_h: "<<alpha_high<<endl;
  cout<<"\nProb f in [-inf,L]: "<<p_iL<<endl;
  cout<<"Prob f in [L,Uc]: "<<p_LU<<endl;
  cout<<"Prob f in [Uc,inf]: "<<p_Ui<<endl;
  cout<<"\nE[f|f in [L,Uc]]: "<<ef_LU<<endl;
  cout<<"\na*Prob f in [L,Uc]: "<<a*p_LU<<endl;
  cout<<"b*E[f|f in [L,Uc]*Prob f in [L,Uc]: "<<b*ef_LU*p_LU<<endl;
  cout<<"1*Prob f in [Uc,inf]: "<<p_Ui<<endl;
  cout<<"\nTotal: "<<t_p<<endl;
  }

  return t_p;
}

double ranvar::sampleMean(){
  double mean=0;
  for(int i=0;i<getNum();i++){
    mean=mean+getValue(i);    
  }
  mean=mean/getNum();
  
  return mean;
}

double ranvar::sampleStDv(){
  double mean=sampleMean();

  double s2=0;
  for(int i=0;i<getNum();i++){
    s2=s2+pow((getValue(i)-mean),2);
  }
  s2=s2/(getNum()-1);
  s2=sqrt(s2);
  
  return s2;
}

void ranvar::createRV(int N,double mean, double stdv){
  _N=N;_mean=mean;_stdv=stdv;
  //Draw 0-1
  double zero_one = (rand() % RAND_MAX);
  zero_one = zero_one/RAND_MAX;

  //Gaussian
  std::default_random_engine generator;
  generator.seed(_seed);
  std::normal_distribution<double> distribution(mean,stdv);

  for(int i =0; i< N; i++){
    random_variable.push_back(distribution(generator));
  }
  
}

double ranvar::phi(double x){ 
  if(x >= std::numeric_limits<double>::max()) return 0;
  else return 1/sqrt(2*M_PI)*exp(-.5*pow(x,2)); 
}

double ranvar::PHI(double x){
  if(x >= std::numeric_limits<double>::max()) return 1;
  if(x <= std::numeric_limits<double>::lowest()) return 0;
  else   return 0.5 * erfc(-x * M_SQRT1_2);
}

double ranvar::testPhi()
{
  //CREDIT TO http://www.johndcook.com/cpp_phi.html
  // Select a few input values
  double x[] =
    {
      -3,
      -1,
      0.0,
      0.5,
      2.1
    };
  
  // Output computed by Mathematica
  // y = Phi[x]
  double y[] =
    {
      0.00134989803163,
      0.158655253931,
      0.5,
      0.691462461274,
      0.982135579437
    };
  
  int numTests = sizeof(x)/sizeof(double);
  
  double maxError = 0.0;
  for (int i = 0; i < numTests; ++i)
    {
      double error = fabs(y[i] - PHI(x[i]));
      if (error > maxError)
	maxError = error;
    }
  
  std::cout << "Maximum error: " << maxError << "\n";
  return maxError;
} 


double ranvar::RationalApproximation(double t)
{
  // Abramowitz and Stegun formula 26.2.23.
  // The absolute value of the error should be less than 4.5 e-4.
  double c[] = {2.515517, 0.802853, 0.010328};
  double d[] = {1.432788, 0.189269, 0.001308};
  return t - ((c[2]*t + c[1])*t + c[0]) / 
    (((d[2]*t + d[1])*t + d[0])*t + 1.0);
}

double ranvar::PHIInverse(double p)
{
  //http://www.johndcook.com/normal_cdf_inverse.html
  if (p <= 0.0 || p >= 1.0)
    {
      std::stringstream os;
      os << "Invalid input argument (" << p 
	 << "); must be larger than 0 but less than 1.";
      throw std::invalid_argument( os.str() );
    }

  // See article above for explanation of this section.
  if (p < 0.5)
    {
      // F^-1(p) = - G^-1(p)
      return -RationalApproximation( sqrt(-2.0*log(p)) );
    }
  else
    {
      // F^-1(p) = G^-1(1-p)
      return RationalApproximation( sqrt(-2.0*log(1-p)) );
    }
}


void ranvar::demoInverse()
{
  std::cout << "\nShow that the NormalCDFInverse function is accurate at \n"
	    << "0.05, 0.15, 0.25, ..., 0.95 and at a few extreme values.\n\n";

    double p[] =
      {
        0.0000001,
        0.00001,
        0.001,
        0.05,
        0.15,
        0.25,
        0.35,
        0.45,
        0.55,
        0.65,
        0.75,
        0.85,
        0.95,
        0.999,
        0.99999,
        0.9999999
      };

    // Exact values computed by Mathematica.
    double exact[] =
      {
        -5.199337582187471,
        -4.264890793922602,
        -3.090232306167813,
        -1.6448536269514729,
        -1.0364333894937896,
        -0.6744897501960817,
        -0.38532046640756773,
        -0.12566134685507402,
	0.12566134685507402,
	0.38532046640756773,
	0.6744897501960817,
	1.0364333894937896,
	1.6448536269514729,
	3.090232306167813,
	4.264890793922602,
         5.199337582187471
      };

    double maxerror = 0.0;
    int numValues = sizeof(p)/sizeof(double);
    std::cout << "p, exact CDF inverse, computed CDF inverse, diff\n\n";
    std::cout << std::setprecision(7);
    for (int i = 0; i < numValues; ++i)
      {
        double computed = PHIInverse(p[i]);
        double error = exact[i] - computed;
	std::cout << p[i] << ", " << exact[i] << ", " 
                  << computed << ", " << error << "\n";
        if (fabs(error) > maxerror)
	  maxerror = fabs(error);
      }

    std::cout << "\nMaximum error: " << maxerror << "\n\n";
}


ostream& operator<<(ostream& os, const ranvar& rv)
{

  int nE=rv.random_variable.size();

  for(int i =0; i<nE; i++){
    os<<rv.random_variable[i]<<endl;
  }
  return os;
}
