 

#include <iostream>
#include <cmath>
#include <vector>
#include <gmp.h>

void big_rand(mpz_t &res)
{
 mpz_set_str(res, "1040793219466439908192524032736408553861526224726670480531911235040360805967336029801223944173232418"
		  "4842421613954281007791383566248323464908139906605677320762924129509389220345773183349661583550472959"
		  "4205476898112116936771475484788669625013844382602917323488853111608285384165850282556046662248318909"
		  "188018470682222031405210266984354887329580288780508697361869007147207", 10); 
}

//A = g^a mod p 
unsigned pow_mod_std(const unsigned g, const unsigned a, const unsigned p)
{
  return unsigned(pow((g), a)) % p;
}

void pow_mod_stdz(mpz_t &res, const mpz_t &g, const mpz_t &a, const mpz_t &p)
{
  mpz_powm(res, g, a, p);
}






//что если ограничить фактаризацию? и выбрать из нескольких корней?
#define MAX_FACT_SIZE 200
long generator (const mpz_t &p) 
{
	std::vector<int> fact;
	int rem;
	mpz_t phi, n, tmp;
	mpz_init_set(phi, p);
	mpz_sub_ui(phi, phi, 1);//phi = p - 1
	mpz_init_set(n, phi);
	mpz_init(tmp);
	
	std::cout << "n = ";
	mpz_out_str(stdout, 10, n);
	   
	for (int i=2; i < MAX_FACT_SIZE; ++i)
	{
	  rem = mpz_fdiv_ui(n, i);
	  /*mpz_out_str(stdout, 10, n);
	  std::cout << "  i" << i << " r" << rem << "\t\t";*/
	      
	  if (rem == 0)//n % i == 0)
	  {
	    fact.push_back (i);
	    //std::cout << "   " << i;
			//while (n % i == 0)
				//n /= i;
	    while(rem == 0)
	    {
	      rem = mpz_fdiv_q_ui(tmp, n, i);
	      if(rem == 0)
		mpz_set(n, tmp);
	      /*mpz_out_str(stdout, 10, n);
	      std::cout << "  i" << i << " r" << rem << "\t\t";
	      if(mpz_cmp_ui(n, 0) == 0)
		 break;*/
	      /*std::cout << "i" << i << " ";
	      mpz_out_str(stdout, 10, n);
	      std::cout << "\n = "; std::cout.flush();*/
	    }
	  }
	}
	/*if (n > 1)
		fact.push_back (n);*/
	
	std::cout << "fact ";
	for(int i = 0; i < fact.size(); i++)
	{
	  std::cout << fact[i] << " ";
	}
	std::cout << "\n";
	
	mpz_t q;
	mpz_init(q);
	for (int res=2; res < MAX_FACT_SIZE; ++res)
	{
	  bool ok = true;
	  for(size_t i=0; i<fact.size() && ok; ++i)
	  {
	    mpz_fdiv_q_ui(q, phi, fact[i]);//phi / fact[i]
	    mpz_set_ui(tmp, res); 
	    mpz_powm(q, tmp, q, p);
	    ok &= mpz_cmp_ui(q, 1) != 0;
	    //ok &= powmod (res, phi / fact[i], p) != 1; //res^phi / fact[i] %p
	  }
	  if(ok)return res;
	}
	return 0;
}


int main()
{
  mpz_t g, p, a, b, A, B, Ka, Kb;
  
  mpz_init_set_str(g, "5", 10);
  mpz_init_set_str(p, "1040793219466439908192524032736408553861526224726670480531911235040360805967336029801223944173232418"
    "4842421613954281007791383566248323464908139906605677320762924129509389220345773183349661583550472959"
    "4205476898112116936771475484788669625013844382602917323488853111608285384165850282556046662248318909"
    "18801847068222203140521026698435488732958028878050869736186900714720710555703168729087", 10);
  //mpz_init_set_str(a, "4", 10);
  mpz_init(a);
  mpz_init_set_str(b, "76464564", 10);
  mpz_inits(A, B, Ka, Kb, 0);
  big_rand(a);
  //mpz_set_ui(integ, 1232);
  //mpz_out_str(stdout, 10, integ);
  pow_mod_stdz(A, g, a, p);
  pow_mod_stdz(B, g, b, p);
  pow_mod_stdz(Ka, B, a, p);
  pow_mod_stdz(Kb, A, b, p);
  
  std::cout << "A = "; std::cout.flush();
  mpz_out_str(stdout, 10, A);
  std::cout << "\n\nB = "; std::cout.flush();
  mpz_out_str(stdout, 10, B);
  std::cout << "\n\nKa = ";std::cout.flush();
  mpz_out_str(stdout, 10, Ka);
  std::cout << "\n\nKb = "; std::cout.flush();
  mpz_out_str(stdout, 10, Kb);
  std::cout << std::endl;
  
  //mpz_set_str(p, "113", 10);
  std::cout << generator(p) << "\n";
  mpz_clears(g, p, a, b, A, B, Ka, Kb, 0);
  //mpz_clear(integ);*/
  
  /*unsigned g = 3;
  unsigned p = 11;
  
  unsigned a = 4;
  unsigned b = 7;
  
  unsigned A, B, Ka, Kb;
  
  A = pow_mod_std(g, a, p);
  B = pow_mod_std(g, b, p);
  Ka = pow_mod_std(B, a, p);
  Kb = pow_mod_std(A, b, p);
  std::cout << "A = " << A << "\t\tB = " << B << "\t\tKa = " << Ka << "\t\tKb = " << Kb << std::endl;*/
  
  return 0;
}


