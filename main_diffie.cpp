 

#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>
#include <gmp.h>

#include <cstring>

#include <cstdio>

#include <time.h>

#include <unistd.h>
unsigned some_pseudo_rnd_detector(unsigned entropy)
{
	static unsigned x1 = 0, x2 = 0, x3 = 0, x4 = 0, x5 = 0;
	int x;
	x = (22695477*x1 +  3*x5*x5 + 1) + entropy;
	x5 = x4;
	x4 = x3;
	x3 = x2;
	x2 = x1;
	x1 = x;
	return x;
}

//size - число десятичных знаков в числе (фактических)
void big_rand(mpz_t &res, const unsigned size)
{
	char *buffer;
	char tmp[10];
	unsigned rnd;
	buffer = new char[size+1];
	for(unsigned i = 0; i < size; )
	{
		sprintf(tmp, "%u", some_pseudo_rnd_detector(0));
		int len = strlen(tmp);
		for(unsigned j = 0; i < size && j < len; j++, i++)
		{
			buffer[i] = tmp[j];
		}
	}
	buffer[size] = 0;
	mpz_set_str(res, buffer, 10);
	delete[] buffer;
}

//size - число десятичных знаков в числе (фактических)
void gen_prime_rnd(mpz_t &res, const unsigned size)
{
	big_rand(res, size);
	mpz_nextprime(res, res);
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


struct pow_mod_pthread_arg
{
	mpz_t *g, *a, *p, *res;
};

static void * pow_mod_pthread(void *arg)
{
	pow_mod_pthread_arg *in = reinterpret_cast<pow_mod_pthread_arg*>(arg);
	mpz_powm(*in->res, *in->g, *in->a, *in->p);
	return (void*)0;
}

//факторизация числа, не более factor_depth чем элементов
void ifactor(std::vector<int> &fact, const mpz_t &a, const int factor_depth)
{
	int rem;
	mpz_t tmp, n;
	mpz_inits(tmp, n, 0);
	mpz_set(n, a);
	for (int i=2; i < factor_depth; i = i + 2)
	{
		if(i == 4)//for sequence 2 3 5 7 9 11, for only odd numbers
		i--;
		rem = mpz_fdiv_ui(n, i);		// rem = n % i
		if (rem == 0)
		{
			fact.push_back(i);
			while(rem == 0)
			{
				rem = mpz_fdiv_q_ui(tmp, n, i);	//rem = n % i, tmp = n / i;
				if(rem == 0)
					mpz_set(n, tmp);
			}
		}
	}
	#if DEBUG
	std::cout << "factors for ";
	std::cout.flush();
	mpz_out_str(stdout, 10, a);
	std::cout <<  ": ";
	for(int i = 0; i < fact.size(); i++)
	{
		std::cout << fact[i] << " ";
	}
	std::cout << "\n";
	#endif
	mpz_clears(tmp, n, 0);
}


void primitive_root (mpz_t &g, const mpz_t &p, const int greater, const bool odd, const int factor_depth) 
{
	std::vector<int> fact;
	mpz_t phi, tmp, q;
	mpz_inits(phi, tmp, q, 0);
	mpz_sub_ui(phi, p, 1);//phi = p - 1
	
	ifactor(fact, phi, factor_depth);

	for (int i = 2; i < factor_depth; i++)
	{
		bool ok = true;
		for(size_t j = 0 ; j < fact.size() && ok; j++)
		{										//(i^(phi / fact[j]) )%p
			mpz_fdiv_q_ui(q, phi, fact[j]);		//q = phi / fact[j]
			mpz_set_ui(tmp, i); 				//tmp = i;
			mpz_powm(q, tmp, q, p);				//q = tmp^q % p;
			ok &= mpz_cmp_ui(q, 1) != 0;		//ok &= q != 0;
		}
		if(ok && i > greater && i & 1 == odd)
		{
			mpz_set_ui(g, i);
			break;
		}
	}
	mpz_clears(phi, tmp, q, 0);
}

/*
void primitive_r1(mpz_t &g, const mpz_t &p)
{
	std::vector<int> roots;
	primitive_root(p, roots);
	for(int i = 0; i < roots.size();i++)
	{	//std::cout << roots[i] << "   ";
		if(roots[i] > 2 && roots[i] & 1 == 1)
		{
			mpz_set_ui(g, roots[i]);
			return;
		}
	}
	if(roots.size() > 0)
		mpz_set_ui(g, roots[0]);//видимо не будет этого случая, есть другие корни, в том числе сам модуль
}
*/

void conditional_diffie_hellman(int rand_size, int mod_size)
{
	mpz_t g, p, a, b, A, B, Ka, Kb;

	mpz_inits(a, b, g, p, A, B, Ka, Kb, 0);
	big_rand(a, rand_size);
	big_rand(b, rand_size);
	gen_prime_rnd(p, mod_size);
	primitive_root(g, p, 2, true, 20);
	std::cout << "p: ";	std::cout.flush();	mpz_out_str(stdout, 10, p);
	std::cout << "\ng: ";	std::cout.flush();	mpz_out_str(stdout, 10, g);
	std::cout << "\n\n";

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
	mpz_clears(a, b, g, p, A, B, Ka, Kb, 0);
}


void create_inverse_m(mpz_t **m, const mpz_t *p, const int size)
{
	mpz_t tmp;
	mpz_init(tmp);
	for(int i = 0; i < size; i++)
	{
		m[i] = new mpz_t[size];//Верхнетреугольная матрица, для простоты все вычислим
		for(int j = 0; j < size; j++)
		{
			mpz_init(m[i][j]);
			mpz_sub_ui(tmp, p[j], 2);
			mpz_powm(m[i][j], p[i], tmp, p[j]);//m[i][j] = p[i]^(p[j]-2) % p[j];
		}
	}
	mpz_clear(tmp);
}
void clear_inverse_m(mpz_t **m, const int size)
{
	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < i; j++)
		{
			mpz_clear(m[j][i]);
		}
		delete[] m[i];
	}
}

//преобразование в ПСС методом Гарнера
void rns2conditional(mpz_t &res, const mpz_t *a, const mpz_t *p, mpz_t **r, const int size)
{
	mpz_t *x;
	mpz_t tmp, prod;
	x = new mpz_t[size];
	mpz_inits(tmp, prod, 0);
	for(int i = 0; i < size; i++)
	{
		mpz_init(x[i]);
		mpz_set(x[i], a[i]);			//x[i] = a[i];
		for(int j = 0; j < i; j++)
		{
			mpz_sub(tmp, x[i], x[j]);	//x[i] = r[j][i]*(x[i] - x[j]);
			mpz_mul(x[i], r[j][i], tmp);
			mpz_sub(tmp, x[i], p[i]);	//tmp = x[i] - p[i];
		}
		mpz_fdiv_r(x[i], x[i], p[i]);
		if(mpz_cmp_ui(x[i], 0) < 0)		//if(x < 0)
			mpz_add(x[i], x[i], p[i]);	//x += p[i]
	}
	mpz_set(res, x[0]);//res = x[0];
	mpz_set_ui(tmp, 1);//tmp = 1;
	for(int i = 0; i < size - 1; i++)
	{
		mpz_mul(tmp, tmp, p[i]);//tmp = tmp*p[i]; 
		mpz_mul(prod, x[i+1], tmp);//res += x[i]*tmp;
		mpz_add(res, res, prod);
	}
	mpz_clears(tmp, prod, 0);
	delete[] x;
}

void rns_diffie_hellman(int rand_size, int mod_size, int size)
{
	mpz_t *g, *p, *a, *A, *B, *Ka;
	mpz_t **inv_p;
	mpz_t condG, P, condA, condb, condB, condKa, condKb;
	double alice_rand, alice_create, alice_calc1, alice_calc2, time_st;
	double bob_rand, bob_calc;
	p = new mpz_t[size];
	g = new mpz_t[size];
	a = new mpz_t[size];
	A = new mpz_t[size];
	B = new mpz_t[size];
	Ka = new mpz_t[size];
	inv_p = new mpz_t*[size];
	/*
	#pragma omp parallel num_threads(5) 
	{
		std::cout << "Hello World!\n";
	}*/
	//Alice
	time_st = omp_get_wtime();
	mpz_init(P);//	mpz_init(condG);	mpz_init(condA);	mpz_init(condKa);	mpz_init(condb);	mpz_init(condB);	mpz_init(condKb);
	mpz_inits(condG, condA, condKa, condb, condB, condKb, 0l);
	
	mpz_set_ui(P, 1);
	for(int i = 0; i < size;i++)	//Генерация простых случайных чисел, Алиса
	{
		mpz_inits(g[i], p[i], a[i], A[i], B[i], Ka[i], 0l);
		//mpz_init(g[i]);	mpz_init(p[i]);	mpz_init(a[i]);	mpz_init(A[i]);	mpz_init(B[i]);	mpz_init(Ka[i]);
		gen_prime_rnd(p[i], mod_size);//300ms
		big_rand(a[i], rand_size);
	}
	alice_rand = omp_get_wtime() - time_st;
	time_st = omp_get_wtime();
	for(int i = 0; i < size;i++)	//вычисление первообразных корней
	{
		primitive_root(g[i], p[i], 2+i, true, 40);//2+i отличие корней
		mpz_mul(P, P, p[i]);			//P *= p[i];
	}
	
	create_inverse_m(inv_p, p, size);		//
	rns2conditional(condG, g, p, inv_p, size);
	alice_create = omp_get_wtime() - time_st;
	time_st = omp_get_wtime();
	
	
	int num_threads = size;
	pow_mod_pthread_arg *th_args = new pow_mod_pthread_arg[num_threads];
	pthread_t *threads = new pthread_t[num_threads];
	for(int i = 0; i < num_threads;i++)
	{
		th_args[i].g = &g[i];
		th_args[i].a = &a[i];
		th_args[i].p = &p[i];
		th_args[i].res = new mpz_t[1];
		mpz_init(*th_args[i].res);
		pthread_create(&threads[i], 0, pow_mod_pthread, &th_args[i]);
		//pthread_join(threads[i], 0);
	}
	
	
	
	int thrmask[5], thrmasknum = size;
	for(int i = 0; i < 5;i++)
	{
		thrmask[i] = 1;
	}
	while(thrmasknum)
	{
	for(int i = 0; i < num_threads; i++)
	{
		//mpz_t *ret;
		if(thrmask[i] == 0)
			continue;
		//pthread_join(threads[i], reinterpret_cast<void**>(&ret));
		int joret = pthread_tryjoin_np(threads[i], 0);//reinterpret_cast<void**>(&ret));
		if(joret == 0)
		{
			thrmask[i] = 0;
			thrmasknum--;
			mpz_set(A[i], *th_args[i].res);
			mpz_clear(*th_args[i].res);
			//delete ret;
		}else
			usleep(1000);
	}}
	delete[] threads;
	delete[] th_args;
	/*
	for(int i = 0; i < size; i++)
	{
		
		mpz_init(res);
		mpz_powm(A[i], g[i], a[i], p[i]);
		mpz_set(A[i], res);
		mpz_clear(res);
	}*/
	alice_calc1 = omp_get_wtime() - time_st;
	rns2conditional(condA, A, p, inv_p, size);
	
	//std::cout << "Hello World!\n";
	//Bob
	time_st = omp_get_wtime();
	big_rand(condb, rand_size*size);
	bob_rand = omp_get_wtime() - time_st;
	mpz_powm(condB, condG, condb, P);
	mpz_powm(condKb, condA, condb, P);
	bob_calc = omp_get_wtime() - time_st - bob_rand;
	
	//Alice
	time_st = omp_get_wtime();

	#pragma omp parallel for num_threads(5)
	for(int i = 0; i < size; i++)
	{
		mpz_fdiv_r(B[i], condB, p[i]);
	}
	
	#pragma omp parallel for num_threads(5)
	for(int i = 0; i < size; i++)
	{
		mpz_powm(Ka[i], B[i], a[i], p[i]);
	}

	rns2conditional(condKa, Ka, p, inv_p, size);
	alice_calc2 = omp_get_wtime() - time_st;
	
	std::cout << "Alice\n*************************************************************\n";
	for(int i = 0; i < size; i++)
	{
		std::cout << "p[" << i << "]: ";	std::cout.flush();	mpz_out_str(stdout, 10, p[i]);
		std::cout << "\n\n";
	}
	for(int i = 0; i < size; i++)
	{
		std::cout << "g[" << i << "]: ";	std::cout.flush();	mpz_out_str(stdout, 10, g[i]);
		std::cout << "\n";
	}
	for(int i = 0; i < size; i++)
	{
		std::cout << "a[" << i << "]: ";	std::cout.flush();	mpz_out_str(stdout, 10, a[i]);
		std::cout << "\n\n";
	}
	for(int i = 0; i < size; i++)
	{
		std::cout << "A[" << i << "]: ";	std::cout.flush();	mpz_out_str(stdout, 10, A[i]);
		std::cout << "\n\n";
	}
	std::cout << "Transmission.....\nBob\n*************************************************************\n";
	std::cout << "\nP: ";		std::cout.flush();
	mpz_out_str(stdout, 10, P);
	std::cout << "\n\ncondG: ";	std::cout.flush();
	mpz_out_str(stdout, 10, condG);
	std::cout << "\n\ncondA: ";	std::cout.flush();
	mpz_out_str(stdout, 10, condA);
	std::cout << "\n\ncondB: ";	std::cout.flush();
	mpz_out_str(stdout, 10, condB);
	std::cout << "\n\ncondKb: ";	std::cout.flush();
	mpz_out_str(stdout, 10, condKb);
	
	std::cout << "\n\nAlice\n*************************************************************\n";
	for(int i = 0; i < size; i++)
	{
		std::cout << "Ka[" << i << "]: ";	std::cout.flush();	mpz_out_str(stdout, 10, Ka[i]);
		std::cout << "\n\n";
	}
	std::cout << "\n\ncondKa: ";	std::cout.flush();
	mpz_out_str(stdout, 10, condKa);
	
	
	std::cout << "\n\nTiming for Alice:\ncreate:\t" << double(alice_create) << "\ncalc1:\t" 
		<< double(alice_calc1)/CLOCKS_PER_SEC << "\ncalc2:\t" << double(alice_calc2)
		<< "\nrand:\t" << double(alice_rand)
		<< "\n\nTiming for Bob:\ncalc:\t" << double(bob_calc) << "\nrand:\t" << double(bob_rand) << "\n"; 

	
	clear_inverse_m(inv_p, size);
	for(int i = 0; i < size;i++)
	{
		mpz_clears(g[i], p[i], a[i], A[i], Ka[i], 0l);
	}
	mpz_clears(condG, condA, condKa, condb, condB, condKb, 0l);
	//mpz_clear(condG), mpz_clear(condA), mpz_clear(condKa), mpz_clear(condb), mpz_clear(condB), mpz_clear(condKb);
	mpz_clear(P);
	delete[] inv_p;
	delete[] A;
	delete[] a;
	delete[] g;
	delete[] p;
}




 
  
int main()
{
	//conditional_diffie_hellman(100, 300);
	rns_diffie_hellman(500, 800, 4);


#if 0
  mpz_t **r;
  mpz_t *p, *a, A, B;
  
  p = new mpz_t[5];
  a = new mpz_t[5];
  mpz_inits(p[0], p[1], p[2], p[3], p[4], 0);
  mpz_inits(a[0], a[1], a[2], a[3], a[4], 0);
  mpz_inits(A, B, 0);
  mpz_set_str(p[0], "2", 10);	mpz_set_str(p[1], "3", 10); mpz_set_str(p[2], "5", 10);	mpz_set_str(p[3], "7", 10);	mpz_set_str(p[4], "11", 10);
  
  r = new mpz_t*[5];
  create_inverse_m(r, p, 5);
  /*
  for(int i = 0; i < 5; i++)
  {
		for(int j = 0; j < 5; j++)
		{
			mpz_out_str(stdout, 10, r[i][j]);
			std::cout << "\t"; std::cout.flush();
		}
		/*for(int j = i; j < 5; j++)
			std::cout << "0\t";*/
		/*std::cout << "\n";
  }*/
	int ai[5] = {0, 0, 0, 0, 0};
	int ap[5] = {2, 3, 5, 7, 11};
	int S = 0;
	for(int i = 0; i < 2310; i++)
	{
		for(int j = 0; j < 5;j++)
		{
			ai[j] = (ai[j] + 1)%ap[j];
		}
	//mpz_set_str(a[0], "0", 10);	mpz_set_str(a[1], "1", 10); mpz_set_str(a[2], "0", 10);	mpz_set_str(a[3], "6", 10);	mpz_set_str(a[4], "10", 10);
		mpz_set_ui(a[0], ai[0]);	mpz_set_ui(a[1], ai[1]); mpz_set_ui(a[2], ai[2]);	mpz_set_ui(a[3], ai[3]);	mpz_set_ui(a[4], ai[4]);
		rns2conditional(A, a, p, r, 5);
		S += mpz_get_ui(A);
		mpz_out_str(stdout, 10, A);
		std::cout << "\t"; std::cout.flush();
	}
	std::cout << "\nS = " << S << "\n";		
			
  mpz_clears(A, B, 0);
  mpz_clears(a[0], a[1], a[2], a[3], a[4], 0);
  mpz_clears(p[0], p[1], p[2], p[3], p[4], 0);
  clear_inverse_m(r, 5);
  delete[] r;
  delete[] a;
  delete[] p;
  
  /*
  std::vector<int> root_list;
  mpz_set_str(g, "17077", 10);
  std::cout << "is prime " << mpz_probab_prime_p(g, 2) << "\n";
  mpz_nextprime(a, g);
  std::cout << "my prime is: ";
  mpz_out_str(stdout, 10, a);
  std::cout.flush(); std::cout << "\n";
  
  mpz_set_str(p, "113", 10);
  primitive_root(p, root_list);
  std::cout << "roots: "; std::cout.flush();
  for( int i = 0; i < root_list.size(); i++)
  {
	std::cout << root_list[i] << "  ";
  }
  mpz_clears(g, p, a, b, A, B, Ka, Kb, 0);
  
  std::cout << "\n\n";
  
  std::FILE *f = std::fopen("1.txt", "w");
  for(int i = 0; i < 256*256*256; i++)
  {
	char str[5];
	*(unsigned*)&str = some_pseudo_rnd_detector(0);
	str[4] = 0;
	std::fprintf(f, "%s", str);
  }
  std::fclose(f);
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
#endif
  return 0;
}


