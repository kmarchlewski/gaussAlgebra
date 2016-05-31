#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <omp.h>

#define DEBUG(x)
 
typedef int ULONG;

const double pi = atan(1)*4;
 
ULONG choose(ULONG n, ULONG k)
{
	ULONG r = 1, d = n - k;
 
	/* choose the smaller of k and n - k */
	if (d > k) { k = d; d = n - k; }
 
	while (n > k) {
		if (r >= INT_MAX / n) return 0; /* overflown */
		r *= n--;
 
		/* divide (n - k)! as soon as we can to delay overflows */
		while (d > 1 && !(r % d)) r /= d--;
	}
	return r;
}
 
 
 
class funmat {
public:
	double * tab;
	int D, o, nf;
	int d, n;
	funmat(double * tab_, int D_, int o_, int nf_): tab(tab_), D(D_), o(o_), nf(nf_), d(0), n(0) {};
	void set_d(int d_) { d=d_; }
	void set_n(int n_) { n=n_; }
	double sigma() { return tab[d+n*(3+o)*D]; }
	double mu() { return tab[d+D+n*(3+o)*D]; }
	double c(int i) { 
		if (i>o) return 0;
		else return tab[d+D*(i+2)+n*(3+o)*D];
	}
	void set_sigma(double v) { tab[d+n*(3+o)*D] = v; }
	void set_mu(double v) { tab[d+D+n*(3+o)*D] = v; }
	void set_c(int i, double v) { 
		if (i>o) {DEBUG(printf("------------ Something is wrong! Beware! -----------\n");)}
		else tab[d+D*(i+2)+n*(3+o)*D] = v;
	}
	void add_c(int i, double v) { 
		if (i>o) {DEBUG(printf("------------ Something is wrong! Beware! -----------\n");)}
		else tab[d+D*(i+2)+n*(3+o)*D] += v;
	}
};

double normal_central_moment(const double & sig, const int & n) {
	if (n % 2) return 0; //odd moments are equal to 0
	double ret = 1;
	for (int i=1; i < n; i+=2) {
		ret *= i * sig;
	}
	return ret;
}

void conv_poly_mult(funmat & F, funmat & G, funmat & FG, double & wf, double & wg, double & sig) {
	double wgk;
	double wfp;
	wgk=1.;
	double mult = sqrt(2*pi*sig);
	for (int k=0; k<=F.o; k++) {
		wfp = 1;
		for (int p=0; p<=G.o; p++) {
			double sum =0;
			for (int i=0; i<=F.o-k; i++) {
				for (int j=0; j<=G.o-p; j++) {
					double ret = 1;
					ret *= choose(i+k,k);
					ret *= choose(j+p,p);
					ret *= F.c(i+k);
					ret *= G.c(j+p);
					if (j % 2) ret = -ret;
					ret *= normal_central_moment(sig, i+j);
					sum += ret;
				}
			}
			sum *= wfp * wgk;
			sum *= mult;
			DEBUG(printf("c[%d] += %f\n",p+k,sum);)
			FG.add_c(p+k, sum);
			wfp *= wf;
		}
		wgk *= wg;
	}		
}
 

void mult_poly_mult(funmat & F, funmat & G, funmat & FG, double  wf, double  wg, double & sig) {
	double wfj;
	double wgi;
	double mult;
	mult = F.mu() - G.mu();
	wf *=  mult;
	wg *= -mult;
	mult = exp(-(mult*mult)/(2*(F.sigma() + G.sigma())));
	for (int k=0; k<=F.o; k++) {
		for (int p=0; p<=G.o; p++) {
			double sum =0;
			wgi = 1;
			for (int i=0; i<=F.o-k; i++) {
				wfj=1;
				for (int j=0; j<=G.o-p; j++) {
					double ret = 1;
					ret *= choose(i+k,k);
					ret *= choose(j+p,p);
					ret *= F.c(i+k);
					ret *= G.c(j+p);
					ret *= wfj * wgi;
					DEBUG(printf("(%f,%f)+ %f\n",wfj,wgi,ret);)
					sum += ret;
					wfj *= wf;
				}
				wgi *= wg;
			}
			sum *= mult;
			DEBUG(printf("c[%d] += %f\n",p+k,sum);)
			FG.add_c(p+k, sum);
		}
	}		
}
 


extern "C" {

void fun(int * sizes, double * v) {
	for (int i=0;i<sizes[0];i++) {
		if (isinf(v[i])) {
			DEBUG(printf("Inf\n");)
		} else {
			DEBUG(printf("%f\n",v[i]);)
		}
	}
}


void conv(int * sizes, double * f, double * g, double * fg) {
	int D = sizes[0], fo = sizes[1], fn=sizes[2];
	int go = sizes[3], gn=sizes[4];
	int mult = sizes[5];
	funmat F(f,D,fo,fn);
	funmat G(g,D,go,gn);
	funmat FG(fg,D,fo+go,fn*gn);
	DEBUG(printf("Dimension: %d\n",D);)
	DEBUG(printf("f: polynomial order %d. number of functions %d\n",fo,fn);)
	DEBUG(printf("g: polynomial order %d. number of functions %d\n",go,gn);)
	DEBUG(if (mult) printf("Multiplication\n"); else printf("Convolution\n");)
	for (int ff = 0; ff < fn; ff++) {
		for (int gf = 0; gf < gn; gf++) {
			F.set_n(ff);
			G.set_n(gf);
			FG.set_n(ff + gf*fn);
			for (int d = 0; d < D; d++) {
				F.set_d(d);
				G.set_d(d);
				FG.set_d(d);
				double wf,wg,sig;
				if (isinf(F.sigma())) {
					wf  = 0.;
					wg  = 1.;
					sig = G.sigma();
				} else if (isinf(G.sigma())) {
					wf  = 1.;
					wg  = 0.;
					sig = F.sigma();
				} else {
					sig = F.sigma() + G.sigma();
					wf  = G.sigma() / sig;
					wg  = F.sigma() / sig;
					sig = F.sigma() * wf;
				}
				DEBUG(printf("sigma_f=%f sigma_g=%f\n",F.sigma(),G.sigma());)
				if (mult) {
					FG.set_sigma(sig);
					FG.set_mu(wf*F.mu() + wg*G.mu());
					mult_poly_mult(F,G,FG, wf, wg, sig);
				} else {
					FG.set_sigma(F.sigma() + G.sigma());
					FG.set_mu(F.mu() + G.mu());
					conv_poly_mult(F,G,FG, wf, wg, sig);
				}
				DEBUG(printf("w_f=%f w_g=%f sigma_{fg}=%f sigma_sum=%f\n",wf,wg,sig,F.sigma() + G.sigma());)
			}
		}
	}
	
}

int threads_max;
  
void init_threads_max() {
  printf("Checking maximum number of threads.\n");
  threads_max = omp_get_max_threads();
  printf("threads_max: %d\n", threads_max);
}
  
  
void calc(int * sizes, double * f, double * x, double * ret) {
	int D = sizes[0], fo = sizes[1], fn=sizes[2];
	int xn = sizes[3];
	int threads_to_use;
	int fn_max, threads_fn;

	fn_max = 1000;
	threads_fn = (int) ceil(fn / ((double) fn_max));
	
	if (threads_fn > threads_max) {
	  threads_to_use = threads_max;
	} else {
	  threads_to_use = threads_fn;
	}

	omp_set_num_threads(threads_to_use);
	
	funmat F(f,D,fo,fn);
	DEBUG(printf("Dimension: %d\n",D);)
	DEBUG(printf("f: polynomial order %d. number of functions %d\n",fo,fn);)
	for (int i = 0; i < xn; i++) {
		double sum=0;
    #pragma omp parallel for reduction(+:sum) firstprivate(F)
	  for (int ff = 0; ff < fn; ff++) {
			double prod=1;
			F.set_n(ff);
			for (int d = 0; d < D; d++) {
				F.set_d(d);
				double px = x[i+d*xn];
				double pxk = 1;
				px -= F.mu();
        
				double ret = 0;
				for (int k = 0; k <= fo; k++) {
					ret += pxk * F.c(k);
					pxk *= px;
				}
				ret *= exp(-(px*px)/(2*F.sigma()));
				prod *= ret;
			}
			sum += prod;
		}
		ret[i] = sum;
    DEBUG(printf("%d:%f\n",i,ret[i]));
	}
}


}///// extern "C"
