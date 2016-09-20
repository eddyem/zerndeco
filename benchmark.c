/*
 * benchmark.c - test of different methods
 *
 * Copyright 2016 Edward V. Emelianov <eddy@sao.ru, edward.emelianoff@gmail.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */
#include <stdio.h>
#include <math.h>
#include <time.h> // for srand48(time(NULL))
#include "benchmark.h"
#include "usefull_macros.h"
#include "cmdlnopts.h"

void free_wavefront(wavefront **F){
	FREE((*F)->zdata);
	FREE((*F)->coordinates);
	FREE(*F);
}
static void free_stat(WF_stat *st){
	if(!st) return;
	FREE(st->ZC_sum);
	FREE(st->ZC_sum2);
	FREE(st->max_dZC);
}

void free_bench(WF_benchmark **bench){
	if(!bench || !*bench) return;
	WF_benchmark *b = *bench;
	for(int i = 0; i < WF_bench_size; ++i)
		free_stat(&(b->stat[i]));
	FREE(*bench);
}

/**
 * calculate Z-coordinate of surface with rectangle nodes by given Zernike coefficients
 * @param sz     - image size
 * @param Zidxs  - Zernike coefficients
 * @param znum   - amount of Zidx
 */
double *calc_surface(int sz, double *Zidxs, int znum){
	if(!Zidxs || znum < 1) return NULL;
	return Zcompose(znum, Zidxs, sz, sz);
}

/**
 * calculate surface with rectangle nodes by given Zernike coefficients using `R`-functions
 * @param sz     - image size
 * @param Zidxs  - Zernike coefficients
 * @param znum   - amount of Zidx
 */
wavefront *calc_surfaceR(int sz, double *Zidxs, int znum){
	int G2 = sz * sz;
	if(!Zidxs || znum < 1) return NULL;
	wavefront *F = MALLOC(wavefront, 1);
	F->size = G2;
	F->coordinates = MALLOC(polar, G2);
	polar *rptr = F->coordinates;
	double Stp = 2./((double)sz - 1.);
	for(int j = 0; j < sz; ++j){
		double y = ((double) j) * Stp - 1.;
		for(int i = 0; i < sz; ++i, ++rptr){
			double x = ((double) i)* Stp - 1.;
			double R2 = x*x + y*y;
			rptr->r = sqrt(R2);
			rptr->theta = atan2(y, x);
		}
	}
	F->zdata = ZcomposeR(znum, Zidxs, G2, F->coordinates);
	return F;
}

/**
 * calculate surface with scattered nodes by given Zernike coefficients using `R`-functions
 * @param sz     - amount of points
 * @param points - given nodes coordinates
 * @param Zidxs  - Zernike coefficients
 * @param znum   - amount of Zidx
 */
wavefront *calc_surfaceRS(int sz, polar *points, double *Zidxs, int znum){
	if(!Zidxs || !points || znum < 1 || sz < 1) ERRX(_("Bad arguments of calc_surfaceRS:%s"));
	wavefront *F = MALLOC(wavefront, 1);
	F->size = sz;
	F->coordinates = MALLOC(polar, sz);
	memcpy(F->coordinates, points, sizeof(polar)*sz);
	F->zdata = ZcomposeR(znum, Zidxs, sz, F->coordinates);
	return F;
}

/**
 * calculate coordinates of wavefront for BTA Hartmann mask
 */
polar *calc_BTA_Hpoints(int *sz){
	// radius on Hartmann mask (in mm) - for normalize coordinates @ [-1, 1]
	double rmax = MIR_R * (1. - HARTMANN_Z/FOCAL_R);
	double r[] = {175., 247., 295., 340., 379., 414., 448., 478.}; // radii of rings
	polar *pts = MALLOC(polar, 258), *pt = pts;
	const double ray_step = M_PI/16., theta0 = M_PI/32; // angle between rays & first ray theta coordinate
	for(int ring = 0; ring < 8; ++ring){ // main rings
		for(int ray = 0; ray < 32; ++ray, ++pt){
			pt->r = r[ring] / rmax;
			pt->theta = ray * ray_step + theta0;
		}
	}
	// and markers: ring 2 ray 24.5 and ring 7 ray 29.5
	pt->r = r[2] / rmax; pt->theta = 24.5 * ray_step + theta0;
	++pt;
	pt->r = r[7] / rmax; pt->theta = 29.5 * ray_step + theta0;
	if(sz) *sz = 258;
	return pts;
}

static char *bench_names[] = {
	N_("Direct decomposition (ZdecomposeR)"),
	N_("Less square decomposition (LS_decompose)"),
	N_("QR-decomposition (QR_decompose)"),
	N_("Direct gradient decomposition (directGradZdecomposeR)"),
	N_("LS gradient decomposition (LS_gradZdecomposeR)")
};

/**
 * calculate difference of two wavefronts
 * @param WF1, WF2 (i) - wavefronts
 * @param sz (i) - size of WF1 and WF2
 * @param std  (o) - RMS of difference
 * @param maxd (o) - maximal difference value
 * @param maxd (o) - maximal relative difference value
 */
void WF_diff(double *WF1, double *WF2, int sz, double *std, double *maxd, double *maxdr){
	int i;
	double _maxd = 0., _maxdr = 0., sum = 0., sum2 = 0.;
	for(i = 0; i < sz; ++i){
		double z = WF1[i];
		double diff = z - WF2[i];
		sum += diff;
		sum2 += diff*diff;
		diff = fabs(diff);
		if(_maxd < diff){
			_maxd = diff;
			z = fabs(z);
			if(z > WF_EPSILON){
				_maxdr = diff / z;
			//	diff /= z;
			//	if(_maxdr < diff) _maxdr = diff;
			}
		}
	}
	sum /= sz; sum2 /= sz;
	if(std) *std = sqrt(sum2 - sum*sum);
	if(maxd) *maxd = _maxd;
	if(maxdr) *maxdr = _maxdr;
}

/**
 * calculate errors and fill WF_stat fields for single measurement
 * @param stat    - output stat structure
 * @param Zidxs   - calculated ZC
 * @param Zidxs0  - original ZC
 * @param znum    - size of 'Zidxs'
 * @param WF0     - original wavefront
 * @param verbose - ==1 for messages in stdout
 */
static void fill_stat(WF_benchmark *bench, WF_bench_type type, double *Zidxs, double *Zidxs0, wavefront *WF0, int verbose){
	double *zdata0 = WF0->zdata;
	polar *P = WF0->coordinates;
	int sz = WF0->size;
	WF_stat *stat = &bench->stat[type];
	int znum = bench->Znum;
	// test errors by wavefront
	double *zdata = ZcomposeR(znum, Zidxs, sz, P);
	//for(i = 0; i < znum; ++i) printf("ZC_%d\t%10.3f\t%10.3f\n", i, Zidxs0[i], Zidxs[i]);
	//for(i = 0; i < sz; ++i) printf("Point_%d\t%10.3f\t%10.3f\n", i, zdata0[i], zdata[i]);
	double maxdr;
	WF_diff(zdata0, zdata, sz, &(stat->WF_std), &(stat->max_dWF), &maxdr);
	stat->ZC_sum = MALLOC(double, znum);
	stat->ZC_sum2 = MALLOC(double, znum);
	stat->max_dZC = MALLOC(double, znum);
	if(verbose){
		printf("%s\n", bench_names[type]);
		printf("Wavefront error: STD=%g, |max err| = %.3f, |max err rel|=%.3f%%; Zernike:\n#\t   val0   \t    val   \tdiff%%\n", stat->WF_std, stat->max_dWF, maxdr);
	}
	for(int p = 0; p < znum; ++p){
		double diff = Zidxs0[p] - Zidxs[p];
		stat->ZC_sum[p] = diff;
		stat->ZC_sum2[p] = diff*diff;
		if(fabs(Zidxs0[p]) > Z_prec) stat->max_dZC[p] = fabs(diff / Zidxs0[p]);
		if(verbose){
			printf("%d\t%10f\t%10f\t", p, Zidxs0[p], Zidxs[p]);
			if(fabs(Zidxs0[p]) > Z_prec) printf("%4.1f\n", stat->max_dZC[p]*100.);
			else printf("  -  \n");
		}
	}
}

/**
 * Process a benchmark step with different methods of wavefront restoration
 * @param Zidxs   - Zernike coefficients
 * @param znum    - length of Zidxs
 * @param verbose - ==1 for detailed information on stdout
 */
static WF_benchmark *benchmark_step(double *Zidxs0, int znum, int verbose){
	#define MESG(...) do{if(verbose){printf(__VA_ARGS__); printf("\n");}}while(0)
	int sz; //G2 = G.imagesize * G.imagesize;
	wavefront *WF;
	polar *P;
	double *Zidxs, *zdata;
	int N, M, Zsz;
	convert_Zidx(znum, &N, &M); // calculate max power of Zc
	WF_benchmark *bench = MALLOC(WF_benchmark, 1);
	bench->Znum = znum;
	MESG(_("Direct decomposition of given wavefront"));
	MESG(_("1. Scattered points for BTA hartmannogram"));
	P = calc_BTA_Hpoints(&sz);
	WF = calc_surfaceRS(sz, P, Zidxs0, znum);
	zdata = WF->zdata;
	// ZdecomposeR
	Zidxs = ZdecomposeR(N, sz, WF->coordinates, zdata, &Zsz, NULL);
	fill_stat(bench, Scatter_direct, Zidxs, Zidxs0, WF, verbose);
	FREE(Zidxs);
	// LS_decompose
	Zidxs = LS_decompose(N, sz, WF->coordinates, zdata, &Zsz, NULL);
	fill_stat(bench, Scatter_LS, Zidxs, Zidxs0, WF, verbose);
	FREE(Zidxs);
	// QR_decompose
	Zidxs = QR_decompose(N, sz, WF->coordinates, zdata, &Zsz, NULL);
	fill_stat(bench, Scatter_QR, Zidxs, Zidxs0, WF, verbose);
	FREE(Zidxs);
	free_wavefront(&WF);
	Zidxs0[0] = 0.; // decomposition by gradients cannot know anything about zero coefficient
	WF = calc_surfaceRS(sz, P, Zidxs0, znum); // recalculate new wavefront
	// directGradZdecomposeR
	point *grads = directGradZcomposeR(znum, Zidxs0, sz, WF->coordinates);
	Zidxs = directGradZdecomposeR(N, sz, WF->coordinates, grads, &Zsz, NULL);
	fill_stat(bench, Scatter_grad, Zidxs, Zidxs0, WF, verbose);
	FREE(Zidxs);
	// LS_gradZdecomposeR
	Zidxs = LS_gradZdecomposeR(N, sz, WF->coordinates, grads, &Zsz, NULL);
	fill_stat(bench, Scatter_LSgrad, Zidxs, Zidxs0, WF, verbose);
	FREE(Zidxs);
	// gradZdecomposeR
	Zidxs = gradZdecomposeR(N, sz, WF->coordinates, grads, &Zsz, NULL);
	fill_stat(bench, Scatter_Zhao, Zidxs, Zidxs0, WF, verbose);
	FREE(grads);
	FREE(P);
	free_wavefront(&WF);
	/*WF = calc_surface(G.imagesize, Zidxs0, znum);
	P = WF->coordinates;
	zdata = WF->zdata;
	grads = gradZcomposeR(znum, Zidxs0, G2, P);
	Zidxs = ZdecomposeR(N, G2, WF->coordinates;, zdata, &Zsz, &lastIdx);
	FREE(grads);
	*/
	#undef MESG
	return bench;
}

/**
 * Update values stored in `total` by values of `current`
 */
static void update_stat_values(WF_benchmark *total, WF_benchmark *current){
	++total->measurements;
	int znum = total->Znum;
	for(int i = 0; i < WF_bench_size; ++i){
		WF_stat *out = &total->stat[i], *in = &current->stat[i];
		if(!out->ZC_sum || !in->ZC_sum) continue; // given field is empty - didn't run benchmark
		out->WF_std += in->WF_std;
		if(in->max_dWF > out->max_dWF)
			out->max_dWF = in->max_dWF;
		for(int j = 0; j < znum; ++j){
			out->ZC_sum[j] += in->ZC_sum[j];
			out->ZC_sum2[j] += in->ZC_sum2[j];
			if(out->max_dZC[i] < in->max_dZC[i])
				out->max_dZC[i] = in->max_dZC[i];
		}
	}
}

/**
 * Process a benchmark with different methods of wavefront restoration
 * @param Zidxs - user defined Zernike coefficients or NULL
 * @param znum  - length of Zidxs
 */
void do_benchmark(double *Zidxs0, int znum0){
	int znum = G.bench_maxP;
	double *Zidxs = MALLOC(double, znum);
	srand48(time(NULL));
	inline void gen_zidxs(){ // generate random Z coefficients in interval (-10..10)/p
		for(int i = 0; i <  znum; ++i)
			Zidxs[i] = (20.*drand48() - 10.) / ((double)i+1);
	}
	WF_benchmark *total = NULL, *current = NULL;
	if(Zidxs0){ // do only one iteration - for given coefficients
		current = benchmark_step(Zidxs0, znum0, 1);
		free_bench(&current);
		return;
	}
	if(G.bench_iter < 2) ERRX(_("Need at least two iterations for benchmark"));
	gen_zidxs();
	total = benchmark_step(Zidxs, znum, 0);
	total->measurements = 1;
	for(int i = 1; i < G.bench_iter; ++i){
		gen_zidxs();
		current = benchmark_step(Zidxs, znum, 0);
		update_stat_values(total, current);
		free_bench(&current);
	}
	// print stat values
	double Nm = total->measurements;
	for(int i = 0; i < WF_bench_size; ++i){
		WF_stat *stat = &total->stat[i];
		if(!stat->ZC_sum) continue;
		printf("%s:\n<STD> = %g\n", bench_names[i], stat->WF_std/Nm);
		printf("max_dWF=%g%%\n", stat->max_dWF);
		printf("     p(n,m)\t   Z_std   \tmax_dZC(%%) \tsum2\tsum\n");
		for(int i = 0; i < znum; ++i){
			double sum2 = stat->ZC_sum2[i] / Nm;
			double sum = stat->ZC_sum[i] / Nm;
			int N,M;
			convert_Zidx(i, &N, &M);
			printf("%3d (%2d, %2d)\t%10g\t%10g\t%10g\t%10g\n", i, N, M, sqrt(sum2 - sum*sum),
				stat->max_dZC[i]*100., sum2, sum);
		}
	}
	free_bench(&total);
	FREE(Zidxs);
}

