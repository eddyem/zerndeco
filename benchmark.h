/*
 * benchmark.h
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

#pragma once
#ifndef __BENCHMARK_H__
#define __BENCHMARK_H__

#include "zernike.h"
#include "usefull_macros.h"

#ifndef WF_EPSILON
#define WF_EPSILON (1e-12)
#endif

typedef struct{
	size_t size;        // size of arrays
	double *zdata;      // Z-coordinate of wavefront
	polar *coordinates; // coordinates on WF (-1..1)
} wavefront;

// errors of wavefront restoration method
typedef struct{
	double WF_std;   // standard deviation for restored wavefront weighted by wavefront values
	double max_dWF;  // max{|WF - WF_0|/WF_0} if WF_0 > WF_EPSILON
	double *ZC_sum;  // data for std calculation: sum of differences
	double *ZC_sum2; //   and sum of their squares
	double *max_dZC; // max|ZC/ZC_0 - 1|
} WF_stat;

typedef enum{
	 Scatter_direct = 0      // ZdecomposeR
	,Scatter_LS              // LS_decompose
	,Scatter_QR              // QR_decompose
	,Scatter_grad            // directGradZdecomposeR
	,Scatter_LSgrad          // LS_gradZdecomposeR
	,WF_bench_size
} WF_bench_type;

// wavefront restoration methods (according to functions' names)
typedef struct{
	WF_stat stat[WF_bench_size];
	int Znum; // amount of Zernike coefficients
	int measurements; // amount of measurements: std=sqrt(ZC_sum2/measurements-(ZC_sum/measurements)^2)
} WF_benchmark;

void free_wavefront(wavefront **F);
double *calc_surface(int sz, double *Zidxs, int znum);
wavefront *calc_surfaceR(int sz, double *Zidxs, int znum);
wavefront *calc_surfaceRS(int sz, polar *points, double *Zidxs, int znum);
polar *calc_BTA_Hpoints(int *sz);
void free_bench(WF_benchmark **bench);
void do_benchmark(double *Zidxs0, int znum0);
#endif // __BENCHMARK_H__
