/*
 * cmdlnopts.h - comand line options for parceargs
 *
 * Copyright 2013 Edward V. Emelianoff <eddy@sao.ru>
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
#ifndef __CMDLNOPTS_H__
#define __CMDLNOPTS_H__

#include "parseargs.h"


typedef struct{
	char **rest_pars;    // names of input files or zernike coefficients
	int rest_pars_num;   // amount of files' names / coefficients
	char *gradname;      // file with pre-computed gradients
	double pixsize;      // CCD pixel size
	int zerngen;         // generate FITS file with Zernike surface
	char *outfile;       // output file name
	int imagesize;       // output image size (1024 by default)
	int benchmark;       // different methods benchmark
	int bench_maxP;      // max ZC number for behchmark (default - 20)
	int bench_iter;      // amount of iterations for benchmark (default - 100)
} glob_pars;


// default & global parameters
extern glob_pars const Gdefault;
extern int rewrite_ifexists, verbose;

glob_pars *parse_args(int argc, char **argv);
extern glob_pars  G;
#endif // __CMDLNOPTS_H__
