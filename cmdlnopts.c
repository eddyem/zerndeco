/*
 * cmdlnopts.c - the only function that parse cmdln args and returns glob parameters
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
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "zernike.h"
#include "cmdlnopts.h"
#include "usefull_macros.h"

#define RAD 57.2957795130823
#define D2R(x) ((x) / RAD)
#define R2D(x) ((x) * RAD)

/*
 * here are global parameters initialisation
 */
int help;
glob_pars  G;

int rewrite_ifexists = 0, // rewrite existing files == 0 or 1
	verbose = 0; // each -v increments this value, e.g. -vvv sets it to 3
//            DEFAULTS
// default global parameters
glob_pars const Gdefault = {
	 .pixsize = 30.e-3
	,.outfile = "wavefront"
	,.imagesize = 1024
	,.bench_maxP = 21  // Nmax = 5
	,.bench_iter = 100
};

/*
 * Define command line options by filling structure:
 *	name	has_arg	flag	val		type		argptr			help
*/
myoption cmdlnopts[] = {
	// set 1 to param despite of its repeating number:
	{"help",	NO_ARGS,	NULL,	'h',	arg_int,	APTR(&help),		N_("show this help")},
	{"verbose",	NO_ARGS,	NULL,	'v',	arg_int,	APTR(&verbose),		N_("be more and more verbose")},
	{"rewrite",	NO_ARGS,	NULL,	'r',	arg_none,	APTR(&rewrite_ifexists),N_("rewrite existant files")},
	{"gradients",NEED_ARG,	NULL,	'g',	arg_string,	APTR(&G.gradname),	N_("file with pre-computed gradients")},
	{"pixsize",	NEED_ARG,	NULL,	'p',	arg_double,	APTR(&G.pixsize),	N_("CCD pixel size (if differs from FITS header)")},
	{"zern-gen",NO_ARGS,	NULL,	'z',	arg_none,	APTR(&G.zerngen),	N_("generate FITS file with wavefront by given Zernike coefficients, in benchmark mode - use this coefficients instead of random")},
	{"output",	NEED_ARG,	NULL,	'o',	arg_string,	APTR(&G.outfile),	N_("output file name (\"wavefront\" by default")},
	{"benchmark",NO_ARGS,	NULL,	'b',	arg_none,	APTR(&G.benchmark),	N_("benchmark: test accuracy of different methods")},
	{"imsize",	NEED_ARG,	NULL,	's',	arg_int,	APTR(&G.imagesize),	N_("size of output image (rectangular)")},
	{"bench-maxP",NEED_ARG,	NULL,	'm',	arg_int,	APTR(&G.bench_maxP),N_("maximum amount of Zernike coefficients for benchmark")},
	{"bench-iter",NEED_ARG,	NULL,	'i',	arg_int,	APTR(&G.bench_iter),N_("amount of iterations for benchmark (if -z is absent)")},
	{"zprec",	NEED_ARG,	NULL,	'Z',	arg_double,	APTR(&Z_prec),		N_("minimal value of non-zero Zernike coefficients")},
	end_option
};

/**
 * Parse command line options and return dynamically allocated structure
 * 		to global parameters
 * @param argc - copy of argc from main
 * @param argv - copy of argv from main
 * @return allocated structure with global parameters
 */
glob_pars *parse_args(int argc, char **argv){
	int i;
	void *ptr;
	ptr = memcpy(&G, &Gdefault, sizeof(G)); assert(ptr);
	// format of help: "Usage: progname [args]\n"
	change_helpstring(_("Usage: %s [args] <input files or Zernike coefficients>\n\n\tWhere args are:\n"));
	// parse arguments
	parseargs(&argc, &argv, cmdlnopts);
	if(help) showhelp(-1, cmdlnopts);
	if(argc > 0){
		G.rest_pars_num = argc;
		G.rest_pars = calloc(argc, sizeof(char*));
		for (i = 0; i < argc; i++)
			G.rest_pars[i] = strdup(argv[i]);
	}
	return &G;
}

