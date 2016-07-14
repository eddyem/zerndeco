/*
 * Z-BTA_test.c - simple test for models of hartmannograms
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
#include <stdio.h>
#include <math.h>
#include <getopt.h>
#include <stdarg.h>
#include <string.h>
#include "usefull_macros.h"
#include "cmdlnopts.h"
#include "zernike.h"
#include "spots.h"
#include "saveimg.h"
#include "benchmark.h"


void signals(int signo){
	exit(signo);
}

int main(int argc, char **argv){
	int i, j;
	//double scale;
	hartmann **images;
	mirror *mir = NULL;
	wavefront *comp;
	double *Zidxs = NULL;
	size_t L;
	polar *coords = NULL;
	point *grads = NULL;
	initial_setup();
	parse_args(argc, argv);
	if(G.zerngen){ // user give his zernike coefficients
		if(G.rest_pars_num < 1) ERRX(_("You should give at least one Zernike coefficient"));
		Zidxs = MALLOC(double, G.rest_pars_num);
		for(i = 0; i < G.rest_pars_num; ++i){
			if(!str2double(&Zidxs[i], G.rest_pars[i])){
				ERRX(_("Bad double number: %s!"), G.rest_pars[i]);
			}
		}
	}
	if(G.benchmark){
		DBG("Benchmark");
		do_benchmark(Zidxs, G.rest_pars_num);
		FREE(Zidxs);
		return 0;
	}
	if(G.zerngen){ // generate Zernike wavefront
		DBG("Generate Zernike surface");
		double *z = calc_surface(G.imagesize, Zidxs, G.rest_pars_num);
		if(z){
			writeimg(G.outfile, G.imagesize, mir, z);
			FREE(z);
		}
		return 0;
	}
	if(!G.gradname){
		DBG("Calculate wavefront by hartmannograms");
		if(G.rest_pars_num < 2) ERR(_("You should give at least two file names: pre- and postfocal!"));
		images = MALLOC(hartmann*, G.rest_pars_num + 1);
		for(i = 0; i < G.rest_pars_num; ++i)
			images[i] = read_spots(G.rest_pars[i]);
		mir = calc_mir_coordinates(images);
		getQ(mir);
		calc_Hartmann_constant(mir);
		spot_diagram *spot_dia = calc_spot_diagram(mir, mir->z07);
		//printf("\nmirror's coordinates should be corrected to (%g, %g)\n",
			//spot_dia->center.x, spot_dia->center.y);
		printf("Projection of center on mirror shifted by (%g, %g), image shifted by (%g, %g)\n",
			images[1]->center.x*HARTMANN_Z/distance, images[1]->center.y*HARTMANN_Z/distance,
			images[1]->center.x*(FOCAL_R-HARTMANN_Z)/distance, images[1]->center.y*(FOCAL_R-HARTMANN_Z)/distance);
		double tr = sqrt(images[1]->center.x*images[1]->center.x+images[1]->center.y*images[1]->center.y);
		printf("Beam tilt is %g''\n", tr/distance*206265.);
		calc_gradients(mir, spot_dia);
		coords = MALLOC(polar, mir->spotsnum);
		grads  = MALLOC(point, mir->spotsnum);
		printf("\nGradients of aberrations (*1e-6):\nray#      x         y         dx         dy         R         phi\n");
		for(j = 0, i = 0; i < 258; ++i){
			if(!mir->got[i]) continue;
			printf("%4d  %8.2f  %8.2f %10.6f %10.6f %10.2f %10.6f\n",
				i, mir->spots[i].x, mir->spots[i].y,
				mir->grads[i].x * 1e6, mir->grads[i].y * 1e6,
				mir->pol_spots[i].r * MIR_R, mir->pol_spots[i].theta * 180. / M_PI);
			memcpy(&coords[j], &mir->pol_spots[i], sizeof(polar));
			memcpy(&grads[j], &mir->grads[i], sizeof(point));
			grads[j].x *= 1e3;
			grads[j].y *= 1e3;
			j++;
		}
		L = mir->spotsnum;
		//scale = mir->Rmax;
		FREE(spot_dia);
		FREE(mir);
		for(i = 0; i < G.rest_pars_num; ++i)
			h_free(&images[i]);
		FREE(images);
	}else{
//		L = read_gradients(gradname, &coords, &grads, &scale);
		return 0;
	}
/*
	// spots information
	printf(GREEN "\nSpots:\n" OLDCOLOR "\n      r\ttheta\tDx(mm/m)\tDy(mm/m)\n");
	for(i = 0; i < L; i++){
		printf("%8.1f%8.4f%8.4f%8.4f\n", coords[i].r, coords[i].theta,
				grads[i].x*1000., grads[i].y*1000.);
	}
	*/
	// gradients decomposition (Less squares, Zhao polinomials)
	int Zsz, lastidx;
	Zidxs = gradZdecomposeR(10, L, coords, grads, &Zsz, &lastidx);
	//double *Zidxs = LS_gradZdecomposeR(10, L, coords, grads, &Zsz, &lastidx);
//	Zidxs[1] = 0.;
//	Zidxs[2] = 0.;
	lastidx++;
	printf("\nGradZ decompose, coefficients (%d):\n", lastidx);
	for(i = 0; i < lastidx; i++) printf("%g, ", Zidxs[i]);
	printf("\n\n");

	int GridSize = 15;
	comp = calc_surfaceR(GridSize, Zidxs, lastidx);
	if(comp){
		double zero_val = 0., N = 0., *zd = comp->zdata;
		polar *rect = comp->coordinates;
		for(j = GridSize-1; j > -1; j--){
			int idx = j*GridSize;
			for(i = 0; i < GridSize; i++,idx++){
				if(zd[idx] > 1e-9){
					zero_val += zd[idx];
					N++;
				}
			}
		}
		zero_val /= N;
		printf("zero: %g\nMirror surface:\n", zero_val);
		for(j = GridSize-1; j > -1; j--){
			int idx = j*GridSize;
			for(i = 0; i < GridSize; i++,idx++){
				if(rect[idx].r > 1. || rect[idx].r < 0.2){
					printf("           ");
				}else{
					printf("%7.3f", zd[idx]  * MIR_R);
				}
			}
			printf("\n");
		}
		free_wavefront(&comp);
	}
	// save matrix 100x100 with mirror deformations in Octave file format
	FILE *f = fopen("file.out", "w");
	if(f){
		GridSize = 100;
		comp = calc_surfaceR(GridSize, Zidxs, lastidx);
		if(comp){
			polar *rect = comp->coordinates;
			double *zd = comp->zdata;
			fprintf(f, "# generated by Z-BTA_test\n"
				"# name: dev_matrix\n# type: matrix\n# rows: %d\n# columns: %d\n",
					GridSize, GridSize);
			for(j = 0; j < GridSize; j++){
				int idx = j*GridSize;
				for(i = 0; i < GridSize; i++,idx++){
					if(rect[idx].r > 1. || rect[idx].r < 0.2)
						fprintf(f, "0. ");
					else
						fprintf(f, "%g ", zd[idx]  * MIR_R);
				}
				fprintf(f,"\n");
			}
			free_wavefront(&comp);
		}
		fclose(f);
	}
	// save image
	double *z = calc_surface(G.imagesize, Zidxs, lastidx);
	if(z){
		writeimg(G.outfile, G.imagesize, mir, z);
		FREE(z);
	}
	FREE(Zidxs);
	FREE(coords);
	FREE(grads);
	return 0;
}
