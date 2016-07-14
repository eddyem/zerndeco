/*
 * saveimg.c - functions to save data in png and FITS formats
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

#include "usefull_macros.h"
#include "cmdlnopts.h" // for flag "-r", which will tell to rewrite existing file
#include "saveimg.h"
#include <png.h>
#include <fitsio.h>
#include <math.h>

#define Stringify(x) #x
#define OMP_FOR(x) _Pragma(Stringify(omp parallel for x))

#ifndef THREAD_NUMBER
	#define THREAD_NUMBER 4
#endif

#ifndef OMP_NUM_THREADS
	#define OMP_NUM_THREADS THREAD_NUMBER
#endif


#define TRYFITS(f, ...)						\
do{int status = 0; f(__VA_ARGS__, &status);	\
	if (status){ ret = 0;					\
		fits_report_error(stderr, status);	\
		goto returning;}					\
}while(0)
#define WRITEKEY(...)						\
do{ int status = 0;							\
	fits_write_key(fp, __VA_ARGS__, &status);\
	if(status)								\
		fits_report_error(stderr, status);	\
}while(0)

/**
 * Create filename as   outfile + number + "." + suffix
 * 		number -- any number from 1 to 9999
 * This function simply returns "outfile.suffix" when "-r" option is set
 *
 * @param outfile - file name
 * @param suffix  - file suffix
 * @return created filename or NULL
 */
char *createfilename(char* outfile, char* suffix){
	FNAME();
	struct stat filestat;
	char buff[256], sfx[32];
	if(suffix) snprintf(sfx, 31, ".%s", suffix);
	else sfx[0] = 0; // no suffix
	if(rewrite_ifexists){ // there was key "-f": simply return copy of filename
		if(snprintf(buff, 255, "%s%s", outfile, sfx) < 1){
			DBG("error");
			return NULL;
		}
		DBG("(-r present) filename: %s", buff);
		return strdup(buff);
	}
	int num;
	if(!outfile) outfile = "";
	for(num = 1; num < 10000; num++){
		if(snprintf(buff, 255, "%s_%04d%s", outfile, num, sfx) < 1){
			DBG("error");
			return NULL;
		}
		if(stat(buff, &filestat)){ // || !S_ISREG(filestat.st_mode)) // OK, file not exists
			DBG("filename: %s", buff);
			return strdup(buff);
		}
	}
	DBG("n: %s\n", buff);
	WARN("Oops! All  numbers are busy or other error!");
	return NULL;
}

typedef struct{
	double *image;
	double min;
	double max;
	double avr;
	double std;
} ImStat;

static ImStat glob_stat;

/**
 * compute basics image statictics
 * @param img - image data
 * @param size - image size W*H
 * @return
 */
void get_stat(double *img, size_t size){
	FNAME();
	if(glob_stat.image == img) return;
	size_t i;
	double pv, sum=0., sum2=0., sz=(double)size;
	double max = -1., min = 1e15;
	for(i = 0; i < size; i++){
		pv = (double) *img++;
		sum += pv;
		sum2 += (pv * pv);
		if(max < pv) max = pv;
		if(min > pv) min = pv;
	}
	glob_stat.image = img;
	glob_stat.avr = sum/sz;
	glob_stat.std = sqrt(fabs(sum2/sz - glob_stat.avr*glob_stat.avr));
	glob_stat.max = max;
	glob_stat.min = min;
	DBG("Image stat: max=%g, min=%g, avr=%g, std=%g", max, min, glob_stat.avr, glob_stat.std);
}

/**
 * Save data to fits file
 * @param filename - filename to save to
 * @param sz  - image size: sz x sz
 * @param imbox - image bounding box
 * @data  image data
 * @return 0 if failed
 */
int writefits(char *filename, size_t sz, mirror *mir, double *data){
	FNAME();
	long naxes[2] = {sz, sz};
	static char* newname = NULL;
	char buf[80];
	int ret = 1;
	time_t savetime = time(NULL);
	fitsfile *fp;
	if(!filename) return 0;
	newname = realloc(newname, strlen(filename + 2));
	sprintf(newname, "!%s", filename); // say cfitsio that file could be rewritten
	TRYFITS(fits_create_file, &fp, newname);
	TRYFITS(fits_create_img, fp, FLOAT_IMG, 2, naxes);
	// FILE / Input file original name
	WRITEKEY(TSTRING, "FILE", filename, "Input file original name");
	WRITEKEY(TSTRING, "DETECTOR", "Wavefront model", "Detector model");
	// IMAGETYP / object, flat, dark, bias, scan, eta, neon, push
	WRITEKEY(TSTRING, "IMAGETYP", "object", "Image type");
	// DATAMAX, DATAMIN / Max,min pixel value
	WRITEKEY(TDOUBLE, "DATAMAX", &glob_stat.max, "Max data value");
	WRITEKEY(TDOUBLE, "DATAMIN", &glob_stat.min, "Min data value");
	// Some Statistics
	WRITEKEY(TDOUBLE, "DATAAVR", &glob_stat.avr, "Average data value");
	WRITEKEY(TDOUBLE, "DATASTD", &glob_stat.std, "Standart deviation of data value");
	// DATE / Creation date (YYYY-MM-DDThh:mm:ss, UTC)
	strftime(buf, 79, "%Y-%m-%dT%H:%M:%S", gmtime(&savetime));
	WRITEKEY(TSTRING, "DATE", buf, "Creation date (YYYY-MM-DDThh:mm:ss, UTC)");
	// DATE-OBS / DATE OF OBS.
	WRITEKEY(TSTRING, "DATE-OBS", buf, "DATE OF OBS. (YYYY-MM-DDThh:mm:ss, local)");
	// OBJECT  / Object name
	WRITEKEY(TSTRING, "OBJECT", "Modeled object", "Object name");
	// BINNING / Binning
	WRITEKEY(TSTRING, "XBIN", "1", "Horizontal binning");
	WRITEKEY(TSTRING, "YBIN", "1", "Vertical binning");
	// PROG-ID / Observation program identifier
	WRITEKEY(TSTRING, "PROG-ID", "Wavefront modeling", "Observation program identifier");
	// AUTHOR / Author of the program
	WRITEKEY(TSTRING, "AUTHOR", "Edward V. Emelianov", "Author of the program");
	if(mir){
		WRITEKEY(TDOUBLE, "ZBESTFOC", &mir->zbestfoc, "Z coordinate of best focus");
		WRITEKEY(TDOUBLE, "Z70",      &mir->z07,      "Z coordinate of best PSF on 70\%");
	}
	TRYFITS(fits_write_img, fp, TDOUBLE, 1, sz * sz, data);
	TRYFITS(fits_close_file, fp);
returning:
	return ret;
}

static uint8_t *rowptr = NULL;
uint8_t *processRow(double *irow, size_t width, double min, double wd){
	FREE(rowptr);
	//float umax = ((float)(UINT16_MAX-1));
	rowptr = MALLOC(uint8_t, width * 3);
	OMP_FOR(shared(irow))
	for(size_t i = 0; i < width; i++){
		double gray = (irow[i] - min) / wd;
		if(gray == 0.) continue;
		int G = (int)(gray * 4.);
		double x = 4.*gray - (double)G;
		uint8_t *ptr = &rowptr[i*3];
		uint8_t r = 0, g = 0, b = 0;
		switch(G){
			case 0:
				g = (uint8_t)(255. * x + 0.5);
				b = 255;
			break;
			case 1:
				g = 255;
				b = (uint8_t)(255. * (1. - x) + 0.5);
			break;
			case 2:
				r = (uint8_t)(255. * x + 0.5);
				g = 255;
			break;
			case 3:
				r = 255;
				g = (uint8_t)(255. * (1. - x) + 0.5);
			break;
			default:
				r = 255;
		}
		ptr[0] = r; ptr[1] = g; ptr[2] = b;
		//ptr[0] = ptr[1] = ptr[2] = gray*255;
	}
	return rowptr;
}

int writepng(char *filename, size_t sz, double *data){
	FNAME();
	int ret = 1;
	FILE *fp = NULL;
	png_structp pngptr = NULL;
	png_infop infoptr = NULL;
	double min = glob_stat.min, wd = glob_stat.max - min;
	double *row;

	if ((fp = fopen(filename, "w")) == NULL){
		perror("Can't open png file");
		ret = 0;
		goto done;
	}
	if ((pngptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
							NULL, NULL, NULL)) == NULL){
		perror("Can't create png structure");
		ret = 0;
		goto done;
	}
	if ((infoptr = png_create_info_struct(pngptr)) == NULL){
		perror("Can't create png info structure");
		ret = 0;
		goto done;
	}
	png_init_io(pngptr, fp);
	png_set_compression_level(pngptr, 1);
	png_set_IHDR(pngptr, infoptr, sz, sz, 8, PNG_COLOR_TYPE_RGB,//16, PNG_COLOR_TYPE_GRAY,
				PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
				PNG_FILTER_TYPE_DEFAULT);
	png_write_info(pngptr, infoptr);
	png_set_swap(pngptr);
	size_t height = sz;
	for(row = &data[sz*(height-1)]; height > 0; row -= sz, height--)
		png_write_row(pngptr, (png_bytep)processRow(row, sz, min, wd));
	png_write_end(pngptr, infoptr);
done:
	if(fp) fclose(fp);
	if(pngptr) png_destroy_write_struct(&pngptr, &infoptr);
	return ret;
}

/**
 * Save data to image file[s] with format t
 * @param name - filename prefix or NULL to save to "outXXXX.format"
 * @param t - image[s] type[s]
 * @param width, height - image size
 * @param imbox - image bounding box (for FITS header)
 * @param mirror - mirror parameters (for FITS header)
 * @param data  image data
 * @return number of saved images
 */
int writeimg(char *name, size_t sz, mirror *mir, double *data){
	FNAME();
	char *filename = NULL;
	int ret = 0;
	get_stat(data, sz*sz);
	filename = createfilename(name, "fits");
	if(filename){
		ret = writefits(filename, sz, mir, data);
		FREE(filename);
	}
	filename = createfilename(name, "png");
	if(filename){
		ret += writepng(filename, sz, data);
		FREE(filename);
	}
	return ret;
}
