/*
 * main.cpp
 *
 *  Created on: Jan 19, 2016
 *      Author: vova
 */
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include "cspline.h"


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define MAX_HEADER_LENGTH       2048
#define MAX_INPUT_LINE_LENGTH   2048
#define PI                      M_PI

typedef
enum compression_t {
    compression_unknown = -1,
    compression_none,
    compression_bzip,
    compression_gzip
} compression_t;
//using namespace std;

/** uses fopen() if input file seems to be uncompresssed, and popen() if input file seems compressed */
static FILE * open_file( const char * fname, compression_t * compression )
{
  char cmd[PATH_MAX] = {0};
  FILE * input = NULL;

  if ( *compression == compression_unknown )
  {
    const char * suffix;

    if ( (suffix = strstr(fname, ".bz2")) && *(suffix + 4) == 0 ) {
      *compression = compression_bzip;
    }
    else if ( (suffix = strstr(fname, ".bz")) && *(suffix + 3) == 0 ) {
      *compression = compression_bzip;
    }
    else if ( (suffix = strstr(fname, ".gz")) && *(suffix + 3) == 0 ) {
      *compression = compression_gzip;
    }
    else {
      *compression = compression_none;
    }
  }

  switch ( *compression )
  {
  case compression_bzip:
    snprintf(cmd, sizeof(cmd) - 1, "bzip2 -dc '%s'", fname);
    if ( !(input = popen(cmd, "r"))) {
      fprintf(stderr, "popen('%s') fails: %s\n", cmd, strerror(errno));
    }
    break;

  case compression_gzip:
    snprintf(cmd, sizeof(cmd) - 1, "gzip -dc '%s'", fname);
    if ( !(input = popen(cmd, "r"))) {
      fprintf(stderr, "popen('%s') fails: %s\n", cmd, strerror(errno));
    }
    break;

  case compression_none:
    if ( !(input = fopen(fname, "r")) ) {
      fprintf(stderr, "fopen('%s') fails: %s\n", fname, strerror(errno));
    }
    break;

  default:
    fprintf(stderr, "BUG IN CODE: invalid compression tag=%d\n", *compression);
    break;
  }

  return input;
}

static void close_file( FILE * input, compression_t compression )
{
  if ( input && input != stdin ) {
	if ( compression > compression_none ) {
	  pclose(input);
	}
	else {
	  fclose(input);
	}
  }
}


void skipline(FILE * fr)
{
	while(!feof(fr))
		if(fgetc(fr)=='\n')
			break;
}
//static long int correction_objects( const char * fname, int cmag, int cmura, int cmude,  int print_dx)
static long int correction_objects( const char * fname, int cmag, int cmura, int cmude,int print_dx,cubic_spline SMURA,cubic_spline SMUDE, double maxmagcor)
{
	FILE * corfile = fopen(fname,"r");
	if(!corfile )
	{
		fprintf(stderr,"cannot load data %s file\n",fname);
	}

	char line[MAX_INPUT_LINE_LENGTH];
	int ic;
	char * pc;
	double mag,mura,mude;

	/*
	char header[MAX_HEADER_LENGTH];

	header[MAX_HEADER_LENGTH-1]=0;


	  if ( !fgets(header, MAX_HEADER_LENGTH, corfile) ) {
	    fprintf(stderr,"fgets(header) fails: %d (%s)\n", errno, strerror(errno));
	    return -1;
	  }

	  if ( header[MAX_HEADER_LENGTH - 1] != 0 ) {
	    fprintf(stderr,"too long header line in this file\n");
	    return -1;
	  }
*/
	  if(print_dx==1)
	  {
		  printf("mag\tfit_dmura\tfit_dmude\tcormura\tcormude\n");
	  }
	  else{
		  //printf("%s\tcormura\tcormude\n",header);
	  }

	long int ii=0;

	while ( !feof(corfile) )
	{

	  //printf("%.10f\t%.10f\t%.2f\t%.2f\n",m[i].x,m[i].y,m[i].f,m[i].ef);

      line[MAX_INPUT_LINE_LENGTH-1] = 0;
      if ( !fgets(line, MAX_INPUT_LINE_LENGTH, corfile) ) {
        break;
      }


      if ( line[MAX_INPUT_LINE_LENGTH - 1] != 0 ) {
        fprintf(stderr,"too long input line in this file\n");
        return -1;
      }

      /* remove trailing new line */
      line[strlen(line) - 1] = 0;


      for ( pc = line, ic = 1; ic < cmag; ++ic ) {
        if ( !(pc = strchr(pc + 1, '\t')) ) {
          break;
        }
      }
      if ( ic != cmag || sscanf(pc, " %lf",&mag) != 1 ) {
        continue;
      }

      for ( pc = line, ic = 1; ic < cmura; ++ic ) {
        if ( !(pc = strchr(pc + 1, '\t')) ) {
          break;
        }
      }
      if ( ic != cmura || sscanf(pc, " %lf", &mura) != 1 ) {
        continue;
      }

      for ( pc = line, ic = 1; ic < cmude; ++ic ) {
		  if ( !(pc = strchr(pc + 1, '\t')) ) {
			break;
		  }
		}
		if ( ic != cmude || sscanf(pc, " %lf", &mude) != 1 ) {
		  continue;
		}

		double cor_mura=SMURA.fspline(mag);
		double cor_mude=SMUDE.fspline(mag);
		//double cor_mura=0;
		//double cor_mude=0;

		if(maxmagcor>0 && mag>maxmagcor)
		{
			cor_mura=0;
			cor_mude=0;
		}

      if(print_dx==1)
          {
    	  printf("%6.3f\t%10.3f\t%10.3f\t%10.3f\t%10.3f\n",mag,cor_mura,cor_mude,mura+cor_mura,mude+cor_mude);
          }
      else{printf("%s\t%0.2f\t%.2f\n",strdup(line),mura+cor_mura,mude+cor_mude);}

      ii++;

	}

	fclose(corfile);
return ii;
}


static long int load_objects(FILE * fr, int mag, int dmura, int dmude, double * magdata,double * dmuradata,double * dmudedata)
{

		char line[MAX_INPUT_LINE_LENGTH];
		int ic;
		char * pc;

		int count=0;

		while ( !feof(fr) )
		{

	      line[MAX_INPUT_LINE_LENGTH-1] = 0;
	      if ( !fgets(line, MAX_INPUT_LINE_LENGTH, fr) ) {
	        break;
	      }


	      if ( line[MAX_INPUT_LINE_LENGTH - 1] != 0 ) {
	        fprintf(stderr,"too long input line in this file\n");
	        return -1;
	      }

	      /* remove trailing new line */
	      line[strlen(line) - 1] = 0;


	      for ( pc = line, ic = 1; ic < mag; ++ic ) {
	        if ( !(pc = strchr(pc + 1, '\t')) ) {
	          break;
	        }
	      }
	      if ( ic != mag || sscanf(pc, " %lf",&magdata[count]) != 1 ) {
	        continue;
	      }

	      for ( pc = line, ic = 1; ic < dmura; ++ic ) {
			if ( !(pc = strchr(pc + 1, '\t')) ) {
			  break;
			}
		  }
		  if ( ic != dmura || sscanf(pc, " %lf",&dmuradata[count]) != 1 ) {
			continue;
		  }

	      for ( pc = line, ic = 1; ic < dmude; ++ic ) {
			if ( !(pc = strchr(pc + 1, '\t')) ) {
			  break;
			}
		  }
		  if ( ic != dmude || sscanf(pc, " %lf",&dmudedata[count]) != 1 ) {
			continue;
		  }
		  count++;

		}

		return count;

}


static void show_usage( FILE * output, int argc, char * argv[] )
{
  fprintf(output, "Cubic_spline regresion UTILITY\n");
  fprintf(output, "USAGE:\n");
  fprintf(output, "  %s OPTIONS DATA FILE\n", basename(argv[0]));
  fprintf(output, "  %s OPTIONS CORRECTED_FILE\n", basename(argv[1]));
  fprintf(output, "OPTIONS:\n");
  fprintf(output, "	mag=	<integer>	colum number for mag column in data file\n");
  fprintf(output, "	dmura=	<integer>	colum number for dmura column in data file\n");
  fprintf(output, "	dmude=	<integer>	colum number for dmude column in data file\n");
  fprintf(output, "	cmag=	<integer>	colum number for mag column in corrected file\n");
  fprintf(output, "	cmura=	<integer>	colum number for mura column in corrected file\n");
  fprintf(output, "	cmude=	<integer>	colum number for mude column in corrected file\n");
  fprintf(output, "	npoint=	<integer>	count used point for spline, default npoint=1000\n");
  fprintf(output, "	-maxmag= <double>	max magnitude for corection, default maxmag = 0 unused\n");
  fprintf(output, "	-v					Print some diagnostic messages to stderr (verbose mode)\n");
  fprintf(output, "	-d					Print print only dmu\n");
  fprintf(output, "	\n");
}

int main(int argc, char *argv[])
{
		const char * fname_in = {NULL};
		const char * fname_out = {NULL};

		FILE * fp; //FILE

		compression_t compression ={ compression_unknown};

		int mag=-1;		//	colum number for mag column in data file
		int dmura=-1;	//	colum number for dmura column in data file
		int dmude=-1;	//	colum number for dmude column in data file
		int cmag=-1;	//	colum number for mag column in corrected file
		int cmura=-1;	//	colum number for mura column in corrected file
		int cmude=-1;	//	colum number for mude column in corrected file

		double maxmag=0;

		int npoint=100;

		int print_dmu=-1;

		int beverbose 	= 0; 		//	Print some diagnostic messages to stderr (verbose mode)


		int usedfiles=0;
		//double max_x = 0,min_x=1000000,max_y=0,min_y=1000000;

		/* parse command line */

		int i;
		for ( i = 1; i < argc; ++i )
		{
			 if ( strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-help") == 0 ) {
				  show_usage(stdout, argc, argv);
				  return 0;
				}

			/* read used columns from file1*/
				if ( strncmp(argv[i], "mag=", 4) == 0 )
				{
				  if ( sscanf(argv[i] + 4, "%d", &mag) != 1 )
				  {
					fprintf(stderr, "Invalid value of %s\n", argv[i]);
					return 1;
				  }
				}
				else if ( strncmp(argv[i], "dmura=", 6) == 0 )
				{
				  if ( sscanf(argv[i] + 6, "%d", &dmura) != 1 )
				  {
					fprintf(stderr, "Invalid value of %s\n", argv[i]);
					return 1;
				  }
				}
				else if ( strncmp(argv[i], "dmude=", 6) == 0 )
				{
				  if ( sscanf(argv[i] + 6, "%d", &dmude) != 1 )
				  {
					fprintf(stderr, "Invalid value of %s\n", argv[i]);
					return 1;
				  }
				}
				else if ( strncmp(argv[i], "cmag=", 5) == 0 )
				{
				  if ( sscanf(argv[i] + 5, "%d", &cmag) != 1 )
				  {
					fprintf(stderr, "Invalid value of %s\n", argv[i]);
					return 1;
				  }
				}
				else if ( strncmp(argv[i], "cmura=", 6) == 0 )
				{
				  if ( sscanf(argv[i] + 6, "%d", &cmura) != 1 )
				  {
					fprintf(stderr, "Invalid value of %s\n", argv[i]);
					return 1;
				  }
				}
				else if ( strncmp(argv[i], "cmude=", 6) == 0 )
				{
				  if ( sscanf(argv[i] + 6, "%d", &cmude) != 1 )
				  {
					fprintf(stderr, "Invalid value of %s\n", argv[i]);
					return 1;
				  }
				}
			else if ( strncmp(argv[i], "npoint=", 7) == 0 )
				{
				  if ( sscanf(argv[i] + 7, "%d", &npoint) != 1 || npoint < 1 )
				  {
					fprintf(stderr, "Invalid value of %s\n", argv[i]);
					return 1;
				  }
				}
			else if ( strncmp(argv[i], "-maxmag=", 8) == 0 )
				{
				  if ( sscanf(argv[i] + 8, "%lf", &maxmag) != 1 || maxmag < 1 )
				  {
					fprintf(stderr, "Invalid value of %s\n", argv[i]);
					return 1;
				  }
				}

			else if ( strcmp(argv[i], "-v") == 0 ) {
				  beverbose = 1;
				}
			else if ( strcmp(argv[i], "-d") == 0 ) {
				  print_dmu = 1;
				}
			else if ( !fname_in ) {
				  fname_in = argv[i];
				  usedfiles = 1;
				}
			else if ( !fname_out ) {
				  fname_out = argv[i];
				  usedfiles = 2;
				}
			else
				{
				  fprintf(stderr, "Invalid argument %s. Try %s --help\n", argv[i], argv[0]);
				  //printf("load data from file %s",fname);
				  return 1;
				}


		  }

		/* check command line inputs */
		  if ( !fname_in) {
			fprintf(stderr,"Input file1 names expected\n");
			show_usage(stderr, argc, argv);
			return -1;
		  }
		  if ( !fname_out && usedfiles==2) {
		  		fprintf(stderr,"Input file3 names expected\n");
		  		show_usage(stderr, argc, argv);
		  		return -1;
		  	  }

		  if ( mag < 1 && usedfiles==1) {
			  fprintf(stderr,"mag argument is mandatory\n");
			  show_usage(stderr, argc, argv);
			  return -1;
			}
		  if ( dmura < 1 && usedfiles==1) {
			  fprintf(stderr,"dmura argument is mandatory\n");
			  show_usage(stderr, argc, argv);
			  return -1;
		  }
		  if ( dmude < 1 && usedfiles==1) {
			  fprintf(stderr,"dmude argument is mandatory\n");
			  show_usage(stderr, argc, argv);
			  return -1;
		  }
		  if ( cmag < 1 && usedfiles==2) {
			  fprintf(stderr,"cmag argument is mandatory\n");
			  show_usage(stderr, argc, argv);
			  return -1;
			}
		  if ( cmura < 1 && usedfiles==2) {
			  fprintf(stderr,"cmura argument is mandatory\n");
			  show_usage(stderr, argc, argv);
			  return -1;
		  }
		  if ( cmude < 1 && usedfiles==2) {
			  fprintf(stderr,"cmude argument is mandatory\n");
			  show_usage(stderr, argc, argv);
			  return -1;
		  }
		  if ( beverbose ) {
			  if(usedfiles==2)
			  {
			  fprintf(stderr,"Used file_data '%s' and out_file '%s '\n",fname_in,fname_out);
			  }else{
				  fprintf(stderr,"Used file_data '%s'\n",fname_in);
			  }
			  fprintf(stderr,"Used mag=%d dmura=%d dmude=%d npoint=%d in '%s' file and cmag=%d cmura=%d cmude=%d in '%s '\n",mag,dmura,dmude,npoint,fname_in,cmag,cmura,cmude,fname_out);
		  }

		  if (usedfiles>0)
		  	{
		  		/* check if input files are readable */
		  		    if ( access(fname_in, R_OK) != 0 ) {
		  		      fprintf(stderr, "Can't read %s: %s\n", fname_in, strerror(errno));
		  		      return -1;
		  		    }


		  		  if ( beverbose ) {
		  			fprintf(stderr,"loading %s....\n", fname_in);
		  		  }


		  		if ( !(fp = open_file(fname_in, &compression)) ) {
		  					fprintf(stderr, "Can't read '%s': %s\n", fname_in, strerror(errno));
		  					return -1;
		  				  }

		  		double * mag_data 	= 	new double [npoint];
		  		double * dmura_data	= 	new double [npoint];
		  		double * dmude_data	=	new double [npoint];


				  int nstars=load_objects(fp,mag,dmura,dmude,mag_data,dmura_data,dmude_data);

				  if ( nstars <1){
					fprintf(stderr, "Can't load data from %s file\n", fname_in);
					return -1;
				 }
				close_file( fp, compression );

		  		if (beverbose)
		  		{
		  			fprintf(stderr,"load %d rows from %s \n", nstars,fname_in);
		  		}


				/*for ( int n=0; n < nstars; ++n )
				  {
					fprintf(stderr,"%lf\t%lf\t%lf\n",mag_data[n],dmura_data[n],dmude_data[n]);
				  }
				  */


		  		cubic_spline SP_mura;
		  		cubic_spline SP_mude;


				SP_mura.build_spline(mag_data,dmura_data,nstars);
				SP_mude.build_spline(mag_data,dmude_data,nstars);


				int ncorection =correction_objects( fname_out, cmag, cmura, cmude,print_dmu,SP_mura, SP_mude, maxmag);


		  		if (beverbose)
		  		{
		  			fprintf(stderr,"corrected %d objects from %s \n", ncorection,fname_out);
		  		}


		  		delete [] mag_data;
		  		delete [] dmura_data;
		  		delete [] dmude_data;

		  		//SP_mura.~cubic_spline();
		  		//SP_mude.~cubic_spline();
		  	}

	return 0;
}
