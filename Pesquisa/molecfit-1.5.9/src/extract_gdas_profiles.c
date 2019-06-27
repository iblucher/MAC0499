/*
 * extract_gdas_profiles.c
 *
 * Programme to extract site-related records of ARL packed meteorological
 * files provided by a list.
 * The time-specific output files (3h resolution; date/time format: YYMMDDHH)
 * contain pressure, geopotential height, temperature, dew point
 * (approximation), wind direction, and wind speed.
 *
 * Author:  Stefan Noll, Wolfgang Kausch, University of Innsbruck, Austria
 * Created: 28 Oct 2013
 * History: based on a FORTRAN procedure for GDAS file checking provided by
 *          the NOAA
 *
 * Compilation:
 * gcc gdas_extract.c -o gdas_extract -lm
 *
 * Input parameters:
 * (1) list of file names (creation by, e.g., "ls gdas1.* >gdas.list")
 * (2) site (P = Paranal; for other sites programme enquires coordinates)
 * (3) output folder (default: "outdat"; must exist!)
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Number of characters for parameter values in GDAS files */
#define LENLINE 80
#define LENLAB 50
#define LENTIME 8
#define LENLEV 2
#define LENID 4
#define LENEXP 4
#define LENVAR 14
#define NX 360
#define NY 181
#define LENPACK NX * NY
#define NLEV 24

/* Position of parameters in GDAS record label */
#define POSTIME 0
#define POSLEV 10
#define POSID 14
#define POSEXP 18
#define POSVAR 36

/* Coordinates of Paranal */
#define LON_P -70.4
#define LAT_P -24.6
#define LON_A -70.2  /* Cerro Armazones */
#define LAT_A -24.6
#define LON_L -70.7  /* La Silla */
#define LAT_L -29.3
#define LON_C -67.8  /* Cajnantor */
#define LAT_C -23.0

/* Default input parameters */
char *site0;
char *outdir0;

/* Level-specific pressure values */
/*double press[NLEV] = {1013., 1000., 975., 950., 925., 900., 850., 800.,
                      750., 700., 650., 600., 550., 500., 450., 400., 350.,
                      300., 250., 200., 150., 100., 50., 20.};*/
double press[NLEV] = {1013., 1000., 975., 950., 925., 900., 850., 800.,
                      750., 700., 650., 600., 550., 500., 450., 400., 350.,
                      300., 250., 200., 150., 100., 50., 20.};

/* Declaration of functions */
void split_label(char [LENLAB+1], char [LENTIME+1], int *, char [LENID+1],
                 int *, double *);
double unpack_site(unsigned char [LENPACK], int, double, double, double);
double DewPoint(double, double);


int main(int argc, char *argv[])
{

    FILE *list=NULL, *indata=NULL, *outdata=NULL;
    char *filelist=NULL, *site=NULL, *outdir=NULL;
    char line[LENLINE], infile[LENLINE], outfile[LENLINE];

    int count = 0;
    int check, krec;
    size_t namelen;

    char label[LENLAB+1];
    unsigned char cpack[LENPACK];

    char time[LENTIME+1], kvar[LENID+1], yr[5],mn[3],dy[3],hr[3];
    int kl, nexp;
    double var1, lon, lat, val;

    double relhum[NLEV] = {0.};
    double temp[NLEV], hgeopot[NLEV];

    /* Get parameters */
    if (argc != 4 && argc != 6) {
        printf("Wrong number of arguments!\n");
        printf("Usage: extract_gdas_profiles <listoffiles> <site> "
               "<output dir>\n");
        printf("   <listoffiles> = ASCII file of downloaded GDAS archives\n");
        printf("   <site>        = P (Paranal), L (La Silla), C (Cajnantor), "
                                  "A (Armazones)\n");
        printf("   <output dir>  = abs path to output dir\n");
        exit(EXIT_FAILURE);
    } else {
        filelist = (char *) argv[1];
        site = (char *) argv[2];
        if (argc == 4) {
            outdir = (char *) argv[3];
        }
        else {
            outdir = (char *) argv[5];
        }
    }

//     printf("Command-line parameters: %s %s %s\n",filelist, site, outdir);

    if (strlen(outdir) + strlen(site) + LENTIME + 5 > LENLINE - 1) {
        printf("Path + name of output files too long! -> Exit!\n");
        exit(EXIT_FAILURE);
    }

    /* Longitude and latitude of site */
    if (strcmp(site, "custom") == 0) {
        lon = NAN;
        lat = NAN;
        sscanf(argv[3], "%lf", &lon);
        sscanf(argv[4], "%lf", &lat);
        if (lon < -180 || lon > 180 || lat < -90 || lat > 90) {
            printf("Input coordinates are not valid! -> Exit!\n");
            exit(EXIT_FAILURE);
        }
    }
    else if (site[0] == 'P' || site[0] == 'p') {
        /* Coordinates of Paranal */
        lon = LON_P;
        lat = LAT_P;
    } else if (site[0] == 'C' || site[0] == 'c') {
        /* Coordinates of Paranal */
        lon = LON_C;
        lat = LAT_C;
    } else if (site[0] == 'A' || site[0] == 'a') {
        /* Coordinates of Paranal */
        lon = LON_A;
        lat = LAT_A;
    } else if (site[0] == 'L' || site[0] == 'l') {
        /* Coordinates of Paranal */
        lon = LON_L;
        lat = LAT_L;
    } else {
        printf("Site does not exist!\n");
        printf("Provide longitude (W negative) and latitude (S negative):\n");
        check = scanf("%lf %lf", &lon, &lat);
        getchar();
        if (check != 2 || lon < -180 || lon > 180 || lat < -90 || lat > 90) {
            printf("Input coordinates are not valid! -> Exit!\n");
            exit(EXIT_FAILURE);
        }
    }

    printf("Coordinates: %+.1f %+.1f\n", lon, lat);

    /* Read file names */

    if ((list = fopen(filelist, "r")) == NULL) {
        printf("Cannot open %s!\n", filelist);
        exit(EXIT_FAILURE);
    }

    while (fgets(line, LENLINE, list) != NULL) {

        count++;
        namelen = strlen(line);
        strncpy(infile, line, namelen);
        infile[namelen-1] = '\0';
        printf("%s\n", infile);

        /* Open GDAS file */

        if ((indata = fopen(infile, "r")) == NULL) {
            printf("Cannot open %s! -> Skip!\n", infile);
            continue;
        }

        krec = 0;

        while (!feof(indata)) {

            /* Read GDAS record */
            if (fgets(label, LENLAB+1, indata)) {};
            if (fread(cpack, 1, LENPACK, indata)) {};

            /* Extract global parameters of record */
            split_label(label, time, &kl, kvar, &nexp, &var1);
            /* Extract data value for given site */
            val = unpack_site(cpack, nexp, var1, lon, lat);

            /* Collect relevant data and calculate output properties */

            if (!strncmp(kvar,"INDX",LENID)) {
                if (krec != 0) {
                    for (kl = 1; kl < NLEV; kl++) {
                            fprintf(outdata, " %5.0f\t %6.3f\t%5.1f\t%5.1f\n",
                                    press[kl], hgeopot[kl], temp[kl],
                                    relhum[kl]);
                    }
                    fclose(outdata);
                }
            } else if (!strncmp(kvar,"PRSS",LENID)) {
                press[kl] = val;
                /* Open output file and write header */
                sprintf(yr,"%c%c",time[0],time[1]);
                sprintf(mn,"%c%c",time[2],time[3]);
                sprintf(dy,"%c%c",time[4],time[5]);
                sprintf(hr,"%c%c",time[6],time[7]);

                sprintf(outfile,"%s/C%+2.1f%+2.1fD20%s-%s-%sT%s.gdas",
                        outdir, lon, lat, yr, mn, dy, hr);
                if ((outdata = fopen(outfile,"w")) == NULL) {
                    printf("Output directory %s does not exist! -> Exit!\n",
                           outdir);
                    exit(EXIT_FAILURE);
                }
                fprintf(outdata,"# P[hPa] HGT[m]         T[K]    "
                                "RELHUM[%%]\n");
            } else if (!strncmp(kvar,"RH2M",LENID) ||
                       !strncmp(kvar,"RELH",LENID)) {
                relhum[kl] = val;
            } else if (!strncmp(kvar,"T02M",LENID) ||
                       !strncmp(kvar,"TEMP",LENID)) {
                temp[kl] = val;
            } else if (!strncmp(kvar,"SHGT",LENID) ||
                       !strncmp(kvar,"HGTS",LENID)) {
              hgeopot[kl] = val;
            }

            krec++;

        }

        fclose(indata);

        for (kl = 1; kl < NLEV; kl++) {
                fprintf(outdata, " %5.0f\t %6.3f\t%5.1f\t%5.1f\n",
                        press[kl], hgeopot[kl], temp[kl], relhum[kl]);
        }
        fclose(outdata);

    }

    printf("%s: %d files decoded\n",filelist,count);

    fclose(list);

    return EXIT_SUCCESS;

}


void split_label(char label[LENLAB+1], char time[LENTIME+1], int *kl,
                 char kvar[LENID+1], int *nexp, double *var1)
{
    /*
     * Split label of GDAS record.
     * Extract date/time, pressure level, variable ID, exponent for decoding,
     * and latitude-specific start value for decoding.
     */

    char level[LENLEV+1], expc[LENEXP+1], var1c[LENVAR+1];
    int i;

    for (i = 0; i < LENLAB; i++) {
        if (i >= POSTIME && i < POSTIME+LENTIME) {
            if (label[i] == ' ') {
                time[i-POSTIME] = '0';
            } else {
                time[i-POSTIME] = label[i];
            }
        } else if (i >= POSLEV && i < POSLEV+LENLEV) {
            level[i-POSLEV] = label[i];
        } else if (i >= POSID && i < POSID+LENID) {
            kvar[i-POSID] = label[i];
        } else if (i >= POSEXP && i < POSEXP+LENEXP) {
            expc[i-POSEXP] = label[i];
        } else if (i >= POSVAR && i < POSVAR+LENVAR) {
            var1c[i-POSVAR] = label[i];
        }
    }

    time[LENTIME] = '\0';
    level[LENLEV] = '\0';
    *kl = atoi(level);
    kvar[LENID] = '\0';
    expc[LENEXP] = '\0';
    *nexp = atoi(expc);
    var1c[LENVAR] = '\0';
    *var1 = atof(var1c);
}


double unpack_site(unsigned char cpack[LENPACK], int nexp, double var1,
                   double lon, double lat)
{
    /* Unpack GDAS record for location given by longitude and latitude */

    int i0, j0, indx, j, i;
    double di, dj, scale, vold;
    double rdata[NX][NY];

    /* Derivation of grid coordinates */
    if (lon < 0) lon = 360. + lon;
    i0 = (int) lon;
    di = lon - i0;
    lat = 90. + lat;
    j0 = (int) lat;
    dj = lat - j0;

    scale = pow(2.0, (7 - nexp));
    vold = var1;
    indx = 0;

    for (j = 0; j < NY; j++) {
        for (i = 0; i < NX; i++) {
            rdata[i][j] = ((double) cpack[indx] - 127.) / scale + vold;
            vold = rdata[i][j];
            indx++;
        }
        vold = rdata[0][j];
    }

    /* Bilinear interpolation (approximation) */
    return rdata[i0][j0] * (1 - di) * (1 - dj)
      + rdata[i0+1][j0] * di * (1 - dj)
      + rdata[i0][j0+1] * (1 - di) * dj
      + rdata[i0+1][j0+1] * di * dj;
}
