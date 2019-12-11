#include <stdio.h>
#include <math.h>
#include <string.h>

/* ephem.c
 * Ephemeris and Almanac calculations.
 *
 * (c) Tom Trebisky  10-4-2017 10-13-2017
 */

/* TODO
 * -- split into multiple files
 * -- rename time to "day" structure
 * -- add mjd to date and time routine
 * -- add command line options to drive this thing
 * -- Anewm calculation for any year
 */

/* Note that the international convention for longitude is positive to the east,
 * so all of the longitude coordinates in Arizona ought to be negative.
 * We maintain the same convention for timezones (so it is -7 hours in Arizona).
 */

struct site_data {
	double lat_deg;
	double long_deg;
	double elevation;
	double tz;
	char *name;
};

/* elevation in meters */
struct site_data site_mmt = { 31.68877778, -110.88456, 2606.0, -7.0, "mmt" };
struct site_data site_castellon = { 32.2627640, -111.0485611, 736.092, -7.0, "castellon" };
struct site_data site_tucson = { 32.195, -110.892, 700.0, -7.0, "tucson" };
struct site_data site_greenwich = { 51.4826, 0.0077, 46.0, 0.0, "greenwich" };

/* These are for the MMT (Mt. Hopkins) */
#define LATITUDE        31.68877778     /* site latitude in degrees */
#define LONGITUDE       7.392303704     /* site longitude in hours west */
#define ELEVATION       2.606           /* site elevation in kilometers */
#define TIMEZONE	7.0		/* hours from greenwich */

/* These are for the earth */
#define F               (1.0/298.257)   /* oblateness constant */
#define RADIUS          6378.160        /* radius of the earth in kilometers */

#ifdef notdef
/* These are for the VATT (Mt. Graham) */
#define LATITUDE        32.70133333     /* site latitude in degrees */
#define LONGITUDE       7.326103704     /* site longitude in hours west */
#define ELEVATION       3.171           /* site elevation in kilometers */
#define TIMEZONE	7.0		/* hours from greenwich */

/* These are for Greenwich, England */
#define LATITUDE        51.476852
#define LONGITUDE       (0.000500 * 24.0 / 360.0)	/* 1.8 arc-seconds West */
#define ELEVATION       0.046				/* 46 meters */
#define TIMEZONE	0.0

/* These are for Castellon, Tucson, Arizona */
#define LATITUDE        32.2627640
#define LONGITUDE       111.0485611	/* West */
#define ELEVATION       0.736092	/* 2415 feet */
#define TIMEZONE	7.0		/* hours from greenwich */
/*
#define LATITUDE        32:15:45.909
#define LONGITUDE       111:02:54.82
*/
#endif

#define PI	3.14159265358979324
#define TWOPI	(2.0 * PI)
#define DEGRAD	(PI/180.0)
#define RADDEG	(180.0/PI)
#define RADHR	(12.0/PI)
#define HRRAD	(PI/12.0)

/* Seconds of time in a day.
 * i.e. 24*60*60
 */
#define SECPD	86400.0

/* Multiply by these to perform the indicated conversion */
#define ASEC_TO_RAD	(DEGRAD / 3600.0)
#define RAD_TO_ASEC	(3600.0 * RADDEG)	/* APC: Arcs */

#define JULIAN_CENTURY	36525.0		/* Days in a Julian century */
#define MJD_2000	51544.5		/* MJD for J2000 epoch */

/*
 * The 2003 Almanac page K5 says that
 *  1900 January 0 at 0h is JD = 2415019.5
 * Note that JD is zero at noon, whereas
 *  MJD is zero at midnight.
 * These are standard conventions.
 */
#define MJD_OFFSET	2400000.5
#define JD_1900		2415020.0

/* Arc seconds in a full circle.
 * i.e. 360*60*60
 */
#define ARCSEC_360	1296000.0

struct site {
	double lat_deg;
	double long_deg;
	double long_hours;
	double latitude;
	double longitude;
	double elevation;
	double tz;
	double sin_lat;
	double cos_lat;
	char name[64];
};

struct day {
	int *dpm;
	int days_in_month;
	int year;
	int month;
	int day;
	double mjd0;
	double jd0;
};

struct site site_info;

void init_site ( struct site *, struct site_data * );
void set_day ( struct day *, int, int, int );
double calc_mjd ( struct day * );
double apc_st ( double );

char * s_dms ( char *, double );
char * s_dms_b ( char *, double );
char * s_dm ( char *, double );
char * s_dm_b ( char *, double );

double hms ( int, int, double );
double mmt_st ( double );
double mmt_lst ( double );

void sun_apc_ll ( double, double *, double *, int );
void sun_apc ( double, double *, double *, int );

void sun_wiki ( double, double * );

void moon_apc_ll ( double, double *, double * );
void moon_apc ( double, double *, double * );
// void MiniMoon_ll ( double, double *, double * );
// void MiniMoon ( double, double *, double * );

void test_jd ( void );
void test_st ( void );
void test_sun1 ( void );
void test_sun2 ( void );
void test_sun3a ( void );
void test_sun3b ( void );
void test_almanac ( void );
void test_anew ( int );
void anew_table ( void );
void show_old_newm ( void );
char * get_today_str ( void );

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */


/* Days per month */
int dpm_normal[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
int dpm_leap[] = { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

char *month_name[] = { "Jan", "Feb", "Mar", "Apr", "May", "Jun",
			"Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };

char *month_name_full[] = { "January", "February", "March", "April", "May", "June",
			"July", "August", "September", "October", "November", "December" };

/*
 * We have leap years because a tropical year
 *  is actually 365.242190 days in length.
 */
int
is_leap ( int year )
{
	if ( (year % 4) != 0 ) 
	    return 0;
	if ( (year % 100) == 0 && (year % 400) != 0 )
	    return 0;
	return 1;
}

/* Calculate the MJD for a given year, month, and day
 *  For years after 1582, we are using the Gregorian calendar,
 *  so for my purposes we could eliminate the test below.
 * This expects month as 1 .. 12
 * This expects day as 1 ... N
 *
 * MJD is zero at midnight.
 * JD is zero at noon.
 */
double 
calc_mjd ( struct day *tp )
{
	int m;
	int y;
	int b;
	int a;
	double mjd;

	m = tp->month;
	y = tp->year;
	if ( m <= 2 ) {
	    m += 12;
	    y--;
	}

	if ( y * 10000 + m * 100 + tp->day <= 15821004 ) {
	    /* Julian before 1582 */
	    b = -2 + ((y+4716)/4) - 1179;	/* Julian calendar */
	} else {
	    /* Gregorian after 1582, We do this */
	    b = y/400 - y/100 + y/4;		/* Gregorian calendar */
	}

	a = 30.6001*(m+1);
	mjd = 365 * y -679004 + b + a + tp->day;
	return mjd;
}

void
set_day ( struct day *tp, int y, int m, int d )
{
	tp->year = y;
	if ( is_leap ( y ) )
	    tp->dpm = dpm_leap;
	else
	    tp->dpm = dpm_normal;
	tp->month = m;
	tp->days_in_month = tp->dpm[tp->month-1];
	tp->day = d;
	tp->mjd0 = calc_mjd ( tp );
	tp->jd0 = tp->mjd0 + MJD_OFFSET;
	// printf ( "Set time, JD = %.5f\n", tp->jd0 );
}

/* Iterate through all the days in a year.
 * - this is happy to run into the next year!
 * (not what you usually want).
 */
void
next_day ( struct day *tp )
{
	if ( tp->day >= tp->days_in_month ) {
	    if ( tp->month == 12 ) {
		++tp->year;
		if ( is_leap ( tp->year ) )
		    tp->dpm = dpm_leap;
		else
		    tp->dpm = dpm_normal;
		tp->month = 1;
		tp->days_in_month = tp->dpm[tp->month-1];
		tp->day = 1;
		tp->mjd0 = calc_mjd ( tp );
		tp->jd0 = tp->mjd0 + MJD_OFFSET;
		return;
	    }
	    ++tp->month;
	    tp->days_in_month = tp->dpm[tp->month-1];
	    tp->day = 1;
	    tp->mjd0 = calc_mjd ( tp );
	    tp->jd0 = tp->mjd0 + MJD_OFFSET;
	    return;
	}
	++tp->day;
	tp->mjd0 = calc_mjd ( tp );
	tp->jd0 = tp->mjd0 + MJD_OFFSET;
	return;
}

/* Iterate through all the days in a year.
 * returns 0 at end of year.
 * Note that even when it returns 0, it has properly
 * set everything up for the next year.
 * This allows the structure to be used to calculate
 *  "next day" values for the last day in December.
 */
int
next_day_y ( struct day *tp )
{
	if ( tp->day >= tp->days_in_month ) {
	    if ( tp->month == 12 ) {
		++tp->year;
		if ( is_leap ( tp->year ) )
		    tp->dpm = dpm_leap;
		else
		    tp->dpm = dpm_normal;
		tp->month = 1;
		tp->days_in_month = tp->dpm[tp->month-1];
		tp->day = 1;
		tp->mjd0 = calc_mjd ( tp );
		tp->jd0 = tp->mjd0 + MJD_OFFSET;
		return 0;
	    }
	    ++tp->month;
	    tp->days_in_month = tp->dpm[tp->month-1];
	    tp->day = 1;
	    tp->mjd0 = calc_mjd ( tp );
	    tp->jd0 = tp->mjd0 + MJD_OFFSET;
	    return 1;
	}
	++tp->day;
	tp->mjd0 = calc_mjd ( tp );
	tp->jd0 = tp->mjd0 + MJD_OFFSET;
	return 1;
}

#ifdef notdef
int
next_day_y ( struct day *tp )
{
	if ( tp->day >= tp->days_in_month ) {
	    if ( tp->month == 12 )
		return 0;
	    ++tp->month;
	    tp->days_in_month = tp->dpm[tp->month-1];
	    tp->day = 1;
	    tp->mjd0 = calc_mjd ( tp );
	    tp->jd0 = tp->mjd0 + MJD_OFFSET;
	    return 1;
	}
	++tp->day;
	tp->mjd0 = calc_mjd ( tp );
	tp->jd0 = tp->mjd0 + MJD_OFFSET;
	return 1;
}
#endif

void
prior_day ( struct day *cur, struct day *prior )
{
	// structure copy !
	*prior = *cur;

	if ( prior->day == 1 ) {
	    if ( prior->month == 1 ) {
		--prior->year;
		if ( is_leap ( prior->year ) )
		    prior->dpm = dpm_leap;
		else
		    prior->dpm = dpm_normal;
		prior->month = 12;
		prior->days_in_month = prior->dpm[prior->month-1];
		prior->day = prior->days_in_month;
	    } else {
		--prior->month;
		prior->days_in_month = prior->dpm[prior->month-1];
		prior->day = prior->days_in_month;
	    }
	} else
	    --prior->day;

	prior->mjd0 = calc_mjd ( prior );
	prior->jd0 = prior->mjd0 + MJD_OFFSET;
	return;
}

/* Iterate through all the days in a month.
 * returns 0 at end of month.
 */
int
next_day_m ( struct day *tp )
{
	if ( tp->day >= tp->days_in_month )
	    return 0;
	++tp->day;

	tp->mjd0 = calc_mjd ( tp );
	tp->jd0 = tp->mjd0 + MJD_OFFSET;
	return 1;
}

/* ---------------------------------------------------------- */


#define MAX_DAYS	366

struct ephem_data {
	int month;
	int day;
	double st_midnight;
	double sun_rise;
	double sun_set;
	double moon_rise;
	double moon_set;
	double moon_age;
	double civil_rise;
	double naut_rise;
	double astro_rise;
	double civil_set;
	double naut_set;
	double astro_set;
};

struct ephem_data ephem_info[MAX_DAYS];

/* Let's test our mjd routine.
 * I have a bunch of JD values from the Astronomical Almanac
 *  pages B12 or B13 (depending on the year).
 * These values are for January 0
 *
 * The ancient fortran MMT almanac program used its own
 *  flavor of MJD (but calling it JD) which was referenced to
 *  the JD_1900 date.
 */

static void
test_jd_one ( int year, double aa_jd )
{
	struct day now;
	double err;

	set_day ( &now, year, 1, 0 );

	err = aa_jd - (now.mjd0 + MJD_OFFSET);
	if ( fabs(err) < 0.0001 )
	    printf ( "jd (%d) = %.2f (%.2f) OK\n", year, now.mjd0 + MJD_OFFSET, aa_jd );
	else
	    printf ( "jd (%d) = %.2f (%.2f) Error !!\n", year, now.mjd0 + MJD_OFFSET, aa_jd );
}

void
test_jd ( void )
{
	test_jd_one ( 1900, JD_1900 - 0.5 );
	test_jd_one ( 2007, 2454100.5 );
	test_jd_one ( 2009, 2454831.5 );
	test_jd_one ( 2010, 2455196.5 );
	test_jd_one ( 2011, 2455561.5 );
	test_jd_one ( 2012, 2455926.5 );
	test_jd_one ( 2013, 2456292.5 );
	test_jd_one ( 2014, 2456657.5 );
	test_jd_one ( 2015, 2457022.5 );
	test_jd_one ( 2016, 2457387.5 );
	test_jd_one ( 2017, 2457753.5 );
	test_jd_one ( 2018, 2458118.5 );
}

/* And let's test our sidereal time calculation.
 * Once again, I have a bunch of ST values from the Astronomical Almanac.
 * I looked up  Sideral Time at Midnight, Jan 1.
 *        almanac page B12 in 2006
 *        almanac page B13 in 2009, ... 2016
 *        G. Sideral Time for Jan 1, 0h UT
 * Note that these are NOT local, but at Greenwich.
 */
static void
test_st_one ( int year, int h, int m, double s )
{
	struct day now;
	double st_apc;
	double st_aa;
	double st_mmt;
	char apc_buf[32];
	char mmt_buf[32];

	st_aa = hms ( h, m, s );

	set_day ( &now, year, 1, 1 );
	st_apc = apc_st ( now.mjd0 ) * RADHR;
	st_mmt = mmt_st ( now.mjd0 + MJD_OFFSET );

	if ( fabs ( st_aa - st_apc ) < 1.0 )
	    printf ( "ST (%d): %s %s (%d:%d:%.3f) OK\n", year, s_dms(apc_buf,st_apc), s_dms(mmt_buf,st_mmt), h, m, s );
	else
	    printf ( "ST (%d): %s %s (%d:%d:%.3f) -- ERROR !!\n", year, s_dms(apc_buf,st_apc), s_dms(mmt_buf,st_mmt), h, m, s );
}

void
test_st ( void )
{
	test_st_one ( 2005, 6, 42, 58.4748 );
	test_st_one ( 2006, 6, 42, 1.5159 );
	test_st_one ( 2007, 6, 41, 4.5504 );
	test_st_one ( 2008, 6, 40, 7.5881 );
	test_st_one ( 2009, 6, 43, 7.1394 );
	test_st_one ( 2010, 6, 42, 10.0357 );
	test_st_one ( 2011, 6, 41, 12.8088 );
	test_st_one ( 2012, 6, 40, 15.4861 );
	test_st_one ( 2013, 6, 43, 14.6094 );
	test_st_one ( 2014, 6, 42, 17.058 );
	test_st_one ( 2015, 6, 41, 19.4297 );
	test_st_one ( 2016, 6, 40, 21.7893 );
	test_st_one ( 2017, 6, 43, 20.7109 );
	test_st_one ( 2018, 6, 42, 23.1082 );
}

/* Here is the routine in use at the MMT for converting
 * from HA and Dec to Alt/Az
 * All arguments and results in radians.
 */
void
toaltaz (double ha, double dec, double *alt, double *az)
{
        double sinh, cosh, sinz, cosz;
        double sindec, cosdec, sinha, cosha;

        sindec = sin (dec);
        cosdec = cos (dec);
        sinha = sin (ha);
        cosha = cos (ha);

        // sinh = site_info.sinlat * sindec + site_info.coslat * cosdec * cosha;
        sinh = site_info.sin_lat * sindec + site_info.cos_lat * cosdec * cosha;

        *alt = asin (sinh);
        cosh = cos (*alt);

        if (fabs(cosh) > 0.0017) {
            sinz = -cosdec * sinha / cosh;
            // cosz = (sindec - site_info.sinlat * sinh) / site_info.coslat /cosh;
            cosz = (sindec - site_info.sin_lat * sinh) / site_info.cos_lat /cosh;
            *az = atan2 (sinz, cosz);
        } else
            *az = PI;

        if (*az < 0.0)
            *az += TWOPI;
}

/* This takes ha in hours, dec in degrees.
 * yields alt/az in degrees, 0.0
 */
void
hd_to_aa ( double ha, double dec, double *alt, double *az )
{
	double sinalt;

	sinalt = site_info.sin_lat * sin(dec*DEGRAD) + site_info.cos_lat * cos(dec*DEGRAD) * cos(ha*HRRAD);
	*alt = asin ( sinalt ) * RADDEG;
	*az = 0.0;	/* XXX */
}

struct sun_info {
	double alt;
	double lst;
};

/* Returns the suns elevation above the horizon in degrees */
static struct sun_info *
sun_alt ( struct day *now, double hour, int verbose )
{
	double jd;
	double lst;
	double ha;
	char buf[32];
	double mjd;
	double t;
	double ra, dec;
	double sinalt;
	double alt, az;
	static struct sun_info info;

	// printf ( " MJD = %.2f\n", now->mjd );
	/*
	 * On October 4, using the MMT mount computer, I get:
	 *  time WED OCT 04 16:22:50 2017
	 *  uttime 23:22:50
	 *  jd 2458031.47420
	 *  mjd 58030.97420
	 *  lst 16:54:38.872
	 * This code calculates (given this JD):
	 *  LST = 16:54:38.779
	 */

	/* The hour argument is a local time, but we want UT to
	 * add to mjd0 to get mjd, so adjust by the time zone
	 */
	mjd = now->mjd0 + (hour - site_info.tz) / 24.0;
	jd = mjd + MJD_OFFSET;
	// jd = 2458031.47420;
	// lst =  mmt_st ( jd ) - site_info.long_hours;
	lst =  mmt_st ( jd ) + site_info.long_hours;
	if ( lst < 0.0 ) lst += 24.0;
	if ( lst > 24.0 ) lst -= 24.0;
	// printf ( "Hour: %6.2f, JD = %.2f, LST = %s\n", hour, jd, s_dms(buf,lst) );
	/* OK to here XXX XXX */

	sun_apc ( mjd, &ra, &dec, verbose );

	ha = lst - ra;
	// printf ( "Hour: %6.2f, JD = %.2f, LST = %s HA = %.2f\n", hour, jd, s_dms(buf,lst), ha );

	/* Convert to alt, az */
	hd_to_aa ( ha, dec, &alt, &az );

	//sinalt = site_info.sin_lat * sin(dec*DEGRAD) + site_info.cos_lat * cos(dec*DEGRAD)*cos(ha*HRRAD);
	// alt = asin ( sinalt ) * RADDEG;

	// if ( verbose && fabs(alt) < 1.0 )
	if ( verbose  )
	    printf ( "Hour: %10.4f, JD = %.2f, LST = %s RA = %.3f HA = %.2f -- alt: %.3f\n", hour, jd, s_dms(buf,lst), ra, ha, alt );

	info.alt = alt;
	info.lst = lst;
	return &info;
	// return alt;
}

/* Returns the elevation of the moon above the horizon in degrees */
static double
moon_alt ( struct day *now, double hour, int verbose )
{
	double jd;
	double lst;
	double ha;
	char buf[32];
	double mjd;
	double t;
	double ra, dec;
	double sinalt;
	double alt, az;

	/* The hour argument is a local time, but we want UT to
	 * add to mjd0 to get mjd, so adjust by the time zone
	 */
	mjd = now->mjd0 + (hour - site_info.tz) / 24.0;
	jd = mjd + MJD_OFFSET;

	lst =  mmt_st ( jd ) + site_info.long_hours;
	if ( lst < 0.0 ) lst += 24.0;
	if ( lst > 24.0 ) lst -= 24.0;

	moon_apc ( mjd, &ra, &dec );

	ha = lst - ra;

	/* Convert to alt, az */
	hd_to_aa ( ha, dec, &alt, &az );

	if ( verbose  )
	    printf ( "Hour: %10.4f, JD = %.2f, LST = %s RA = %.3f HA = %.2f -- alt: %.3f\n", hour, jd, s_dms(buf,lst), ra, ha, alt );

	return alt;
}

/* Quadratic interpolation - APC pages 48-49 */
static int
quad ( double ya, double yb, double yc, double *xe, double *ye, double *r1, double *r2 )
{
	double a, b, c;
	double disc;
	double dx;
	double lxe;
	int nr = 0;

	a = 0.5 * ( ya + yc ) - yb;
	b = 0.5 * ( yc - ya );
	c = yb;

	*xe = lxe = -b / (2.0 * a);
	*ye = (a*lxe + b) * lxe + c;
	disc = b * b - 4.0 * a * c;

	if ( disc < 0.0 )
	    return 0;

	dx = 0.5 * sqrt(disc) / fabs(a);
	*r1 = lxe - dx;
	*r2 = lxe + dx;
	if ( fabs(*r1) <= 1.0 ) ++nr;
	if ( fabs(*r2) <= 1.0 ) ++nr;
	if ( *r1 < -1.0 ) *r1 = *r2;
	return nr;
}

struct sun_event_info {
	double lt_rise;		/* local time of rise */
	double lt_set;		/* local time of set */
	double lst_rise;	/* LST of rise */
	double lst_set;		/* LST of set */
	int rises;
	int sets;
	int above;
};

/* Find the sunrise and sunset times for a given day.
 */
static struct sun_event_info *
sun_events ( struct day *now, double horizon )
/*
	double *lt_rise, double *lt_set,
	int *arises, int *asets, int *aabove )
    */
{
	double ya, yb, yc;
	double hour;
	double xe, ye, r1, r2;
	int nr;
	static struct sun_event_info info;
	struct sun_info *si;

	info.rises = 0;
	info.sets = 0;
	info.above = 0;

	hour = 1.0;

	si = sun_alt ( now, hour - 1.0, 0 );
	ya = si->alt - horizon;
	if ( ya > 0.0 )
	    info.above = 1;

	do {
	    si = sun_alt ( now, hour, 0 );
	    yb = si->alt - horizon;
	    si = sun_alt ( now, hour + 1.0, 0 );
	    yc = si->alt - horizon;

	    nr = quad ( ya, yb, yc, &xe, &ye, &r1, &r2 );
	    // printf ( "abc = %.3f %.3f %.3f  %d\n", ya, yb, yc, nr );
	    if ( nr == 1 ) {
		if ( ya < 0.0 ) {
		    info.lt_rise = hour + r1;
		    info.rises = 1;
		    si = sun_alt ( now, hour + r1, 0 );
		    info.lst_rise = si->lst;
		} else {
		    info.lt_set = hour + r1;
		    info.sets = 1;
		    si = sun_alt ( now, hour + r1, 0 );
		    info.lst_set = si->lst;
		}
	    }
	    if ( nr == 2 ) {
		if ( ye < 0.0 ) {
		    info.lt_rise = hour + r2;
		    info.lt_set = hour + r1;
		    si = sun_alt ( now, hour + r2, 0 );
		    info.lst_rise = si->lst;
		    si = sun_alt ( now, hour + r1, 0 );
		    info.lst_set = si->lst;
		} else {
		    info.lt_rise = hour + r1;
		    info.lt_set = hour + r2;
		    si = sun_alt ( now, hour + r1, 0 );
		    info.lst_rise = si->lst;
		    si = sun_alt ( now, hour + r2, 0 );
		    info.lst_set = si->lst;
		}
		info.rises = 1;
		info.sets = 1;
	    }

	    ya = yc;
	    hour += 2.0;
	} while ( hour < 25.0 && info.rises + info.sets < 2 );

	return &info;
}

/* Find the sunrise and sunset times for a given day.
 * Same code as above with different xx_alt() routine.
 * Also, we scan from 12 to 36 hours, to get afternoon and
 * evening ris/set events.
 */
static void
moon_events ( struct day *now, double horizon,
	double *lt_rise, double *lt_set,
	int *arises, int *asets, int *aabove )
{
	double ya, yb, yc;
	double hour;
	double xe, ye, r1, r2;
	int nr;
	int rises, sets, above;

	rises = 0;
	sets = 0;
	above = 0;

	/* With these, the above flags may be unnecessary */
	*lt_rise = 99.0;
	*lt_set = 99.0;

	hour = 13.0;

	ya = moon_alt ( now, hour - 1.0, 0 ) - horizon;
	if ( ya > 0.0 )
	    above = 1;

	do {
	    yb = moon_alt ( now, hour, 0 ) - horizon;
	    yc = moon_alt ( now, hour + 1.0, 0 ) - horizon;
	    nr = quad ( ya, yb, yc, &xe, &ye, &r1, &r2 );
	    // printf ( "abc = %.3f %.3f %.3f  %d\n", ya, yb, yc, nr );
	    if ( nr == 1 ) {
		if ( ya < 0.0 ) {
		    *lt_rise = hour + r1;
		    rises = 1;
		} else {
		    *lt_set = hour + r1;
		    sets = 1;
		}
	    }
	    if ( nr == 2 ) {
		if ( ye < 0.0 ) {
		    *lt_rise = hour + r2;
		    *lt_set = hour + r1;
		} else {
		    *lt_rise = hour + r1;
		    *lt_set = hour + r2;
		}
		rises = 1;
		sets = 1;
	    }

	    ya = yc;
	    hour += 2.0;
	} while ( hour < 37.0 && rises + sets < 2 );

	if ( rises ) {
	    if (*lt_rise > 24.0 )
		*lt_rise -= 24.0;
	}

	if ( sets ) {
	    if (*lt_set > 24.0 )
		*lt_set -= 24.0;
	}

	*arises = rises;
	*asets = sets;
	*aabove = above;
}


/* Horizon target values in degrees.
 * We call it sunrise (for example) when the calculated geocentric
 *  altitude of the sun matches this target value.
 * The Sun and Moon values are the angular half size
 *  corrected for parallax and refraction.
 * The apparent angular size of the Sun varies (due to the eccentricity
 *  of the earths orbit from 0.524 to 0.542 degrees.
 *  A value of 0.53 is often used.
 *  The value 50/60 is 0.8333 and is the apparent radius (0.265)
 *  along with a ballpark value for parallax and refraction.
 *  This is a hoary historical value used by tradition.
 *
 * For the purposes of civil, nautical, and astronomical twilight,
 *  we don't fret over the correction for refraction.
 *
 *  The value for the moon is curious.
 *  The moons apparent diameter is only slightly less than that of the sun
 *  (roughly 0.5 degrees), so this value of 8/60 positive is another hoary
 *  traditional value that I have no explanation for thus far.
 *
 * Note that the 0.8333 degree offset translates to 3.33 minutes of time.
 */
#define SUN_HORIZON	(-50.0/60.0)
#define MOON_HORIZON	(8.0/60.0)

#define CIVIL_HORIZON	-6.0
#define NAUT_HORIZON	-12.0
#define ASTRO_HORIZON	-18.0

static void
test_sun_one ( struct day *now, double hour, int verbose )
{
	    (void) sun_alt ( now, hour, verbose );
}

/* Generate table of sun altitudes by hour */
void
test_sun1 ( void )
{
	double hour;
	struct day now;

	// set_day ( &now, 2017, 10, 5 );
	set_day ( &now, 2017, 10, 4 );

	for ( hour = 0.0; hour < 23.5; hour += 1.0 ) {
	    test_sun_one ( &now, hour, 1 );
	}
}

/* Generate table of sun altitudes by minute */
void
test_sun2 ( void )
{
	double hour;
	struct day now;
	double del = 1.0 / (24.0 * 60.0 );

	set_day ( &now, 2017, 10, 4 );

	for ( hour = 0.0; hour < 24.0; hour += del ) {
	    test_sun_one ( &now, hour, 1 );
	}
}

/* Here the aim is to get a correct sunrise/sunset time for a specific date
 */

#define SUN3A_YEAR	2017
#define SUN3A_MONTH	10
#define SUN3A_DAY	4

/* From NOAO website - 
 * https://www.esrl.noaa.gov/gmd/grad/solcalc/sunrise.html
 * using 32:15:46 and 11:02:55
 * The following are times in hours.
 * Note that:
 *	1 minute = .01667 hours
 *	1 second = .00027778 hours
 */
#define SUN3A_RISE	6.333	/* 6:20 */
#define SUN3A_SET	18.0833	/* 6:05 */

#define SUN3_DEL	0.02	/* a bit over 1 minute */

void
test_sun3a ( void )
{
	struct day now;
	double rise_hour;
	double set_hour;
	// int rises, sets, above;
	double horizon;
	char buf[32];
	struct sun_event_info *se;

	init_site ( &site_info, &site_castellon );
	set_day ( &now, SUN3A_YEAR, SUN3A_MONTH, SUN3A_DAY );

	horizon = SUN_HORIZON;

	// sun_events ( &now, horizon, &rise_hour, &set_hour, &rises, &sets, &above );
	se = sun_events ( &now, horizon );
	rise_hour = se->lt_rise;
	set_hour = se->lt_set;

/*
 *  I get:
 *	Sun rise 6.339 6:20:21.237
 *	Sun  set 18.077 18:4:38.193
 */

	printf ( "\n" );
	printf ( "Location: %s\n", site_info.name );
	printf ( "Day: %d-%d-%d\n", SUN3A_MONTH, SUN3A_DAY, SUN3A_YEAR );

	if ( fabs(rise_hour-SUN3A_RISE) < SUN3_DEL )
	    printf ( "Sun rise %.3f %s OK\n", rise_hour, s_dms(buf,rise_hour) );
	else
	    printf ( "Sun rise %.3f %s --- ERR\n", rise_hour, s_dms(buf,rise_hour) );

	if ( fabs(set_hour-SUN3A_SET) < SUN3_DEL )
	    printf ( "Sun set %.3f %s OK\n", set_hour, s_dms(buf,set_hour) );
	else
	    printf ( "Sun set %.3f %s --- ERR\n", set_hour, s_dms(buf,set_hour) );

	// double mjd;
	// mjd = 2458031.47420 - MJD_OFFSET;
	// sun_wiki ( mjd, &xx );
}

#define SUN3B_YEAR	2017
#define SUN3B_MONTH	12
#define SUN3B_DAY	1

/* These are from the Fortran code MMT almanac program */
#define SUN3B_RISE	7.11667	/* 7:07 */
#define SUN3B_SET	17.3333	/* 17:20 */

void
test_sun3b ( void )
{
	struct day now;
	double rise_hour;
	double set_hour;
	// int rises, sets, above;
	double horizon;
	char buf[32];
	struct sun_event_info *se;

	init_site ( &site_info, &site_mmt );
	set_day ( &now, SUN3B_YEAR, SUN3B_MONTH, SUN3B_DAY );

	horizon = SUN_HORIZON;

	// sun_events ( &now, horizon, &rise_hour, &set_hour, &rises, &sets, &above );
	se = sun_events ( &now, horizon );
	rise_hour = se->lt_rise;
	set_hour = se->lt_set;

	printf ( "\n" );
	printf ( "Location: %s\n", site_info.name );
	printf ( "Day: %d-%d-%d\n", SUN3B_MONTH, SUN3B_DAY, SUN3B_YEAR );

	if ( fabs(rise_hour-SUN3B_RISE) < SUN3_DEL )
	    printf ( "Sun rise %.3f %s OK\n", rise_hour, s_dms(buf,rise_hour) );
	else
	    printf ( "Sun rise %.3f %s --- ERR\n", rise_hour, s_dms(buf,rise_hour) );

	if ( fabs(set_hour-SUN3B_SET) < SUN3_DEL )
	    printf ( "Sun set %.3f %s OK\n", set_hour, s_dms(buf,set_hour) );
	else
	    printf ( "Sun set %.3f %s --- ERR\n", set_hour, s_dms(buf,set_hour) );

	// double mjd;
	// mjd = 2458031.47420 - MJD_OFFSET;
	// sun_wiki ( mjd, &xx );
}

static void
rise_set ( struct day *now, double h, double *rh, double *sh, double *lst_r, double *lst_s )
{
	// int rises, sets, above;
	struct sun_event_info *se;

	// sun_events ( now, h, rh, sh, &rises, &sets, &above );
	se = sun_events ( now, h );
	*rh = se->lt_rise;
	*sh = se->lt_set;
	*lst_r = se->lst_rise;
	*lst_s = se->lst_set;
}

static void
rise_set_moon ( struct day *now, double *rh, double *sh )
{
	int rises, sets, above;

	moon_events ( now, MOON_HORIZON, rh, sh, &rises, &sets, &above );
	// moon_events ( now, SUN_HORIZON, rh, sh, &rises, &sets, &above );
}

#ifdef notdef
# ----> ANEWM - First new moon at Greenwich.
#    almanac, page A1
#    subtract 7h for longitude.
#    subtract 1d for midnight.
# (2008) newmoon = dhm( 8.0, 11.0, 37.0 ) - dhm( 1.0, timezone, 0.0 )
# (2009) newmoon = dhm( 26.0, 7.0, 55.0 ) - dhm( 1.0, timezone, 0.0 )
# (2010) newmoon = dhm( 15.0, 7.0, 11.0 ) - dhm( 1.0, timezone, 0.0 )
# (2011) newmoon = dhm( 4.0, 9.0, 3.0 ) - dhm( 1.0, timezone, 0.0 )
# (2012) newmoon = dhm( 23.0, 7.0, 39.0 ) - dhm( 1.0, timezone, 0.0 )
# (2013) newmoon = dhm( 11.0, 19.0, 44.0 ) - dhm( 1.0, timezone, 0.0 )
# (2014) newmoon = dhm( 1.0, 11.0, 14.0 ) - dhm( 1.0, timezone, 0.0 )
# (2015) newmoon = dhm( 20.0, 13.0, 14.0 ) - dhm( 1.0, timezone, 0.0 )
# (2016) newmoon = dhm( 10.0, 1.0, 31.0 ) - dhm( 1.0, timezone, 0.0 )
# (2017) newmoon = dhm( 28.0, 0.0, 7.0 ) - dhm( 1.0, timezone, 0.0 )
# (2018) newmoon = dhm( 17.0, 2.0, 17.0 ) - dhm( 1.0, timezone, 0.0 )
#endif

static void
show_old_one ( int y, double d, double h, double m )
{
	double val;
	
	h += m / 60.0;
	h -= 7.0;

	val = d - 1.0;
	val += h/24.0;

	printf ( "ANEWM: %d  %10.4f\n", y, val );
}

void
show_old_newm ( void )
{
	show_old_one ( 2008, 8.0, 11.0, 37.0 );
	show_old_one ( 2009, 26.0, 7.0, 55.0 );
	show_old_one ( 2010, 15.0, 7.0, 11.0 );
	show_old_one ( 2011, 4.0, 9.0, 3.0 );
	show_old_one ( 2012, 23.0, 7.0, 39.0 );
	show_old_one ( 2013, 11.0, 19.0, 44.0 );
	show_old_one ( 2014, 1.0, 11.0, 14.0 );
	show_old_one ( 2015, 20.0, 13.0, 14.0 );
	show_old_one ( 2016, 10.0, 1.0, 31.0 );
	show_old_one ( 2017, 28.0, 0.0, 7.0 );
	show_old_one ( 2018, 17.0, 2.0, 17.0 );
}

#ifdef notdef
/* This demonstrates that the solar longitude is 280 degrees on
 * January 1 of every year, within about +-1 degree
 */
void
xx_solar ( void )
{
	struct day now;
	double s_lat, s_long;
	int y;

	init_site ( &site_info, &site_mmt );

	for ( y=2000; y<2020; y++ ) {
	    set_day ( &now, y, 1, 1 );
	    sun_apc_ll ( now.mjd0, &s_lat, &s_long, 0 );
	    printf ( "Sun long (%d) = %.5f\n", y, s_long * RADDEG );
	}
}
#endif

/* The ANEW_STR value is from the Blue book (Astronomical Almanac)
 * as per these comments in my old mkdata program:
 * These values were important for verification and testing
 * of this program, but are now only of historic interest.
 *
# ----> ANEWM - First new moon at Greenwich.
#    almanac, page A1
#    subtract 7h for longitude.
#    subtract 1d for midnight.
# (2008) newmoon = dhm( 8.0, 11.0, 37.0 ) - dhm( 1.0, timezone, 0.0 )
# (2009) newmoon = dhm( 26.0, 7.0, 55.0 ) - dhm( 1.0, timezone, 0.0 )
# (2010) newmoon = dhm( 15.0, 7.0, 11.0 ) - dhm( 1.0, timezone, 0.0 )
# (2011) newmoon = dhm( 4.0, 9.0, 3.0 ) - dhm( 1.0, timezone, 0.0 )
# (2012) newmoon = dhm( 23.0, 7.0, 39.0 ) - dhm( 1.0, timezone, 0.0 )
# (2013) newmoon = dhm( 11.0, 19.0, 44.0 ) - dhm( 1.0, timezone, 0.0 )
# (2014) newmoon = dhm( 1.0, 11.0, 14.0 ) - dhm( 1.0, timezone, 0.0 )
# (2015) newmoon = dhm( 20.0, 13.0, 14.0 ) - dhm( 1.0, timezone, 0.0 )
# (2016) newmoon = dhm( 10.0, 1.0, 31.0 ) - dhm( 1.0, timezone, 0.0 )
# (2017) newmoon = dhm( 28.0, 0.0, 7.0 ) - dhm( 1.0, timezone, 0.0 )
# (2018) newmoon = dhm( 17.0, 2.0, 17.0 ) - dhm( 1.0, timezone, 0.0 )
newmoon = dhm( 6.0, 1.0, 28.0 ) - dhm( 1.0, timezone, 0.0 )
 */

// #define ANEW_YEAR	2012
// #define ANEW_STR	"23:7:39"
// #define ANEW_YEAR	2017
// #define ANEW_STR	"28:0:7"
#define ANEW_YEAR	2019
// #define ANEW_STR	"6:1:28"

/* The solar longitude increases about 1 degree per day.
 * The lunar longitude increases about 13 degrees per day.
 * On January 1, the sun is at an ecliptic longitude of about 280
 * and will advance to about 311 by the end of the month.
 * For the moon, anything is possible, but we need to go to
 * some pains to ensure a monotonic value.
 */

/* Helper function for anew calculation */
double
anew_diff ( double mjd )
{
	double s_lat, s_long;
	double m_lat, m_long;
	double dif;

	sun_apc_ll ( mjd, &s_lat, &s_long, 0 );
	moon_apc_ll ( mjd, &m_lat, &m_long );

	/* Values in radians */
	// printf ( "MJD %.2f, moon, sun %.2f %.2f\n", mjd, m_long, s_long );

	dif = m_long - s_long;
	if ( dif > PI )
	    dif -= TWOPI;
	if ( dif < -PI )
	    dif += TWOPI;

	return dif;
}

/* Comparing the numbers generated by test_anew() below to
 *  the values in the almanac (using the 2012 copy I have).
 *  I see discrepencies up to 5 minutes using the linear
 *  interpolation across 24 hours.
 * Compare results to values in the table "Phenomena" on page A1
 *  This refinement yields improvement, but not an exact match.
 *
 * f is time within the day in units of hours.
 *   This routine returns an adjusted "f"
 * 
 * XXX - caller needs to beware that f could get adjusted into
 *  the prior or the following day in rare cases.
 */
double
refine_anew_internal ( double mjd, double f )
{
	double f1, f2;
	double d1, d2;
	double fnew;

	/* Use 15 minutes on either side */
	f1 = f - 0.25;
	f2 = f + 0.25;

	d1 = anew_diff ( mjd + f1 / 24.0 );
	d2 = anew_diff ( mjd + f2 / 24.0 );

	// printf ( "  %.3f  %.3f\n", d1, d2 );
	// (usually .002, .003)

	fnew = f1 - d1 * (f2 - f1) / (d2 - d1);

	return fnew;
}

double
refine_anew_mjd ( double mjd, double f )
{
	return mjd + refine_anew_internal ( mjd, f ) / 24.0;
}

double
print_anew ( struct day *d, double f )
{
	char buf[32];
	double fnew;

	fnew = refine_anew_internal ( d->mjd0, f );

	printf ( "%d %s %2d -- New moon: %s\n", d->year, month_name[d->month-1], d->day, s_dms(buf,fnew) );
}

#define TWENTY	(20.0 * DEGRAD)

/* This finds ALL the lunations for a year.
 * 
 * for 2019, it prints the first one as follows:
 *   2019 Jan  6 -- New moon: 1:31:56
 * however, the routine below (print_anew_all()) prints this:
 *   2019 Jan  5 -- New moon: 18:31:56 (ANEW =   4.7722)
 *
 * Indeed in local time the new moon is at 18:31 on 1/5
 * and in UT time it is 1:31 on 1/6
 *
 * So, this routine prints times in UT, and that is due to
 * the "site_greenwich" being specified as the argument to
 * init_site()
 */
void
test_anew ( int year )
{
	struct day now;
	struct day prior;
	struct day prior2;
	double dif;
	double last_dif;
	double f;
	char buf[32];
	int first = 1;
	int count = 0;

	// init_site ( &site_info, &site_mmt );
	init_site ( &site_info, &site_greenwich );

	set_day ( &now, year, 1, 1 );

	// printf ( "ANEW search for %d -- %s\n", year, ANEW_STR );
	// printf ( "ANEW search for %d\n", year );
	printf ( "Lunations for %d\n", year );

	for ( ;; ) {
	    dif = anew_diff ( now.mjd0 );
		
	    // printf ( "%s %2d  %.5f", month_name[now.month-1], now.day, ep->st_midnight );
	    /*
	    printf ( "%s %2d", month_name[now.month-1], now.day );
	    printf ( "  %s (%10.4f)", s_dms(buf,s_long * RADDEG), s_long * RADDEG );
	    printf ( "  %s (%10.4f)", s_dms(buf,m_long * RADDEG), m_long * RADDEG );
	    printf ( "  %10.4f", dif * RADDEG );
	    */
	    if ( ! first && dif > 0.0 && last_dif < 0.0 ) {
		// printf ( " ***\n" );
		prior_day ( &now, &prior );
		f = -last_dif / (dif - last_dif) * 24.0 + site_info.tz;
		if ( f < 0.0 ) {
		    prior_day ( &prior, &prior2 );
		    f += 24.0;
		    print_anew ( &prior2, f );
		    // printf ( "%d %s %2d -- New moon: %s\n", year, month_name[prior2.month-1], prior2.day, s_dms(buf,f) );
		} else {
		    print_anew ( &prior, f );
		    // printf ( "%d %s %2d -- New moon: %s\n", year, month_name[prior.month-1], prior.day, s_dms(buf,f) );
		}
		++count;
	    } else {
		// printf ( "\n" );
	    }

	    first = 0;
	    last_dif = dif;
	    // m_long_last = m_long;
	    // if ( ! next_day_m ( &now ) )
	    if ( ! next_day_y ( &now ) )
		break;
	}

	printf ( "%d lunations\n", count );
}

/* Here are the UT Times for all the lunations in 2019 as yielded by the above.
2019 Jan  6 -- New moon: 1:31:56
2019 Feb  4 -- New moon: 21:08:39
2019 Mar  6 -- New moon: 16:06:54
2019 Apr  5 -- New moon: 8:51:08
2019 May  4 -- New moon: 22:49:23
2019 Jun  3 -- New moon: 10:04:05
2019 Jul  2 -- New moon: 19:18:13
2019 Aug  1 -- New moon: 3:12:29
2019 Aug 30 -- New moon: 10:37:26
2019 Sep 28 -- New moon: 18:30:29
2019 Oct 28 -- New moon: 3:46:01
2019 Nov 26 -- New moon: 15:13:57
2019 Dec 26 -- New moon: 5:19:04
 */

/* Note that the Synodic month (the time between "lunations" varies.
 * The time can vary from 29.18 to 29.93 days because the speed of the earth
 * around the sun varies in the earth's eliptical orbit.
 * An average value of 29.53 days is often quoted.
 * One spot check (the first two in 2012) yielded 29.62 days and I was
 * concerned that this did not match 29.53 until I learned more.
 */

/* This walks forward to the next new moon from the given day
 */

#ifdef notdef
/* XXX - now completely broken due to redefinitions of refine_anew()
 * Active code now uses next_newm_mjd()
 */
double
// next_newm ( int year, int month, int day )
next_newm ( struct day *arg )
{
	struct day now;
	struct day prior;
	struct day prior2;
	double dif;
	double last_dif;
	double f;
	char buf[32];
	int first = 1;
	double val;

	// init_site ( &site_info, &site_mmt );

	// set_day ( &now, year, month, day );
	now = *arg;	/* structure assignment */

	// printf ( "ANEW search for %d -- %s\n", year, ANEW_STR );
	// printf ( "ANEW search for %d\n", year );
	// printf ( "Lunations for %d\n", year );
	printf ( "Next new moon after: %s %d, %d\n", month_name[now.month-1], now.day, now.year );

	for ( ;; ) {
	    dif = anew_diff ( now.mjd0 );
		
	    // printf ( "%s %2d  %.5f", month_name[now.month-1], now.day, ep->st_midnight );
	    /*
	    printf ( "%s %2d", month_name[now.month-1], now.day );
	    printf ( "  %s (%10.4f)", s_dms(buf,s_long * RADDEG), s_long * RADDEG );
	    printf ( "  %s (%10.4f)", s_dms(buf,m_long * RADDEG), m_long * RADDEG );
	    printf ( "  %10.4f", dif * RADDEG );
	    */
	    if ( ! first && dif > 0.0 && last_dif < 0.0 ) {
		// printf ( " ***\n" );
		prior_day ( &now, &prior );
		f = -last_dif / (dif - last_dif) * 24.0 + site_info.tz;
		if ( f < 0.0 ) {
		    prior_day ( &prior, &prior2 );
		    f += 24.0;
		    val = refine_anew ( &prior2, f, 0 );
		    // printf ( "%d %s %2d -- New moon: %s\n", prior2.year, month_name[prior2.month-1], prior2.day, s_dms(buf,f) );
		} else {
		    val = refine_anew ( &prior, f, 0 );
		    // printf ( "%d %s %2d -- New moon: %s\n", prior.year, month_name[prior.month-1], prior.day, s_dms(buf,f) );
		}
		return val;
	    } else {
		// printf ( "\n" );
	    }

	    first = 0;
	    last_dif = dif;
	    // m_long_last = m_long;

	    next_day ( &now );
	}
}
#endif

/* The above routine, but working only with mjd values */
double
next_newm_mjd ( double now )
{
	double prior;
	double prior2;
	double dif;
	double last_dif;
	double f;
	int first = 1;
	double val;

	// printf ( "Next new moon after: %.3f\n", now );

	for ( ;; ) {
	    dif = anew_diff ( now );
		
	    if ( ! first && dif > 0.0 && last_dif < 0.0 ) {
		prior = now - 1.0;
		// prior_day ( &now, &prior );
		f = -last_dif / (dif - last_dif) * 24.0 + site_info.tz;
		if ( f < 0.0 ) {
		    prior2 = prior - 1.0;
		    // prior_day ( &prior, &prior2 );
		    f += 24.0;
		    val = refine_anew_mjd ( prior2, f );
		} else {
		    val = refine_anew_mjd ( prior, f );
		}
		return val;
	    }

	    first = 0;
	    last_dif = dif;

	    now += 1.0;
	}
}

/* Copied and pruned from the above
 * The idea here is to calculate the time of the first new moon
 *  for the specified year.
 * The value returned is in units of days (with 0.x being in January 1)
 * This was used to generate the big anew table given to the MMT in early 2019
 */
void
print_anew_year ( int year )
{
	struct day now;
	struct day prior;
	struct day prior2;
	double dif;
	double last_dif;
	double f;
	char buf[32];
	int first = 1;
	int aday;
	int amon;
	double anew;

	init_site ( &site_info, &site_mmt );

	set_day ( &now, year, 1, 1 );

	for ( ;; ) {
	    dif = anew_diff ( now.mjd0 );

	    if ( ! first && dif > 0.0 && last_dif < 0.0 ) {
		prior_day ( &now, &prior );
		f = -last_dif / (dif - last_dif) * 24.0 + site_info.tz;
		if ( f < 0.0 ) {
		    prior_day ( &prior, &prior2 );
		    f += 24.0;
		    aday = prior2.day;
		    amon = prior2.month-1;

		    // anew = prior2.day - 1.0 + f / 24.0;
		    // printf ( "%d %s %2d -- New moon1: %8s %8.4f\n", year, month_name[prior2.month-1], prior2.day, s_dms(buf,f), anew );
		} else {
		    aday = prior.day;
		    amon = prior.month-1;

		    // anew = prior.day - 1.0 + f / 24.0;
		    // printf ( "%d %s %2d -- New moon2: %8s %8.4f\n", year, month_name[prior.month-1], prior.day, s_dms(buf,f), anew );
		}
		anew = (aday - 1) + f / 24.0;
		printf ( "%d %s %2d -- New moon: %8s (ANEW = %8.4f)\n", year, month_name[amon], aday, s_dms(buf,f), anew );
		// return anew;
		return;
	    }

	    first = 0;
	    last_dif = dif;

	    if ( ! next_day_y ( &now ) )
		break;
	}

	/* i.e. we rolled into the next year */
	printf ( "BIG trouble in %d\n", year );
	// return 0.0;
}

void
anew_table ( void )
{
	int year;

	for ( year = 2000; year < 2100; year++ ) {
	    // calc_anew ( year );
	    print_anew_year ( year );
	}
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

#define OLD_APC
#ifdef OLD_APC
double
Frac ( double x )
{
	return x - floor(x);
}

double
Modulo ( double x, double y )
{
	return y * Frac(x/y);
}
#else
/* linux libraries have these ... */
#define Modulo(x,y)	fmod(x,y)
#define Frac(x,y)	remainder(x,y)
#endif

double
hms ( int h, int m, double s )
{
	double rv;

	rv = (double) m + s / 60.0;
	rv = (double) h + rv / 60.0;
	return rv;
}

/* encode into a d:m:s string
 * (also works fine for h:m:s)
 */
char *
s_dms_delim ( char *buf, double deg, char delim )
{
	int d, m;
	int neg = 0;
	int s;

	if ( deg < 0.0 ) {
	    neg = 1;
	    deg = -deg;
	}

	d = floor ( deg );
	deg -= d;
	deg *= 60.0;
	m = floor ( deg );
	deg -= m;
	deg *= 60.0;

	s = deg + 0.5;


	if ( neg ) {
	    // sprintf ( buf, "-%d%c%02d%c%.3f", d, delim, m, delim, deg );
	    // sprintf ( buf, "-%d%c%02d%c%02d %.3f", d, delim, m, delim, s, deg );
	    sprintf ( buf, "-%2d%c%02d%c%02d", d, delim, m, delim, s );
	} else {
	    // sprintf ( buf, "%d%c%02d%c%.3f", d, delim, m, delim, deg );
	    // sprintf ( buf, "%d%c%02d%c%02d %.3f", d, delim, m, delim, s, deg );
	    sprintf ( buf, " %2d%c%02d%c%02d", d, delim, m, delim, s );
	}
	return buf;

	/*
	printf ( "D = %d\n", d );
	printf ( "M = %d\n", m );
	printf ( "S = %.3f\n", deg );
	*/
}

/* encode into a d:m string
 * (also works fine for h:m)
 */
char *
s_dm_delim ( char *buf, double deg, char delim )
{
	int d, m;
	int neg = 0;
	int s;

	if ( deg < 0.0 ) {
	    neg = 1;
	    deg = -deg;
	}

	d = floor ( deg );
	deg -= d;
	deg *= 60.0;

	m = deg + 0.5;

	if ( neg ) {
	    sprintf ( buf, "-%2d%c%02d", d, delim, m );
	} else {
	    sprintf ( buf, " %2d%c%02d", d, delim, m );
	}
	return buf;
}


char *
s_dms ( char *buf, double deg )
{
	return s_dms_delim ( buf, deg, ':' );
}

char *
s_dms_b ( char *buf, double deg )
{
	return s_dms_delim ( buf, deg, ' ' );
}

char *
s_dm ( char *buf, double deg )
{
	return s_dm_delim ( buf, deg, ':' );
}

char *
s_dm_b ( char *buf, double deg )
{
	return s_dm_delim ( buf, deg, ' ' );
}

void
init_site ( struct site *sp, struct site_data *dp )
{
	sp->lat_deg = dp->lat_deg;
	sp->long_deg = dp->long_deg;
	sp->long_hours = sp->long_deg * 24.0 / 360.0;

#ifdef notdef
	sp->long_hours = LONGITUDE;	/* hours west */
	sp->long_deg = sp->long_hours * -360.0 / 24.0;	/* degrees E */
#endif

	sp->latitude = sp->lat_deg * DEGRAD;
	sp->longitude = sp->long_deg * DEGRAD;

	sp->elevation = dp->elevation;	/* site elevation, now in meters (never used) */
	sp->tz = dp->tz;		/* hours from greenwich, also never used */

	sp->sin_lat = sin(sp->latitude);
	sp->cos_lat = cos(sp->latitude);
	strcpy ( sp->name, dp->name );
}

/* Sidereal time at Greenwich
 * This is from the book "Astronomy on the Personal Computer"
 *  Fourth edition, page 40.
 *  It yields results accurate to about a second, but gives
 *   inferior results compared with the MMT routine.
 */
double
// gmst ( double mjd )
apc_st ( double mjd )
{
	double mjd0;
	double ut, gmst;
	double t, t0;

	mjd0 = floor ( mjd );
	ut = (mjd - mjd0) * SECPD;

	t0 = (mjd0 - 51544.5) / 36525.0;
	t = (mjd - 51544.5) / 36525.0;

	gmst = 24110.54841 + 8640184.812866*t0 + 1.0027379093*ut
	    + (0.093104-6.2e-6*t)*t*t;

	return (TWOPI/SECPD) * Modulo(gmst,SECPD);
}

#define MMT_LST
#ifdef MMT_LST
/* The following is from the mount code used at the MMT Observatory.
 * (as used from 2000 to 10/2017 ...
 * ---
 * routine to compute the local siderial time (in tics)
 * from UT date and time
 * Added DUT fixup 1/6/2002
 */
#define SECPERDAY  	86400.0
#define JD_1970		2440587.5

/* The mount computer at the MMT keeps track of UT time as
 * seconds since 1970, so in theory it could use this
 * interface (but it doesn't).
 * It actually calls something called synch_lst which
 * computes LST based on values in its time_info structure
 * and placing those LST values into the time_info structure.
 * It also keeps these time values as integer seconds along
 * with a fractional value kept as nanosecond counts.
 */

// void synch_lst ( void )
double
mmt_lst_ticks ( double sec )
{
        double utc, ut1;
	double lst_sec;
	double jd;

        /* utc since 1970 in seconds */
        // sec = time_info.ut.tv_sec;
        // sec += time_info.ut.tv_nsec / 1.e9;

        /* Calculate utc for the current day in hours */
        utc = fmod ( sec, SECPERDAY) / 3600.0;
        // ut1 = utc + time_info.dut / 3600.0;
        ut1 = utc;

        /* JD right now */
        jd = sec / SECPERDAY + JD_1970;

        lst_sec = mmt_lst ( jd ) * 3600.0;
        // time_info.lst.tv_sec = lst_sec;
        // time_info.lst.tv_nsec = (lst_sec - (double) time_info.lst.tv_sec) * 1.e9;
	return lst_sec;
}

/* Return the Greenwich ST in hours 
 */
double
mmt_st ( double jd )
{
        double t0, t;
        double jd0, ut1;
	double gmst, omega, e;
	double sec;

	sec = (jd - JD_1970) * SECPERDAY;
        ut1 = fmod ( sec, SECPERDAY ) / 3600.0;
        jd0 = jd - ut1 / 24.0;
	// printf ( "MMT jd, sec, ut1, jd0: %.6f %.6f %.6f %.6f\n", jd, sec, ut1, jd0 );

#ifdef notdef
        /* Calculate utc within the current day in hours
	 *  (i.e. a fraction of a day)
	 */
        utc = fmod ( sec, SECPERDAY ) / 3600.0;
        // ut1 = utc + time_info.dut / 3600.0;
        ut1 = utc;

        /* JD right now */
        jd = sec / SECPERDAY + JD_1970;

        /* JD today at 0h UTC */
        jd0 = jd - ut1 / 24.0;
#endif

        /* 2451545 is JD for Jan 1.5, 2000 */
        t0 = (jd0 - 2451545.0) / 36525.0;
        t  = (jd  - 2451545.0) / 36525.0;

        gmst = 6.69737456 + 2400.051336 * t0 + 2.58622e-5 *
               t0 * t0 + 1.002737909 * ut1;

        omega = 2.182438586 - 33.75704592 * t + 36.145769e-6 * t * t;

        // e = -2.9e-4 * sin (omega);

	/* The following in units of hours */
        // last = gmst - site_info.geodetic_long_h + e;
        // last = gmst - site_info.long_hours + e;
        // last = gmst + site_info.long_hours + e;
        // lst = gmst + e;

        gmst -= 2.9e-4 * sin (omega);

        gmst = fmod (gmst, 24.0);
        if (gmst < 0.0)
            gmst += 24.0;

	return gmst;
}

/* Pass the time as JD at Greenwich in days
 * get back the LST in hours at current location.
 */
double
mmt_lst ( double jd )
{
	double lst;

	// lst =  mmt_st ( jd ) - site_info.long_hours;
	lst =  mmt_st ( jd ) + site_info.long_hours;

        if (lst < 0.0)
            lst += 24.0;
        if (lst > 24.0)
            lst -= 24.0;
	return lst;
}
#endif

/* Obliquity of the ecliptic
 * t = julian centuries from 2000 January  1 (12h)
 * (Julian century is always 36525 days)
 * XXX - should use this instead of static values below.
 */
double
calc_eps ( double jd )
{
	double rv;
	double t;

	t = ( jd - 2451545.0) / 36525.0;
	rv = 23.43929111 - (46.8150 + ( 0.00059 - 0.001813 * t)*t)*t/3600.0;
	return rv * DEGRAD;
	// return 23.43929111 * DEGRAD;
	// return 23.5 * DEGRAD;
}

/* Transform from ecliptic coordinates
 *  to equatorial (ra, dec) coordinates.
 * Starts with lat and long in radians.
 *
 * Yields ra in hours, dec in degrees.
 */
void
ll_to_rd ( double ll_lat, double ll_long, double *ra, double *dec )
{
	double cb, x, y, z, v, w, rho;
	double coseps = 0.91748;	/* XXX */
	double sineps = 0.39778;	/* XXX */

	// Transform to equatorial coordinates
#ifdef MOON_MATRIX
	e_moon = R_x(-eps) * Vec3D(Polar(ll_long,ll_lat));
	*ra = e_moon[phi];
	*dec = e_moon[theta];
#else
	cb = cos ( ll_lat );
	x = cb * cos ( ll_long );
	v = cb * sin ( ll_long );
	w = sin ( ll_lat );

	y = coseps * v - sineps * w;
	z = sineps * v + coseps * w;
	rho = sqrt ( 1.0 - z*z );

	*dec = RADDEG * atan ( z / rho );
	*ra = 2.0 * RADHR * atan ( y / (x + rho) );
	if ( *ra < 0.0 ) *ra += 24.0;
#endif
}

/* From APC, page 38
 * This is a reduced version to calculate lunar position
 *  suitable for moonrise and moonset times
 * From the APC "MiniMoon" routine.
 *
 */
void
// MiniMoon_ll ( double mjd, double *latp, double *longp )
moon_apc_ll ( double mjd, double *latp, double *longp )
{
	double l0, lm, ls, d, f;
	double dl, s, h, n;
	double t;
	// double l_moon, b_moon;

	/*
	 * t = time in julian centuries since J2000
	 */
	t = ( mjd - MJD_2000 ) / JULIAN_CENTURY;

	l0 = Frac ( 0.606433 + 1336.855225 * t);	/* mean longitude (in degrees/360) */

	lm = TWOPI * Frac ( 0.374897 + 1325.552410 * t);	/* moon - mean anomaly */
	ls = TWOPI * Frac ( 0.993133 + 99.997361 * t);		/* sun - mean anomaly */
	d = TWOPI * Frac ( 0.827361 + 1236.853086 * t);	/* diff - long, moon-sun */
	f = TWOPI * Frac ( 0.259086 + 1342.227825 * t);	/* dist from ascending node */

	// Perturbations in long and lat (yields arc seconds)
	dl = 22640.0 * sin(lm) - 4586.0 * sin(lm-2.0*d) + 2370.0 * sin(2.0*d) + 769.0 * sin(2.0*lm)
	    - 668.0*sin(ls) - 412.0 * sin(2.0*f) - 212.0 * sin(2.0*lm-2.0*d) - 206.0 * sin(lm+ls-2.0*d)
	    + 192.0 * sin(lm+2.0*d) - 165.0 * sin(ls-2.0*d) - 125.0 * sin(d) - 110.0 * sin(lm+ls)
	    + 148.0 * sin(lm-ls) - 55.0 * sin(2.0*f-2.0*d);

	s = f + (dl + 412.0*sin(2.0*f) + 541.0 * sin(ls)) * ASEC_TO_RAD;
	h = f - 2.0*d;
	n = -526.0 * sin(h) + 44.0 * sin(lm+h) - 31.0 * sin(-lm+h) -23.0 * sin(ls+h)
	    + 11.0 * sin(-ls+h) - 25.0 * sin(-2.0*lm+f) + 21.0 * sin(-lm + f);

	// l_moon is ecliptic longitude
	// b_moon is ecliptic latitude
	//    This yields both in radians.
	//
	// l_moon = TWOPI * Frac ( l0 + dl / ARCSEC_360 );
	// b_moon = ( 18520.0 * sin(s) + n ) * ASEC_TO_RAD;
	*longp = TWOPI * Frac ( l0 + dl / ARCSEC_360 );
	*latp = ( 18520.0 * sin(s) + n ) * ASEC_TO_RAD;
}

/* Equations from Wikipedia article */
void
sun_wiki ( double mjd, double *rv )
{
	// double jd;
	double n, l, g, g2, ll;

	// jd = 2458031.47420;
	// n = jd - 2451545.0;
	// printf ( "W: %.2f %.2f\n", jd, n );
	// mjd = jd - MJD_OFFSET;

	n = mjd - MJD_2000;
	// printf ( "W: %.2f %.2f\n", mjd, n );

	l = 280.460 + 0.9856474 * n;	/* mean long - degrees */
	g = 357.528 + 0.9856003 * n;	/* anomaly - degrees */

	// These are big values (nominal 6000.0 or so) */
	// printf ( "Check, l, g = %.2f %.2f\n", l, g );
	l = fmod ( l, 360.0 );
	g = fmod ( g, 360.0 );
	// printf ( "Check, l, g = %.2f %.2f\n", l, g );

	g2 = fmod ( g+g, 360.0 );

	ll = l + 1.915 * sin(g*DEGRAD) + 0.020 * sin(g2*DEGRAD);
	// printf ( "Check, solar long (deg) = %.3f\n", ll );

	*rv = ll;
}

/* From APC, page 39  "MiniSun"
 * This is a reduced version to calculate solar positions
 *  adequate for sunrise and sunset times.
 * Note that the sun always is in the ecliptic plane, so the
 *  solar ecliptic latitude is always zero.
 * (it never exceeds 0.00033 degrees).
 *
 * t = time in julian centuries since J2000
 * Returns:
 *  RA in hours
 *  Dec in degrees
 */
void
sun_apc_ll ( double mjd, double *latp, double *longp, int verbose )
{
	double t;
	double m;
	double sun_long;
	char buf[32];

	t = ( mjd - MJD_2000 ) / JULIAN_CENTURY;
	// printf ( "mjd = %.2f, T = %.7f\n", mjd, t );

	/* This gives the suns ecliptic longitude in radians */
	m = TWOPI * Frac ( 0.993133 + 99.997361 * t );
	sun_long = TWOPI * Frac ( 0.7859453 + m / TWOPI +
	    (6893.0 * sin(m) + 72.0 * sin(2.0*m) + 6191.2 * t ) / ARCSEC_360 );

	if ( verbose > 1 ) {
	    printf ( "Sun, solar long (Mini) = %s\n", s_dms(buf,sun_long * RADDEG) );
	}

	*latp = 0.0;
	*longp = sun_long;
}

/* Solar position in equatorial coordinates (ra/dec)
 */
void
sun_apc ( double mjd, double *ra, double *dec, int verbose )
{
	double s_lat, s_long;
	// double w_long;
	char buf[32];

	sun_apc_ll ( mjd, &s_lat, &s_long, verbose );

	// sun_wiki ( mjd, &w_long );
	// printf ( "Sun, solar long (Mini) = %.4f long (Wiki) = %.4f\n", sun_long * RADDEG, w_long );
	// ll_to_rd ( 0.0, w_long * DEGRAD, ra, dec );

	ll_to_rd ( s_lat, s_long, ra, dec );

	if ( verbose > 1 ) {
	    printf ( "Sun, RA  (Mini) = %.4f %s\n", *ra, s_dms(buf,*ra) );
	    printf ( "Sun, Dec (Mini) = %.4f %s\n", *dec, s_dms(buf,*dec) );
	}
}

/* "MiniMoon - ra and dec of the moon at a given time.
 */
void
moon_apc ( double mjd, double *ra, double *dec )
{
	double m_lat, m_long;

	moon_apc_ll ( mjd, &m_lat, &m_long );
	ll_to_rd ( m_lat, m_long, ra, dec );
}

/* The time between new moons is at most 29.93 days, so this should
 * bump us back before the prior new moon.
 */
#define SYNODIC_JOG	31.0

double
moon_age ( struct day *d )
{
	double newm;
	double newm_prior;
	double age;
	double age_prior;

	//printf ( "MJD cur : %14.4f\n", d->mjd0 );
	newm = next_newm_mjd ( d->mjd0 );
	//printf ( "MJD newm: %14.4f\n", newm );
	newm_prior = next_newm_mjd ( newm - SYNODIC_JOG );
	//printf ( "MJD newm: %14.4f\n", newm_prior );

	age = newm - d->mjd0;
	age_prior =  d->mjd0 - newm_prior;

	if ( age < age_prior )
	    return -age;
	else
	    return age_prior;
}

/* ---------------------------------------------------------- */

#ifdef notdef
/* #define EPHEM_YEAR	2017	*/
#define EPHEM_YEAR	2019
#define EPHEM_MONTH	1

/* This actually grinds out the almanac for the year.
 * XXX - we don't really need to put ep in an array.
 */
void
test_almanac ( void )
{
	struct day now;
	struct ephem_data *ep;
	double ut;
	char buf[32];
	int s;

	init_site ( &site_info, &site_mmt );

	printf ( "    %d\n", EPHEM_YEAR );
	// set_day ( &now, EPHEM_YEAR, 1, 1 );
	set_day ( &now, EPHEM_YEAR, EPHEM_MONTH, 1 );

	ep = ephem_info;
	for ( ;; ) {
	    /* We want lst for night after, not night before,
	     *  so start with 24.0
	     */
	    ut = 24.0 - site_info.tz;
	    ep->st_midnight = mmt_lst ( now.mjd0 + MJD_OFFSET + ut / 24.0 );

	    rise_set ( &now, SUN_HORIZON, &ep->sun_rise, &ep->sun_set );
	    // rise_set_moon ( &now, &ep->moon_rise, &ep->moon_set );
	    rise_set ( &now, CIVIL_HORIZON, &ep->civil_rise, &ep->civil_set );
	    rise_set ( &now, NAUT_HORIZON, &ep->naut_rise, &ep->naut_set );
	    rise_set ( &now, ASTRO_HORIZON, &ep->astro_rise, &ep->astro_set );

	    // printf ( "%s %2d  %.5f", month_name[now.month-1], now.day, ep->st_midnight );
	    printf ( "%s %2d", month_name[now.month-1], now.day );
	    printf ( "   %s", s_dm_b(buf,ep->sun_set) );
	    printf ( "   %s", s_dm_b(buf,ep->civil_set) );
	    printf ( "   %s", s_dm_b(buf,ep->naut_set) );
	    printf ( "   %s", s_dm_b(buf,ep->astro_set) );

	    printf ( "  %s", s_dms_b(buf,ep->st_midnight) );

	    printf ( "   %s", s_dm_b(buf,ep->astro_rise) );
	    printf ( "   %s", s_dm_b(buf,ep->naut_rise) );
	    printf ( "   %s", s_dm_b(buf,ep->civil_rise) );
	    printf ( "   %s", s_dm_b(buf,ep->sun_rise) );
	    printf ( "\n" );

	    // s = next_day_y ( &now );
	    s = next_day_m ( &now );
	    if ( s == 0 )
		break;
	    ep++;
	}
}
#endif

/* Main header */
char mh1[] = "                                         Observer's Almanac -- %4d -- MMT Observatory\n";
char mh2[] = "                                 Latitude=+31d 41m 19.6s  Longitude= 110d 53m 04.4s\n";
char mh3[] = "                                           Elevation=2600 m (8585 ft)  Time Zone=+7h\n";
// char mh4[] = "                                                Tom Trebisky, November 27, 2019\n";
char mh4[] = "                                                %s, %s\n";

/* Column headers */
char ch1[] = "     %4d          Sun     Sun     Sun    RA 3H                RA 3H    Sun     Sun     Sun              Moon    Moon    Moon Age\n";
char ch2[] = " Date    Sunset   6 Deg   12 Deg  18 Deg  West at  Sid Time    East at 18 Deg  12 Deg  6 Deg  Sunrise    rise    set       at\n";
char ch3[] = "                  W Hrz   W Hrz   W Hrz   18 Deg   Midnight    18 Deg  E Hrz   E Hrz   E Hrz                            Midnight\n";

// #define FORM_FEED 0x0c

static char person[] = "Tom Trebisky";

static void
print_header ( int year, int first )
{
	char *timestamp;

	timestamp = get_today_str ();

	if ( first )
	    printf ( "\n" );
	else
	    printf ( "\f\n" );


	printf ( "\n" );
	printf ( "\n" );
	printf ( "\n" );
	printf ( mh1, year );
	printf ( mh2 );
	printf ( mh3 );
	printf ( mh4, person, timestamp );
	printf ( "\n" );
	printf ( "\n" );
	printf ( "\n" );
	printf ( "\n" );
	printf ( "\n" );
	printf ( ch1, year );
	printf ( ch2 );
	printf ( ch3 );
	printf ( "\n" );
	printf ( "\n" );
}

/* I personally prefer to see the moonrise/set times even in daylight,
 * to get it that way, I need to set this to zero.
 */
int dark_only = 1;
// int dark_only = 0;

void
generate_almanac ( int year )
{
	struct day now;

	double sun_rise;
	double sun_set;
	double civil_rise;
	double civil_set;
	double naut_rise;
	double naut_set;
	double astro_rise;
	double astro_set;

	double bogus_sun_set;

	double next_sun_rise;
	double next_sun_set;
	double moon_rise;
	double moon_set;

	double lst_rise_bogus;
	double lst_set_bogus;

	double lst_rise_18;
	double lst_set_18;

	double sun_3w;
	double sun_3e;

	double st_midnight;

	char buf[32];
	double ut;
	double age;
	int cur_month = 0;
	int status;
	int first = 1;
	int blank_rise, blank_set;

	init_site ( &site_info, &site_mmt );

	set_day ( &now, year, 1, 1 );

	do {
	    if ( now.month != cur_month ) {
		print_header ( year, first );
		first = 0;
		cur_month = now.month;
	    }

	    /* We want lst for night after, not night before,
	     *  so start with 24.0
	     */
	    ut = 24.0 - site_info.tz;
	    st_midnight = mmt_lst ( now.mjd0 + MJD_OFFSET + ut / 24.0 );

	    /* We will overwrite these rise times below */
	    rise_set ( &now, SUN_HORIZON, &sun_rise, &sun_set, &lst_rise_bogus, &lst_set_bogus );
	    rise_set ( &now, CIVIL_HORIZON, &civil_rise, &civil_set, &lst_rise_bogus, &lst_set_bogus );
	    rise_set ( &now, NAUT_HORIZON, &naut_rise, &naut_set, &lst_rise_bogus, &lst_set_bogus );
	    rise_set ( &now, ASTRO_HORIZON, &astro_rise, &astro_set, &lst_rise_bogus, &lst_set_18 );

	    sun_3w = lst_set_18 - 3.0;
	    if ( sun_3w > 24.0 ) sun_3w -= 24.0;
	    if ( sun_3w < 0.0 ) sun_3w += 24.0;

	    rise_set_moon ( &now, &moon_rise, &moon_set );
	    age = moon_age ( &now );

	    /* Print this info before we bump the day */
	    printf ( "%s %2d", month_name[now.month-1], now.day );

	    /* We need the next day sunrise time to decide if moonrises
	     * and moonsets are in day or night.
	     * As a matter of fact, we need next day sunrise times overall.
	     */
	    status = next_day_y ( &now );
	    rise_set ( &now, SUN_HORIZON, &sun_rise, &bogus_sun_set, &lst_rise_bogus, &lst_set_bogus );
	    rise_set ( &now, CIVIL_HORIZON, &civil_rise, &bogus_sun_set, &lst_rise_bogus, &lst_set_bogus );
	    rise_set ( &now, NAUT_HORIZON, &naut_rise, &bogus_sun_set, &lst_rise_bogus, &lst_set_bogus );
	    rise_set ( &now, ASTRO_HORIZON, &astro_rise, &bogus_sun_set, &lst_rise_18, &lst_set_bogus );
	    next_sun_rise = sun_rise;

	    sun_3e = lst_rise_18 + 3.0;
	    if ( sun_3e > 24.0 ) sun_3e -= 24.0;
	    if ( sun_3e < 0.0 ) sun_3e += 24.0;

	    printf ( "   %s", s_dm_b(buf,sun_set) );
	    printf ( "  %s", s_dm_b(buf,civil_set) );
	    printf ( "  %s", s_dm_b(buf,naut_set) );
	    printf ( "  %s", s_dm_b(buf,astro_set) );
	    printf ( "  %s", s_dm_b(buf,sun_3w) );

	    printf ( "   %s", s_dms_b(buf,st_midnight) );

	    printf ( "   %s", s_dm_b(buf,sun_3e) );
	    printf ( "  %s", s_dm_b(buf,astro_rise) );
	    printf ( "  %s", s_dm_b(buf,naut_rise) );
	    printf ( "  %s", s_dm_b(buf,civil_rise) );
	    printf ( "  %s", s_dm_b(buf,sun_rise) );

	    /* We only print moon rise/set times that happen in dark hours.
	     * Our moon events are found in the interval from this noon to
	     * the next, so we need to check against this days sunset and
	     * the next days sunrise.
	     * Also!  We may not even have a valid rise or set in our
	     * 24 hour period (this doesn't happen for the sun at our latitude.
	     * Check for our 99.0 sentinel for no event.
	     */

	    blank_rise = 0;
	    if ( moon_rise > 50.0 )
		blank_rise = 1;
	    if ( dark_only && moon_rise > next_sun_rise && moon_rise < sun_set )
		blank_rise = 1;

	    if ( blank_rise )
		printf ( "          " );
	    else
		printf ( "    %s", s_dm_b(buf,moon_rise) );

	    blank_set = 0;
	    if ( moon_set > 50.0 )
		blank_set = 1;
	    if ( dark_only && moon_set > next_sun_rise && moon_set < sun_set )
		blank_set = 1;

	    if ( blank_set )
		printf ( "        " );
	    else
		printf ( "  %s", s_dm_b(buf,moon_set) );

	    // printf ( "    %s", s_dm_b(buf,moon_rise) );
	    // printf ( "  %s", s_dm_b(buf,moon_set) );

	    printf ( "   %5.1f", age );

	    printf ( "\n" );

	} while ( status );
}

/* ---------------------------------------------------------- */

#include <time.h>

void
set_today ( struct day *today )
{
	time_t now;
	struct tm *info;

	time( &now );
	info = localtime( &now );

	// printf("Current local time and date: %s", asctime(info));

	set_day ( today, 1900 + info->tm_year, 1 + info->tm_mon, info->tm_mday );
}

char *
get_today_str ( void )
{
	time_t now;
	struct tm *info;
	static char buf[100];

	time( &now );
	info = localtime( &now );

	sprintf ( buf, "%s %d, %d\n", month_name_full[info->tm_mon], info->tm_mday, 1900 + info->tm_year );
	return buf;
}

void
just_today ( int verbose )
{
	struct day now;
	struct ephem_data *ep;
	char buf1[32];
	char buf2[32];
	double age;

	/* values never used */
	double lst_rise;
	double lst_set;

	init_site ( &site_info, &site_mmt );

	// int year = 2019;
	// int month = 12;
	// int day = 9;
	// set_day ( &now, year, month, day );

	set_today ( &now );

	printf ("%s %d, %d\n", month_name_full[now.month-1], now.day, now.year );

	ep = ephem_info;
	rise_set ( &now, SUN_HORIZON, &ep->sun_rise, &ep->sun_set, &lst_rise, &lst_set );
	rise_set_moon ( &now, &ep->moon_rise, &ep->moon_set );

	printf ( "Sunrise: %s  (%s)\n", s_dms(buf1,ep->sun_rise), s_dm(buf2,ep->sun_rise) );
	printf ( "Sunset:  %s  (%s)\n", s_dms(buf1,ep->sun_set), s_dm(buf2,ep->sun_set) );

	printf ( "Moonrise: %s  (%s)\n", s_dms(buf1,ep->moon_rise), s_dm(buf2,ep->moon_rise) );
	printf ( "Moonset:  %s  (%s)\n", s_dms(buf1,ep->moon_set), s_dm(buf2,ep->moon_set) );

	age = moon_age ( &now );
	printf ( "Moon age: %.2f\n", age );
}

/* ---------------------------------------------------------- */

int
main ( int argc, char **argv )
{
	int m;

	// init_site ( &site_info, &site_mmt );
	// init_site ( &site_info, &site_castellon );
	init_site ( &site_info, &site_tucson );

#ifdef notdef
	test_jd ();
	printf ( "\n" );
	test_st ();
	printf ( "\n" );

	// test_sun1 ();
	// printf ( "\n" );
	// test_sun2 ();
	// printf ( "\n" );

	test_sun3a ();
	test_sun3b ();
#endif

#ifdef notdef
	printf ( "\n" );
	test_anew ( 2012 );

	printf ( "\n" );
	test_almanac ();
#endif

#ifdef notdef
	printf ( "\n" );
	test_anew ( 2019 );

	printf ( "\n" );
	show_old_newm ();

	printf ( "\n" );
	anew_table ();

	printf ( "\n" );
	print_anew_all ( 2019 );

	printf ( "\n" );
	test_anew ( 2019 );
#endif

	just_today ( 0 );

	// generate_almanac ( 2019 );
	// generate_almanac ( 2020 );
	// generate_almanac ( 2021 );

	// printf ( "Done\n" );
}

/* THE END */
