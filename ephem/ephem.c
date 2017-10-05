#include <stdio.h>
#include <math.h>

/* (c) Tom Trebisky  10-4-2017
 */

/* Note that the international convention for longitude is positive to the east,
 * so all these longitude coordinates in Arizona ought to be negative.
 * XXX XXX
 */
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

/* Arc seconds in a full circle.
 * i.e. 360*60*60
 */
#define ARCSEC_360	1296000.0

struct time {
	int year;
	int month;
	int day;
	double mjd0;
};

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
};

void init_site ( struct site * );
void set_time ( struct time *, int, int, int );
double apc_st ( double );
char * s_dms ( char *, double );
double hms ( int, int, double );
double mmt_st ( double );
double mmt_lst ( double );

void MiniSun ( double, double *, double * );
void MiniMoon ( double, double *, double * );

void test_jd ( void );
void test_st ( void );
void test_sun1 ( void );
void test_sun2 ( void );

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

struct site site_info;

int
main ( int argc, char **argv )
{
	init_site ( &site_info );

	// test_jd ();
	// printf ( "\n" );
	// test_st ();

	test_sun1 ();
	printf ( "\n" );
	test_sun2 ();

	// printf ( "Done\n" );
}

/* Let's test our mjd routine.
 * I have a bunch of JD values from the Astronomical Almanac
 *  pages B12 or B13 (depending on the year).
 *
 * Note that JD is zero at noon, whereas
 *  MJD is zero at midnight.
 * These are standard conventions.
 *
 * The ancient fortran MMT almanac program used its own
 *  flavor of MJD (but calling it JD) which was referenced to
 *  the JD_1900 date, the 2003 Almanac page K5 says that
 *  1900 January 0 at 0h is JD = 2415019.5
 *  I add the 0.5 to reckon JD from noon as usual.

 */
#define MJD_OFFSET	2400000.5
#define JD_1900		2415020.0

static void
test_jd_one ( int year, double aa_jd )
{
	struct time now;
	double err;

	set_time ( &now, year, 1, 0 );

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
	struct time now;
	double st_apc;
	double st_aa;
	double st_mmt;
	char apc_buf[32];
	char mmt_buf[32];

	st_aa = hms ( h, m, s );

	set_time ( &now, year, 1, 1 );
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

static double
sun_alt ( struct time *now, double hour, int verbose )
{
	double jd;
	double lst;
	double ha;
	char buf[32];
	double mjd;
	double t;
	double ra, dec;
	double sinalt;
	double alt;

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

	mjd = now->mjd0 + (hour + site_info.long_hours) / 24.0;
	jd = mjd + MJD_OFFSET;
	// jd = 2458031.47420;
	lst =  mmt_st ( jd ) - site_info.long_hours;
	if ( lst < 0.0 ) lst += 24.0;
	if ( lst > 24.0 ) lst -= 24.0;
	// printf ( "Hour: %6.2f, JD = %.2f, LST = %s\n", hour, jd, s_dms(buf,lst) );
	/* OK to here XXX XXX */

	t = ( mjd - MJD_2000 ) / JULIAN_CENTURY;
	// printf ( "mjd = %.2f, T = %.7f\n", mjd, t );
	MiniSun ( t, &ra, &dec );

	/* Convert to alt, az */
	ha = lst - ra;
	// printf ( "Hour: %6.2f, JD = %.2f, LST = %s HA = %.2f\n", hour, jd, s_dms(buf,lst), ha );

	sinalt = site_info.sin_lat * sin(dec*DEGRAD) + site_info.cos_lat * cos(dec*DEGRAD)*cos(ha*HRRAD);
	alt = asin ( sinalt ) * RADDEG;
	if ( verbose )
	    printf ( "Hour: %6.2f, JD = %.2f, LST = %s RA = %.3f HA = %.2f -- alt: %.3f\n", hour, jd, s_dms(buf,lst), ra, ha, alt );

	return alt;
}

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

/* Find the sunrise and sunset times for a given day.
 */
static void
sun_events ( struct time *now,
	double *lt_rise, double *lt_set,
	int *rises, int *sets, int *above )
{
	double ya, yb, yc;
	double hour;
	double xe, ye, r1, r2;
	int nr;

	*rises = 0;
	*sets = 0;
	*above = 0;

	ya = sun_alt ( now, hour - 1.0, 0 );
	if ( ya > 0.0 )
	    *above = 1;

	do {
	    yb = sun_alt ( now, hour, 0 );
	    yc = sun_alt ( now, hour + 1.0, 0 );
	    nr = quad ( ya, yb, yc, &xe, &ye, &r1, &r2 );
	    if ( nr == 1 ) {
		if ( ya < 0.0 ) {
		    *lt_rise = hour + r1;
		    *rises = 1;
		} else {
		    *lt_set = hour + r1;
		    *sets = 1;
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
		*rises = 1;
		*sets = 1;
	    }

	    ya = yc;
	    hour += 2.0;
	} while ( hour < 25.0 && *rises == 0 && *sets == 0 );
}

static void
test_sun_one ( struct time *now, double hour )
{
	    (void) sun_alt ( now, hour, 1 );
}

void
test_sun1 ( void )
{
	double hour;
	struct time now;

	// set_time ( &now, 2017, 10, 5 );
	set_time ( &now, 2017, 10, 4 );

	for ( hour = 0.0; hour < 23.5; hour += 1.0 ) {
	    test_sun_one ( &now, hour );
	}
}

void
test_sun2 ( void )
{
	struct time now;
	double lt_rise, lt_set;
	int rises, sets, above;

	// set_time ( &now, 2017, 10, 5 );
	set_time ( &now, 2017, 10, 4 );

	sun_events ( &now, &lt_rise, &lt_set, &rises, &sets, &above );
	printf ( "Sun rise and set: %.3f %.3f\n", lt_rise, lt_set );
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
s_dms ( char *buf, double deg )
{
	int d, m;

	d = floor ( deg );
	deg -= d;
	deg *= 60.0;
	m = floor ( deg );
	deg -= m;
	deg *= 60.0;

	sprintf ( buf, "%d:%d:%.3f", d, m, deg );
	return buf;

	/*
	printf ( "D = %d\n", d );
	printf ( "M = %d\n", m );
	printf ( "S = %.3f\n", deg );
	*/
}

void
init_site ( struct site *sp )
{
	sp->lat_deg = LATITUDE;		/* degrees N */
	sp->long_hours = LONGITUDE;	/* hours west */
	sp->long_deg = sp->long_hours * -360.0 / 24.0;	/* degrees E */

	sp->latitude = sp->lat_deg * DEGRAD;
	sp->longitude = sp->long_deg * DEGRAD;
	sp->elevation = ELEVATION;      /* site elevation in kilometers (never used) */
	sp->tz = TIMEZONE;		/* hours from greenwich */

	sp->sin_lat = sin(sp->latitude);
	sp->cos_lat = cos(sp->latitude);
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
        // lst = gmst + e;

        gmst -= 2.9e-4 * sin (omega);

        gmst = fmod (gmst, 24.0);
        if (gmst < 0.0)
            gmst += 24.0;

	return gmst;
}

double
mmt_lst ( double jd )
{
	double lst;

	lst =  mmt_st ( jd ) - site_info.long_hours;
        if (lst < 0.0)
            lst += 24.0;
        if (lst > 24.0)
            lst -= 24.0;
	return lst;
}
#endif

/* Calculate the MJD for a given year, month, and day
 *  For years after 1582, we are using the Gregorian calendar,
 *  so for my purposes we could eliminate the test below.
 */
double 
calc_mjd ( struct time *tp )
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

/* Obliquity of the ecliptic
 * t = julian centuries from 2000 January  1 (12h)
 * (Julian century is always 36525 days)
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

/* Transform from ecliptic
 *  to equatorial (ra, dec) coordinates.
 * Starts with lat and long in radians.
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
 *
 * t = time in julian centuries since J2000
 */
void
MiniMoon ( double t, double *ra, double *dec )
{
	double l0, lm, ls, d, f;
	double dl, s, h, n;
	double l_moon, b_moon;

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
	l_moon = TWOPI * Frac ( l0 + dl / ARCSEC_360 );
	b_moon = ( 18520.0 * sin(s) + n ) * ASEC_TO_RAD;

	ll_to_rd ( b_moon, l_moon, ra, dec );
}

void
checksun ( void )
{
	double jd;
	double n, l, g, g2, ll;

	jd = 2458031.47420;
	n = jd - 2451545.0;
	l = 280.460 + 0.9856474 * n;	/* mean long - degrees */
	g = 357.528 + 0.9856003 * n;	/* anomaly - degrees */

	// printf ( "Check, l, g = %.2f %.2f\n", l, g );
	l = fmod ( l, 360.0 );
	g = fmod ( g, 360.0 );
	printf ( "Check, l, g = %.2f %.2f\n", l, g );

	g2 = g + g;
	g2 = fmod ( g2, 360.0 );

	ll = l + 1.915 * sin(g*DEGRAD) + 0.020 * sin(g2*DEGRAD);
	printf ( "Check, solar long (deg) = %.3f\n", ll );
}

/* From APC, page 39
 * This is a reduced version to calculate solar positions
 *  suitable for sunrise and sunset times.
 * Note that the sun always is in the ecliptic plane, so the
 *  solar ecliptic latitude is always zero.
 * (it never exceeds 0.00033 degrees).
 *
 * t = time in julian centuries since J2000
 */
void
MiniSun ( double t, double *ra, double *dec )
{
	double m;
	double sun_long;

	m = TWOPI * Frac ( 0.993133 + 99.997361 * t );
	sun_long = TWOPI * Frac ( 0.7859453 + m / TWOPI +
	    (6893.0 * sin(m) + 72.0 * sin(2.0*m) + 6191.2 * t ) / ARCSEC_360 );

	// printf ( "MiniSun, solar long = %.3f\n", sun_long * RADDEG );

	// checksun ();

	ll_to_rd ( 0.0, sun_long, ra, dec );
}

void
set_time ( struct time *tp, int y, int m, int d )
{
	tp->year = y;
	tp->month = m;
	tp->day = d;
	tp->mjd0 = calc_mjd ( tp );
}

/* THE END */
