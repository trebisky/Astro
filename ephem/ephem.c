#include <stdio.h>
#include <math.h>

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
#define SECPD	86400.0		/* seconds per day */

/* Multiply an angle in radians by this to get arc-seconds */
/* Or divide an angle in arc-seconds by this to get radians */
#define ARCSEC	3600.0 * RADDEG

struct time {
	int year;
	int month;
	int day;
	double mjd;
};

struct site {
	double lat_deg;
	double long_deg;
	double long_hours;
	double latitude;
	double longitude;
	double elevation;
	double tz;
};

struct site site_info;

void init_site ( struct site * );
void set_time ( struct time *, int, int, int );
double gmst ( double );
char * s_dms ( double );

void test_jd ( void );

int
main ( int argc, char **argv )
{
	struct time now;
	double st;

	init_site ( &site_info );

	test_jd ();

#ifdef notdef
	set_time ( &now, 2018, 1, 1 );
	st = gmst ( now.mjd ) * RADHR;
	printf ( "STmid = %.7f\n", st );
	s_dms ( st );
#endif

	printf ( "Done\n" );
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

	set_time ( &now, year, 1, 0 );
	printf ( "jd (%d) = %.2f (%.2f)\n", year, now.mjd + MJD_OFFSET, aa_jd );
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

#ifdef notdef

# ----> STMID1 - Sideral Time at Midnight, Jan 1.
#       almanac page B12 in 2006
#       almanac page B13 in 2009, ... 2016
#       G. Sideral Time for Jan 1, 0h UT
# (2005) stmid = hms( 6.0, 42.0, 58.4748 )
# (2006) stmid = hms( 6.0, 42.0, 1.5159 )
# (2007) stmid = hms( 6.0, 41.0, 4.5504 )
# (2008) stmid = hms( 6.0, 40.0, 7.5881 )
# (2009) stmid = hms( 6.0, 43.0, 7.1394 )
# (2010) stmid = hms( 6.0, 42.0, 10.0357 )
# (2011) stmid = hms( 6.0, 41.0, 12.8088 )
# (2012) stmid = hms( 6.0, 40.0, 15.4861 )
# (2013) stmid = hms( 6.0, 43.0, 14.6094 )
# (2014) stmid = hms( 6.0, 42.0, 17.058 )
# (2015) stmid = hms( 6.0, 41.0, 19.4297 )
# (2016) stmid = hms( 6.0, 40.0, 21.7893 )
# (2017) stmid = hms( 6.0, 43.0, 20.7109 )
# (2018) stmid = hms( 6.0, 42.0, 23.1082 )
stmid = hms( 6.0, 42.0, 23.1082 )

#endif

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

char *
s_dms ( double deg )
{
	int d, m;

	d = floor ( deg );
	deg -= d;
	deg *= 60.0;
	m = floor ( deg );
	deg -= m;
	deg *= 60.0;
	printf ( "D = %d\n", d );
	printf ( "M = %d\n", m );
	printf ( "S = %.3f\n", deg );
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
}

/* Sidereal time at Greenwich */
double
gmst ( double mjd )
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

double
mmt_lst ( double sec )
// void synch_lst ( void )
{
        double t0, t, gmst, omega, e, last, lst;
        double utc, ut1;
        double jd, jd0;
        double lst_sec;
        // double sec;

        /* utc since 1970 in seconds */
        // sec = time_info.ut.tv_sec;
        // sec += time_info.ut.tv_nsec / 1.e9;

        /* Calculate utc for the current day in hours */
        utc = fmod ( sec, 86400.0 ) / 3600.0;
        // ut1 = utc + time_info.dut / 3600.0;
        ut1 = utc;

        /* JD right now */
        jd = sec / SECPERDAY + JD_1970;

        /* JD today at 0h UTC */
        jd0 = jd - ut1 / 24.0;

        /* 2451545 is JD for Jan 1.5, 2000 */
        t0 = (jd0 - 2451545.0) / 36525.0;
        t  = (jd  - 2451545.0) / 36525.0;

        gmst = 6.69737456 + 2400.051336 * t0 + 2.58622e-5 *
               t0 * t0 + 1.002737909 * ut1;

        omega = 2.182438586 - 33.75704592 * t + 36.145769e-6 * t * t;

        e = -2.9e-4 * sin (omega);

	/* The following in units of hours */
        // last = gmst - site_info.geodetic_long_h + e;
        last = gmst - site_info.long_hours + e;

        lst = fmod (last, 24.0);
        if (lst < 0.0)
            lst = lst + 24.0;

        lst_sec = lst * 3600.0;
        // time_info.lst.tv_sec = lst_sec;
        // time_info.lst.tv_nsec = (lst_sec - (double) time_info.lst.tv_sec) * 1.e9;
	return lst_sec;
}
#endif


/* Calculate the MJD for a given year, month, and day
 *  For years after 1582, we are using the Gregorian calendar,
 *  so for my purposes we could eliminate the test below.
 */
double 
Mjd ( struct time *tp )
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

/* From APC, page 38
 * This is a reduced version to calculate lunar position
 *  suitable for moonrise and moonset times
 */
void
MiniMoon ( double t, double *ra, double *dec )
{
	static double eps = 23.43929111 * DEGRAD;
	double l0, lm, ls, d, f;
	double dl, s, h, n;
	double l_moon, b_moon;
	double cb, x, y, z, v, w, rho;
	double coseps = 0.91748;	/* XXX */
	double sineps = 0.39778;	/* XXX */

	l0 = Frac ( 0.606433 + 1336.855225 * t);		/* mean longitude */
	lm = TWOPI * Frac ( 0.374897 + 1325.552410 * t);	/* moon - mean anomaly */
	ls = TWOPI * Frac ( 0.993133 + 99.997361 * t);	/* sun - mean anomaly */
	d = TWOPI * Frac ( 0.827361 + 1236.853086 * t);	/* diff - long, moon-sun */
	f = TWOPI * Frac ( 0.259086 + 1342.227825 * t);	/* dist from ascending node */

	// Perturbations in long and lat
	dl = 22640 * sin(lm) - 4586 * sin(lm-2*d) + 2370 * sin(2*d) + 769 * sin(2*lm)
	    -668*sin(ls) - 412 * sin(2*f) - 212 * sin(2*lm-2*d) - 206 * sin(lm+ls-2*d)
	    +192 * sin(lm+2*d) - 165 * sin(ls-2*d) - 125 * sin(d) - 110 * sin(lm+ls)
	    +148 * sin(lm-ls) - 55 * sin(2*f-2*d);

	s = f + (dl + 412*sin(2*f) + 541 * sin(ls)) / ARCSEC;
	h = f - 2*d;
	n = -526 * sin(h) + 44 * sin(lm+h) - 31 * sin(-lm+h) -23 * sin(ls+h)
	    + 11 * sin(-ls+h) - 25 * sin(-2*lm+f) + 21 * sin(-lm + f);

	// Ecliptic long and lat
	l_moon = TWOPI * Frac ( l0 + dl/1296.0e3 );	/* radians */
	b_moon = ( 18520.0 * sin(s) + n ) / ARCSEC;

	// Transform to equatorial coordinates
#ifdef MOON_MATRIX
	e_moon = R_x(-eps) * Vec3D(Polar(l_moon,b_moon));
	*ra = e_moon[phi];
	*dec = e_moon[theta];
#else
	cb = cos ( b_moon );
	x = cb * cos ( l_moon );
	v = cb * sin ( l_moon );
	w = sin ( b_moon );
	y = coseps * v - sineps * w;
	z = sineps * v + coseps * w;
	rho = sqrt ( 1.0 - z*z );
	*dec = RADDEG * atan ( z / rho );
	*ra = RADHR * atan ( y / (x + rho) );
	if ( *ra < 0.0 ) *ra += 24.0;
#endif
}

void
set_time ( struct time *tp, int y, int m, int d )
{
	tp->year = y;
	tp->month = m;
	tp->day = d;
	tp->mjd = Mjd ( tp );
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
	// return 23.5 * DEGRAD;
}

/* THE END */
