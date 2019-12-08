#!/usr/bin/runghc
-- MMT almanac program
-- Tom Trebisky 10-16-2014
--
-- There was once an ancient and extremely contorted fortran
-- program of uncertain origin.  It would occasionally do odd
-- things.  In October, 2014 it was rewritten in Haskell and
-- this is the result.  An effort has been made to clarify many
-- things that were mysterious in the old program as their
-- intent has been deciphered.
--
-- There are places in this code where I could do things more
-- efficiently, but on my current machine (an Intel i7 3.5Ghz)
-- the calculations for a year run in about half a second using
-- the Haskell interpreter.  The Haskell could be compiled and
-- run perhaps 8 times faster, but who cares, why bother.
-- Clarity and correctness are far more important.

-- import Data.List
import Text.Printf
import Debug.Trace

import Params

year = fromIntegral i_year

-- here are the values for 2015
--      DATA TZONE / 7.0 /
--      DATA ALON / 7.392303704 /
--      DATA PHI / 0.5530735082660357 /

--      DATA STMID1 / 6.688730472222222 /
--      DATA DZERO / 42002.5 /
--      DATA TAUSUN / 4.291666666666667 /

--      DATA TSR / 2.291703518522222 /
--      DATA TSS / 1.7133312037000001 /
--      DATA TMR / 1.5992779166638889 /
--      DATA TMS / 2.1791943055583336 /

--      DATA ETASUN / 0.9856065866666668 /
--      DATA OMGSUN / 283.58802268518514 /
--      DATA ECCSUN / 0.016702111000000002 /
--      DATA EPSUN / 23.4372623 /

--      DATA RIMOON / 5.1566898 /

--      DATA ECCM / 0.054900489 /
--      DATA ANEWM / 19.25972222222222 /

-- here are the values for 2015
--      DATA TZONE / 7.0 /
--      DATA ALON / 7.392303704 /
--      DATA PHI / 0.5530735082660357 /

--      DATA STMID1 / 6.704738333333333 /
--      DATA DZERO / 41637.5 /
--      DATA TAUSUN / 4.5 /

--      DATA TSR / 2.291009074077778 /
--      DATA TSS / 1.7133312037000001 /
--      DATA TMR / 1.2775800771638888 /
--      DATA TMS / 2.7779846913555555 /

--      DATA ETASUN / 0.9856183066666666 /
--      DATA OMGSUN / 284.0503277777778 /
--      DATA ECCSUN / 0.016702531 /
--      DATA EPSUN / 23.4373923 /

--      DATA RIMOON / 5.1566898 /

--      DATA ECCM / 0.054900489 /
--      DATA ANEWM / 0.1763888888888887 /

-- =================================
-- year = 2014
-- 
-- i_stmid1 = 6.704738333333333
-- 
-- i_mjd0 = 41637.5
-- i_tausun = 4.5
-- 
-- i_tsr = 2.291009074077778
-- i_tss = 1.7133312037000001
-- i_tmr = 1.2775800771638888
-- i_tms = 2.7779846913555555
-- 
-- i_etasun = 0.9856183066666666
-- i_omgsun = 284.0503277777778
-- i_epsun = 23.4373923
-- i_rimoon = 5.1566898
-- 
-- eccsun = 0.016702531
-- eccmoon = 0.054900489
-- 
-- anewm = 0.1763888888888887

-- =================================
-- year = 2015
-- 
-- i_stmid1 = 6.688730472222222
-- 
-- -- looks to me like we have been a day ahead here
-- -- for quite some time, but we are sticking with it
-- -- for now till we get some things sorted out.
-- i_mjd0 = 42002.5
-- -- i_mjd0 = 42001.5
-- i_tausun = 4.291666666666667
-- 
-- i_tsr = 2.291703518522222
-- i_tss = 1.7133312037000001
-- i_tmr = 1.5992779166638889
-- i_tms = 2.1791943055583336
-- 
-- i_etasun = 0.9856065866666668
-- i_omgsun = 283.58802268518514
-- i_epsun = 23.4372623
-- i_rimoon = 5.1566898
-- 
-- eccsun = 0.016702111000000002
-- eccmoon = 0.054900489
-- 
-- anewm = 19.25972222222222

-- =================================
-- These things are site specific and never change

-- alon = 7.392303704
-- my_lat = 0.5530735082660357

tzone = 7.0
i_my_lat = 31.688777784
i_my_long = 110.88455556

-- ---------------------------------------
-- ---------------------------------------
-- Values above here are "injected" into the program
-- after being extracted from the Astronomical Almanac.
-- Another program exists that processes raw information from
-- the almanac, then generates the above values.
-- Here are some notes on what they are:
--
-- Tausun = 3.83  This is the date of perihelion,
--    i.e the day when the earth is closest to the sun.
--    these days January 4 more or less.
-- Omgsun = 283.06 degrees This is the longitude of the sun
--                  at perihelion
-- Eccsun = 0.0167 The eccentricity of the suns orbit
--                 at mid year.
-- Etasun = 0.9856 degrees The average daily motion of the sun.
--    This should be about 360/365 = 0.986, and it is!
-- Epsun = 23.4383 degrees The obliquity of the ecliptic
--                 at mid year.
-- Dzero = 39080.5 MJD at January 0
-- Eccm = 0.0549 eccentricity of the lunar orbit.
-- Rimoon = 5.14 degrees - mean inclination of lunar orbit to the eccliptic.
-- Anewm = 17.8757 = Date of first new moon at Greenwich.
--
-- Phi = 0.553 MMT latitude (radians)
-- Alon = 7.39 MMT longitude (hours west)
-- Tzone = 7.0 MMT timezone (hours)
--
-- Stmid1 = 6.68 sidereal time at midnight January 1 (hours)
--
-- Tss = 1.71 time of sunset, January 1 (days)
-- Tsr = 2.29 time of sunrise, January 2 (days)
-- Tmr = 1.62 time of moonrise, January 1 (days)
-- Tms = 2.25 time of moonset, January 2 (days)

-- ---------------------------------------
-- ---------------------------------------
-- Note that the original code used some truncated values
-- for the following, such as:
-- pi = 3.14159
-- rad2hours = 3.81972
-- Because of this, we see some 1 second errors in the
-- generated values when compared with the original Fortran code.

twopi = 2.0 * pi
rad2deg = 180.0 / pi
deg2rad = pi / 180.0
rad2hours = 24.0 / twopi
deg2hours = 24.0 / 360.0
arcsec = 3600.0 * rad2deg

fabs :: Double -> Double
fabs x
    | x < 0.0 = -x
    | otherwise  = x

-- The following two functions are to manhandle
-- Haskell's type system

fint :: Double -> Double
fint x = fromIntegral $ floor x

fround :: Double -> Double
fround x = fromIntegral $ round x

-- This is how much (in hours) the sidereal time advances each day.
-- Sidereal time moves at 15.04106699 arc-seconds / solar second
-- Which is 1.002737799 sidereal second / solar second.
-- In 24 solar hours, sidereal time advances 24.0657072 sidereal hours)
-- Hence the 0.0657098 factor (Tom gets 0.0657072, but that is OK)
-- With an advance each day of 0.0657 hours, in 365 days this is 23.98 hours.

stperday::Double
-- stperday = 0.0657098
stperday = 0.0657072

-- conversion factor from sidereal to solar time.
-- sid2solar = 0.99727
sid2solar :: Double
sid2solar = 24.0 / (24.0 + stperday )

-- The Synodic Month - the mean time in days between two
-- consecutive new moons (the actual time varies substantially)
synm::Double
synm = 29.530589

my_lat = i_my_lat * deg2rad
my_long_days = i_my_long / 360.0
my_long_hours = i_my_long * deg2hours

-- Generate a correction factor in units of days for our precise longitude
-- cor_long = (alon - tzone) / 24.0
cor_long = my_long_days - tzone / 24.0

-- Correct the values that need correcting.

stmid1 = i_stmid1 + tzone + 0.00274 * tzone - my_long_hours
tausun = i_tausun - my_long_days

-- Here is a routine to calculate mjd
-- given time in seconds from 1970
--  (like we would get from a call to unix time())
-- MJD begins at Greenwich midnight
-- and is zero for November 17, 1858
--  note that mjd = jd - 2400000.5

get_mjd t = mjd_1970 + t / sec_per_day
    where mjd_1970 = 40587.0
          sec_per_day = 24.0 * 60.0 * 60.0

-- From Meeus page 61
-- This is valid for the Gregorian Calendar
-- (i.e. any date after October 15, 1582)
get_jd y m d = p + q + d + b - 1524.5
    where a = fint ( y / 100.0 )
          b = 2 - a + fint ( a/4 )
	  p = fint ( 365.25 * (y + 4716))
	  q = fint ( 30.6001 * (m + 1))

jd2mjd jd = jd - 2400000.5

-- The original single precision floating point lost precision here.
-- The MJD value is something like 42002.5 and the correction
-- is something like 0.308xxx, only the 3 digits .308 were retained
-- when single precision was used.

mjd0 = i_mjd0 + my_long_days

tss = i_tss + cor_long
tsr = i_tsr + cor_long
tms = i_tms + cor_long
tmr = i_tmr + cor_long

omgsun = i_omgsun * deg2rad
etasun = i_etasun * deg2rad
epsun  = i_epsun * deg2rad
rimoon = i_rimoon * deg2rad

-- Zenith distances in radian units.
--      ZD(1)=1.58534   = 90.83329109326993
--      ZD(5)=1.58534
--      ZD(2)=1.67552   = 96.00022448975969
--      ZD(6)=1.67552
--      ZD(3)=1.78024   = 102.00023852036968
--      ZD(7)=1.78024
--      ZD(4)=1.88496   = 108.00025255097965
--      ZD(8)=1.88496
--      ZD(9)=1.56861   = 89.87473270201606
--      ZD(10)=1.56861

zd_sun = 90.83329109326993 * deg2rad
zd_sun6 = 96.00022448975969 * deg2rad
zd_sun12 = 102.00023852036968 * deg2rad
zd_sun18 = 108.00025255097965 * deg2rad
zd_moon = 89.87473270201606 * deg2rad

-- We have leap years because a tropical year
-- is actually 365.242190 days in length.
-- On a leap year, February has 29 days instead of 28.
-- Basically every year divisible by 4 is a leap year
-- with some exceptions.  2012 is a leap year.

isLeap::Int->Bool
isLeap y
    | mod y 400 == 0 = True
    | mod y 100 == 0 = False
    | mod y 4 == 0 = True
    | otherwise = False

-- days per month
dpm::Int->[Int]
dpm y
    | isLeap y = [ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
    | otherwise = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]

month_days = dpm year

month_names = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]

-- this array gives cummulative day numbers at the start of each month
-- note that I prepend a zero and will never use the last entry.
-- For a non leap year:
--     0,31,59,90,120,151,181,212,243,273,304,334,366

month_dayno = 0 : scanl1 (+) month_days

dayno::Int->Int->Double
dayno m d = fromIntegral dn
    where dn = month_dayno!!(m-1) + d

stmid m d
    | st > 24.0 = st - 24.0
    | otherwise = st
    where dn = dayno m d
          st = stmid1 + dn * stperday

--dpy y
--    | isLeap y = 366
--    | otherwise = 365

-- Output conversion routines
-- Note that rounding is tricky if we do it after
-- converting to HM, such as if we round 59.7 minutes,
-- we would get 60 minutes, but would need to carry
-- this to adding an hour (and setting minutes to 0).
-- We avoid all this by adding the half second or minute
-- to the whole double value before starting the breakdown.

half_minute = 0.5 / 60.0
half_second = 0.5 / (60.0 * 60.0)

hms::Double->String
hms arg = h ++ ":" ++ m ++ ":" ++ s
    where hh = arg + half_second
          hi = floor hh :: Int
          h = printf "%02d" hi
          mm = (hh-(fromIntegral hi)) * 60.0
	  mi = floor mm :: Int
          m = printf "%02d" mi
          ss = (mm-(fromIntegral mi)) * 60.0
          s = printf "%02d" (floor (ss) :: Int)

hm::Double->String
hm arg = h ++ ":" ++ m
    where hh = arg + half_minute
          hi = floor hh :: Int
          h = printf "%02d" hi
          mm = (hh-(fromIntegral hi)) * 60.0
          m = printf "%02d" (floor (mm) :: Int)

-- -----------------------------------------------------------
-- Solve keplers equation via newtons method
-- -----------------------------------------------------------

-- Here is a working general Newtons method
next f fd x = x - ((f x) / (fd x))

newt f fd startx n
    | n > 20 = (xn, n)
    | fabs ( startx - xn) < 0.0001 = (xn, n)
    | otherwise = newt f fd xn (n+1)
    where xn = next f fd startx

-- use this when you want the tuple
-- fst = result, snd = count of iterations
newtt f fd startx = newt f fd startx 0

-- use this when you just want the result
newton f fd startx = fst $ newt f fd startx 0

-- ans = newtt (\x -> x - eccsun * sin(x) - a) (\x -> 1.0 - eccsun * cos(x)) a

kepler a e = newton f fd start
--    where start = trace ("kepler " ++ show a) a
    where start = a
          f x = x - e *sin(x) - a
          fd x = 1.0 - e * cos(x)

-- Convert eccentric anomaly to true anomaly
ecc2tru e ecc = 2.0 * atan(u)
    where ratio = sqrt ((1.0 + ecc) / (1.0 - ecc) )
          u = ratio * tan ( e/2.0 )

-- -----------------------------------------------------------
-- Solving  keplers equation we begin with
-- the time since epoch (t-tau) in days
-- ( etasun is the average daily solar motion,
--   and has units of radians per day).
-- This starting point is the mean anomaly
-- solving Keplers equation gives us the eccentric anomaly,
-- which we then transform into the true anomaly,
-- which is what we really want.

trusun t = ecc2tru e eccsun
    where a = etasun * (t - tausun)
          e = kepler a eccsun

-- my own flavor of atan2 that returns in [0,twopi]
my_atan2 :: Double -> Double -> Double
my_atan2 y x
    | rv < 0.0 = rv + twopi
    | otherwise = rv
    where rv = atan2 y x

-- Calculate RA and Dec for the Sun given the true anomaly
-- Note than in both Fortran and Haskell, the order of arguments
--  for atan2 is Y then X

rdsun tru = (ra, dec)
    where trulng = omgsun + tru
          a1 = cos(epsun) * sin(trulng)
          a2 = cos(trulng)
          ra = my_atan2 a1 a2
          sdelt = sin(trulng) * sin(epsun)
          dec = asin ( sdelt)

-- constrain value to [0,twopi]
tp_constrain v
    | v > twopi = v - twopi
    | otherwise = v

-- Calculate the RA and Dec of the moon for a certain time
-- Tested OK
rdmoon:: Double -> (Double,Double)
rdmoon t = (ra, dec)
    where mjd = mjd0 + t
	  tj = ((fromIntegral year)-1900.0) / 100.0
	  a1 = 270.434164 + 13.1763965268*mjd - 0.001133*tj*tj +0.0000019*tj*tj*tj
	  a3 = a1 - 360.0 * fint (a1/360.0)
	  gmoonl = a3 * deg2rad
	  b1 = 334.329556 + 0.1114040803*mjd - 0.010325*tj*tj - 0.000012*tj*tj*tj
	  b3 = b1 - 360.0 * fint (b1/360.0)
	  perl = b3 * deg2rad

          e = kepler (gmoonl - perl) eccmoon
	  v = ecc2tru e eccmoon

	  appl = tp_constrain ( v + perl )
	  da = t + fround ( 365.25 * ((fromIntegral year)-1970.0))
	  c1 = -791.12 + 11.316506 * da
          evec = 0.02225 * sin(c1*deg2rad)
          d1 = -908.88 + 24.381498 * da
          var = 0.01149 * sin(d1*deg2rad)
          e1 = 356.8 + 0.985600 * da
          ae = -0.00323 * sin(e1*deg2rad)
          rlmoon = appl+evec+var+ae

          f1 = 259.183275 - 0.0529539222*mjd + 0.002078*tj*tj + 0.000002 * tj*tj*tj
          anl = (f1 - 360.0 * fint (f1/360.0) ) * deg2rad

          sinb = sin(rlmoon-anl) * sin(rimoon)
	  bmoon = asin(sinb)
	  sind = sinb * cos(epsun) + cos(bmoon)*sin(epsun)*sin(rlmoon)
	  dec = asin(sind)

	  cosa = cos(bmoon)*cos(rlmoon)/cos(dec)
	  sina = (-sin(bmoon)*sin(epsun) + cos(bmoon)*cos(epsun)*sin(rlmoon)) / cos(dec)
	  ra = atan2 sina cosa

-- ensure a value in hours is in the range [0-24]
hconstrain h
    | h > 24.0 = h - 24.0
    | h < 0.0 = hconstrain (h+24.0)
    | otherwise = h

-- convert from sidereal time in hours
-- to a time in solar hours
frmsd x
    | f < 0.0 = f + 24.0
    | otherwise = f
    where f = x * sid2solar

-- transform from zd and dec to hour angle
zd2ha dec zd = acos ( h1 / h2 )
    where h1 = cos(zd) - sin(my_lat) * sin(dec)
          h2 = cos(my_lat) * cos(dec)

-- Sunrise and sunset calculations are the same,
--  except:
--   sunrise subtracts rather than adds HA
--   we back up a sidereal day for sunset

st_set (ra, dec) zd stm = sidt
    where ha = zd2ha dec zd
	  sidt = (ra + ha) * rad2hours

st_rise (ra, dec) zd stm = sidt
    where ha = zd2ha dec zd
	  sidt = (ra - ha) * rad2hours

calc_sunset rd zd stm = time
    where sidt = st_set rd zd stm
	  stel = hconstrain $ sidt - (stm - stperday)
	  time = frmsd ( stel )

calc_sunrise rd zd stm = time
    where sidt = st_rise rd zd stm
	  stel = hconstrain $ (sidt - stm)
	  time = frmsd ( stel )

calc_moonrise rd zd stm = time
    where sidt = st_rise rd zd stm
	  stel = hconstrain (sidt - stm)
	  time = frmsd ( stel )

calc_moonset rd zd stm = time
    where sidt = st_set rd zd stm
	  stel = hconstrain $ (sidt - stm)
	  time = frmsd ( stel )

--calc_moonset (ra, dec) zd stm = trace ("frmsd (hours) = " ++ show time) time
--    where ha = zd2ha dec zd
--	  sidt = trace ("moonset, ha, ra = " ++ show ha ++ " " ++ show ra) $ (ra + ha) * rad2hours
--	  stel = trace ("moonset, sidt, stm = " ++ show sidt ++ " " ++ show stm) $ hconstrain (sidt - stm)
--	  time = frmsd ( stel )

-- We iterate on the moon rise/set times till we get a time that is
-- accurate to within one minute.
-- This is required because the moon moves almost an hour on the sky
-- each night (so from zenith to horizon is 15 minutes).
--
-- The sun on the other hand moves about 24/365 hours from "night to night",
-- which is about 3.9 minutes on the sky and from zenith to horizon would
-- be 3.9/4.0 = 0.98 minutes (just right!)  No need to iterate when we
-- are printing times to the nearest minute.

minute_as_day :: Double
minute_as_day = 1.0 / 60.0 / 24.0

iter_moonrise dn t zd stm count
    | count >= 10 = mr
    | fabs (tnew-t) < minute_as_day = mr
    | otherwise = iter_moonrise dn tnew zd stm (count+1)
--    where rd = trace ("--moonrise, t = " ++ show t) $ rdmoon t
    where rd = rdmoon t
	  mr = calc_moonrise rd zd stm
	  tnew = dn + mr / 24.0

iter_moonset dn t zd stm count
    | count >= 10 = mr
    | fabs (tnew-t) < minute_as_day = mr
    | otherwise = iter_moonset dn tnew zd stm (count+1)
--    where rd = trace ("--moonset, t = " ++ show t) $ rdmoon t
    where rd = rdmoon t
	  mr = calc_moonset rd zd stm
	  tnew = dn + mr / 24.0

-- Sidereal Time at midnight.

-- specific calculations for a month and day

moon_per_day :: Double
moon_per_day = 1.0 + 1.0 / synm

-- moonrise m d zd = trace ("*** Moonrise: " ++ show mr ++ " " ++ (hm mr)) $ hm mr
moonrise m d zd = hm mr
    where dn = dayno m (d+1)
          t = tmr + (dn-1.0) * moon_per_day
	  mr = iter_moonrise dn t zd (stmid m d) 0

-- moonset m d zd = trace ("*** Moonset: " ++ show ms ++ " " ++ (hm ms)) $ hm ms
moonset m d zd = hm ms
    where dn = dayno m (d+1)
          t = tms + (dn-2.0) * moon_per_day
	  ms = iter_moonset dn t zd (stmid m d) 0

moonage :: Int -> Int -> String
moonage m d
    | age >= 14.0 = printf "%.1f" (age - 28.0)
    | otherwise = printf "%.1f" age
    where dn = dayno m d
          age1 = (dn - anewm + synm) / synm
	  age2 = ( age1 - fint(age1) ) * synm
	  age = age2 * 28.0 / synm

sunrise m d zd = hm sr
    where dn = dayno m d
          t = tsr + (dn-1.0)
          rd = rdsun $ trusun t
	  sr = calc_sunrise rd zd (stmid m d)

sunset m d zd = hm ss
    where dn = dayno m d
          t = tss + (dn-1.0)
          rd = rdsun $ trusun t
	  ss = calc_sunset rd zd (stmid m d)

-- These are modified versions of the sunrise and sunset functions above
-- We repeat the bulk of the sunrise/sunset function, but who cares
-- the whole program runs in less than a second.
ra3w m d
    | ra > 24.0 = ra - 24.0
    | otherwise = ra
    where dn = dayno m d
	  zd = zd_sun18
          t = tss + (dn-1.0)
          rd = rdsun $ trusun t
	  st = st_set rd zd (stmid m d)
	  ra = st - 3.0

ra3e m d
    | ra < 0.0 = ra + 24.0
    | otherwise = ra
    where dn = dayno m d
	  zd = zd_sun18
          t = tsr + (dn-1.0)
          rd = rdsun $ trusun t
	  st = st_rise rd zd (stmid m d)
	  ra = st + 3.0

todate m d = month_names!!(m-1) ++ " " ++ printf "%2d" d

-- all calculations for a day
-- "m" is month number (1-12)
-- "d" is day number (1-N)

do_day m d = date ++ " " ++
	dayn ++ " " ++
	ss   ++ "   " ++
	ss6  ++ "   " ++
	ss12 ++ "   " ++
	ss18 ++ "   " ++
	r3w  ++ "   " ++
	stm  ++ "   " ++
	r3e  ++ "   " ++
	sr18 ++ "    " ++
	sr12 ++ "    " ++
	sr6  ++ "    " ++
	sr   ++ "    " ++
	mr   ++ "    " ++
	ms   ++ "    " ++
	ma
    where date = todate m d
	  dn = dayno m d
	  dayn = printf "%5.0f  " dn

	  sr = sunrise m d zd_sun
	  sr6 = sunrise m d zd_sun6
	  sr12 = sunrise m d zd_sun12
	  sr18 = sunrise m d zd_sun18

	  r3w = hm $ ra3w m d
          stm = hms $ stmid m d
	  r3e = hm $ ra3e m d

	  ss = sunset m d zd_sun
	  ss6 = sunset m d zd_sun6
	  ss12 = sunset m d zd_sun12
	  ss18 = sunset m d zd_sun18

	  mr = moonrise m d zd_moon
	  ms = moonset m d zd_moon
	  ma = moonage m d

-- all calculations for a month

label = "Date    Dayno  Sunset  SS6     SS12    SS18    RA 3W   ST(mid)    RA 3E   SR18     SR12     SR6     Sunrise  Moonrise Moonset  Moon Age"

do_month m = header ++ days
    where header = [" ", label ]
          days = map (do_day m) [1..month_days!!(m-1)]

do_debug m = header ++ days
    where header = [" ", label ]
          days = map (do_day m) [0..2]

-- all calculations for a year

-- results = map (do_day 1) [1]
-- results = do_debug 1
-- results = do_month 1
results = concat $ map do_month [1..2]
-- results = concat $ map do_month [1..12]

-- ans = unlines $ map show month_days
-- ans = unlines $ map show month_dayno
-- ans = unlines $ map date days

ans = unlines results

jdx = get_jd 2015 1 1
mjdx = jd2mjd jdx

-- ans = show jdx ++ " " ++ show mjdx

main = do
    putStr ans
