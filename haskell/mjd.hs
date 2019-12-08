#!/bin/runghc

fint :: Double -> Double
fint x = fromIntegral $ floor x

-- DZERO - Julian Date at Jan 0 relative to
--    Epoch 1900, Jan 0.5
--    almanac, page B12 in 2006 and 2007 and 2008
--    almanac, page B13 in 2009, ... 2013, 2014
-- (note, 2003 Almanac page K5 says that
-- 1900 January 0 at 0h is JD = 2415019.5
-- we add the 0.5 to reckon JD from noon as usual.)

-- (2007) jd = 2454100.5
-- (2009) jd = 2454831.5
-- (2010) jd = 2455196.5
-- (2011) jd = 2455561.5
-- (2012) jd = 2455926.5
-- (2013) jd = 2456292.5
-- (2014) jd = 2456657.5
-- (2015) jd = 2457022.5

i_jd = 2457022.5

-- we called this DZERO in the old days.
dz jd = jd - j1900
    where j1900 = 2415020.0

jd = dz i_jd

-- Here is a routine to calculate mjd
-- given time in seconds from 1970
--   just like we would get from a call to unix time()
-- (from MMT mount code)
-- MJD begins at Greenwich midnight
-- and is zero for November 17, 1858
--  note that mjd = jd - 2400000.5

get_mjd sec_1970 = mjd_1970 + sec_1970 / sec_per_day
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

-- There are Julian and Gregorian calendars.
-- The Julian calendar is ancient (pre 1582)
-- The Gregorian is what we use now.
-- The Julian date counts from noon Jan 1, 4713 BC
-- for my use the dates pre Gregorian reform
--  are not of interest.
-- Modified Julian dates count from midnight (0 hours)
-- on November 17, 1858, so the mjd rolls over at
-- midnight, not noon as for the jd.

calc_mjd y m d
    | m <= 2 = calc_mjd (y-1) (m+12) d
    | otherwise = 365*y - 679004 + b + fint(30.6001*(m+1)) + d
    where b = (y/400) - (y/100) + (y/4)

jd2mjd jd = jd - 2400000.5

jdx = get_jd 2014 1 1
jdc = get_jd 2014 1 1
mjdx = jd2mjd jdx

x = dz jdx

-- ans = show jd
-- ans = show jdx ++ " " ++ show mjdx
ans = show jdx ++ " " ++ show x 

main = do
    putStrLn ans
