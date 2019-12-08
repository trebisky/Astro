#!/usr/bin/runghc
-- MMT almanac program
--
-- Tom Trebisky 10-16-2014

-- import Data.List

stmid1 = 6.688730472222222

-- Output conversion routines
hms::Double->String
hms hh = show h ++ ":" ++ show m ++ ":" ++ show s
    where h = floor hh
          mm = (hh-(fromIntegral h)) * 60.0
	  m = floor mm
          ss = (mm-(fromIntegral m)) * 60.0
	  s = floor (ss + 0.5)

hm::Double->String
hm hh = show h ++ ":" ++ show m
    where h = floor hh
          mm = (hh-(fromIntegral h)) * 60.0
	  m = floor mm

ans = hm stmid1 ++ "\n"

main = do
    putStr ans
