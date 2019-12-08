#!/usr/bin/runghc
-- MMT almanac program
--
-- Tom Trebisky 10-16-2014

year = 2014

d_1970 :: Integer -> Double
d_1970 y = fromIntegral $ round (365.25 * ((fromIntegral y) - 1970.0))

xx = d_1970 year
ans = 2.5 + xx

main = do
    putStrLn $ show ans
