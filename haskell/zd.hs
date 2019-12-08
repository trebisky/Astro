#!/usr/bin/runghc
-- MMT almanac program
--
-- Tom Trebisky 10-16-2014

-- import Data.List
import Text.Printf

rad2deg = 180.0 / pi
deg2rad = pi / 180.0


-- Zenith distances in radian units.


zd = [ 1.58534, 1.67552, 1.78024, 1.88496, 1.56861 ]
ans = map (* rad2deg) zd

main = do
    putStrLn $ show ans
