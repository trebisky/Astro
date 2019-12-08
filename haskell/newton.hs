#!/usr/bin/runghc
-- MMT almanac program
-- newtons method for sqrt
--
-- Tom Trebisky 10-16-2014

xabs :: Double -> Double
xabs x
    | x < 0.0 = -x
    | otherwise  = x

-- Here is a working general Newtons method
next f fd x = x - ((f x) / (fd x))

newton f fd startx = iterate iter startx
    where iter = next f fd

-- test it on the sin function
-- ans = take 10 $ newton sin cos 3.0

-- test it on the square root.
ans = take 10 $ newton (\x -> x*x-4.0) (\x -> 2.0*x) 1.5

-- iter = next sin cos
-- ans2 = iter ( iter 3.0 )
-- ans1 = take 10 $ iterate iter 2.0

main = do
    putStrLn $ show ans
