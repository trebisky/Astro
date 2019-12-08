#!/usr/bin/runghc
-- MMT almanac program
-- newtons method for sqrt
--
-- Tom Trebisky 10-16-2014

xabs :: Double -> Double
xabs x
    | x < 0.0 = -x
    | otherwise  = x

-- test for convergence
-- we can test if f(x) is close enough to zero
-- or we can test if x is changing by some delta

isok :: (Double -> Double) -> Double -> Bool
isok f x = xabs(f x) < 0.000001

isok2 x xn = xabs(xn - x) < 0.0000001

-- Here is a working general Newtons method
next f fd x = x - ((f x) / (fd x))

--    | isok f xn = (xn, n)

newt f fd startx n
    | n > 10 = (xn, n)
    | isok2 startx xn = (xn, n)
    | otherwise = newt f fd xn (n+1)
    where xn = next f fd startx

-- use this when you want the tuple
-- fst = result, snd = count of iterations
newtt f fd startx = newt f fd startx 0

-- use this when you just want the result
newton f fd startx = fst $ newt f fd startx 0

-- test it on the square root.
-- ans = newton (\x -> x*x-100.0) (\x -> 2.0*x) 100.0

-- test it on kepler
ecc = 0.016702
a = -0.03867

f x = x - ecc *sin(x) - a
fd x = 1.0 - ecc * cos(x)

-- ans = newtt (\x -> x - ecc * sin(x) + a) (\x -> 1.0 - ecc * cos(x)) a
ans = newtt f fd a

main = do
    putStrLn $ show ans

-- THE END
