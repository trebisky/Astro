#!/bin/runghc

deg2rad = pi / 180.0

fint :: Double -> Double
fint x = fromIntegral $ floor x

frac :: Double -> Double
frac x = (x - fint x)
-- frac x = (x - fromIntegral $ floor x)

dot :: [Double] -> [Double] -> Double
dot a b = sum $ zipWith (*) a b

-- matrix times vector
mv m v = map (dot v) m

v1 = [ 1.0, 2.0, 3.0 ]
v2 = v1

vx = [ 1.0, 0.0, 0.0 ]
vy = [ 0.0, 1.0, 0.0 ]
vz = [ 0.0, 0.0, 1.0 ]

mm = [ [1.0, 0.0, 0.0],
	[0.0, 1.0, 0.0],
	[0.0, 0.0, 1.0]]

m2 = [ [2.0, 0.0, 0.0],
	[0.0, 2.0, 0.0],
	[0.0, 0.0, 2.0]]

rx ang = [ [1.0, 0.0, 0.0],
	  [0.0, c, s],
	  [0.0, -s, c]]
	  where c = cos angr
	        s = sin angr
		angr = ang * deg2rad

ry ang = [ [c, 0.0, -s],
	  [0.0, 1.0, 0.0],
	  [s, 0.0, c]]
	  where c = cos angr
	        s = sin angr
		angr = ang * deg2rad

rz ang = [ [c, s, 0.0],
	  [-s, c, 0.0],
	  [0.0, 0.0, 1.0]]
	  where c = cos angr
	        s = sin angr
		angr = ang * deg2rad

-- nv = mv mm v1

-- ans = show $ frac 1.234
-- ans = show $ dot v1 v2
-- ans = show nv

-- ans = show $ mv (rx 90.0) vz
-- ans = show $ mv (ry 90.0) vz
ans = show $ mv (rz 90.0) vx

main = do
    putStrLn ans
