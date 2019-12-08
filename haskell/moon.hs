#!/bin/runghc

-- moon.hs
-- calculate moonrise and moonset
-- code from "Astronomy on the Personal Computer"
--  by Montenbruck and Pfleger

fint :: Double -> Double
fint x = fromIntegral $ floor x

frac :: Double -> Double
frac x = (x - fint x)

twopi = 2.0 * pi
rad2deg = 180.0 / pi
deg2rad = pi / 180.0
rad2hours = 24.0 / twopi
deg2hours = 24.0 / 360.0
arcsec = 3600.0 * rad2deg

-- t is time in julian centuries since J2000

-- XXX the following needs careful checking for typos
-- then I need to implement the matrix routines in Haskell

miniMoon t = (ra, dec)
    where eps = 23.43929111 * deg2rad
          l0 = frac ( 0.606433 + 1336.855225 * t)
          l = twopi * frac ( 0.374897 + 1325.55241 * t)
	  ls = twopi * frac ( 0.993133 + 99.997361 * t)
	  d = twopi * frac ( 0.827361 + 1236.853086 * t)
	  f = twopi * frac ( 0.259086 + 1342.227825 * t)
	  dl = 22640.0*sin(l) - 4586.0 * sin(l-2.0*d) + 2370.0 * sin(2.0*d) +
	    769.0 * sin(2.0*l) - 668.0 * sin(ls) -412.0 * sin(2.0*f) -
	    212.0 * sin(2.0*l-2.0*d) - 206.0 * sin(l+ls-2.0*d) +
	    192.0 * sin(l+2.0*d) - 165.0 * sin(ls-2.0*d) - 125.0 * sin(d) -
	    110.0 * sin(l+ls) + 148.0 * sin(l-ls) - 55.0*sin(2.0*f-2.0*d)
	  s = f + (dl + 412.0 * sin(2.0*f) + 541.0 * sin(ls)) / arcsec
	  h = f - 2.0*d
	  n = -526.0 * sin(h) + 44.0 * sin(l+h) - 31.0 * sin(-l+h) - 23.0 * sin(ls+h) +
	    11.0 * sin(-ls+h) - 25.0 * sin(-2.0*l+f) + 21.0 * sin(-l+f)
	  -- ecliptic longitude and latitude
	  l_moon = twopi * frac ( l0 + dl / 1296.0e3 )
	  b_moon = ( 18520.0 * sin(s) + n ) / arcsec
	  -- equatorial coordinates
--	  e_moon = R_x(-eps) * Vec3D(Polar(l_moon,b_moon))
--	  ra = e_moon(phi)
--	  dec = e_moon(theta)
	  ra = 0.0
	  dec = 0.0

ans = show $ frac 1.234

main = do
    putStrLn ans
