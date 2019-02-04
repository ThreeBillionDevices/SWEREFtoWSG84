"""
Converted from Javascript to Python by Jacob Möller

Original Author: Arnold Andreasson, info@mellifica.se
Original File: http://latlong.mellifica.se/geodesi/gausskruger.js
For more see: http://latlong.mellifica.se/
License: MIT License as follows:

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

=============================================================================
Python-implementation of "Gauss Conformal Projection
(Transverse Mercator), Krügers Formulas".
- Parameters for SWEREF99 to SWEREF99 lat-long
  coordinates (SWEREF99 is used in Swedish maps).
Source: http://www.lantmateriet.se/geodesi/
"""

import math

# To find correct central meridian:
# https://www.lantmateriet.se/sv/Kartor-och-geografisk-information/GPS-och-geodetisk-matning/Referenssystem/Tvadimensionella-system/SWEREF-99-projektioner/
central_meridian = 15.00  			# For SWEREF99-1500

axis = 6378137.0  					# Semi-major axis of the ellipsoid. // GRS 80.
flattening = 1.0 / 298.257222101  	# Flattening of the ellipsoid. 		// GRS 80.
lat_of_origin = 0.0  				# Latitude of origin.
scale = 1.0  						# Scale on central meridian.
false_northing = 0.0  				# Offset for origo.
false_easting = 150000.0  			# Offset for origo.


# Conversion from grid coordinates to geodetic coordinates.
def grid_to_geodetic(x, y):
    """
	Input
		x - SWEREF X coordinate
		y - SWEREF Y coordinate
	Output
		lat_long - [latitude in degrees, longitude in degrees]
	"""
    e2 = flattening * (2.0 - flattening)
    n = flattening / (2.0 - flattening)
    a_roof = axis / (1.0 + n) * (1.0 + n ** 2 / 4.0 + n ** 4 / 64.0)
    delta1 = n / 2.0 - 2.0 * n ** 2 / 3.0 + 37.0 * n ** 3 / 96.0 - n ** 4 / 360.0
    delta2 = n ** 2 / 48.0 + n ** 3 / 15.0 - 437.0 * n ** 4 / 1440.0
    delta3 = 17.0 * n ** 3 / 480.0 - 37 * n ** 4 / 840.0
    delta4 = 4397.0 * n ** 4 / 161280.0
    astar = e2 + e2 ** 2 + e2 ** 3 + e2 ** 4
    bstar = -(7.0 * e2 ** 2 + 17.0 * e2 ** 3 + 30.0 * e2 ** 4) / 6.0
    cstar = (224.0 * e2 ** 3 + 889.0 * e2 ** 4) / 120.0
    dstar = -(4279.0 * e2 ** 4) / 1260.0

    lambda_zero = math.radians(central_meridian)
    xi = (y - false_northing) / (scale * a_roof)
    eta = (x - false_easting) / (scale * a_roof)

    xi_prim = xi - \
              delta1 * math.sin(2.0 * xi) * math.cosh(2.0 * eta) - \
              delta2 * math.sin(4.0 * xi) * math.cosh(4.0 * eta) - \
              delta3 * math.sin(6.0 * xi) * math.cosh(6.0 * eta) - \
              delta4 * math.sin(8.0 * xi) * math.cosh(8.0 * eta)
    eta_prim = eta - \
               delta1 * math.cos(2.0 * xi) * math.sinh(2.0 * eta) - \
               delta2 * math.cos(4.0 * xi) * math.sinh(4.0 * eta) - \
               delta3 * math.cos(6.0 * xi) * math.sinh(6.0 * eta) - \
               delta4 * math.cos(8.0 * xi) * math.sinh(8.0 * eta)
    phi_star = math.asin(math.sin(xi_prim) / math.cosh(eta_prim))
    delta_lambda = math.atan(math.sinh(eta_prim) / math.cos(xi_prim))

    lon_radian = lambda_zero + delta_lambda
    lat_radian = phi_star + math.sin(phi_star) * math.cos(phi_star) * \
                 (astar + bstar * math.pow(math.sin(phi_star), 2) +
                  cstar * math.pow(math.sin(phi_star), 4) +
                  dstar * math.pow(math.sin(phi_star), 6))

    # Conversion from radians to degrees
    return math.degrees(lat_radian), math.degrees(lon_radian)
