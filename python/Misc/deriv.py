"""
# numerical differentiation using finite differences, 2nd order accuracy.
# symmetric inside and asymmetric at edge or periodic boundary
# can do non-equidistant spacing in x
# y = f(x), returns dy/dx; all are numpy arrays
"""
import numpy as np


def deriv(y, x, periodic = False):
	n = y.size
	if(n < 3):
		raise RuntimeError('deriv: arrays must have at least 3 points')

	d = np.zeros(n)
	for i in range(1, n-1):  # inside the array
		d[i] = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1])

	if periodic:
		d[0] = 0.5*(y[1] - y[-1]) / (x[1] - x[0])
		d[-1] = 0.5*(y[0] - y[-2]) / (x[-1] - x[-2])
	else:
		if(abs(x[2] - 2*x[1] + x[0]) < 1e-12):		# equidistant x only
			d[0] = (-y[2] + 4*y[1] - 3*y[0]) / (x[2] - x[0])
		else:
			d[0] = (y[1] - y[0]) / (x[1] - x[0])

		if(abs(x[-1] - 2*x[-2] + x[-3]) < 1e-12):  # equidistant x only
			d[-1] = (y[-3] - 4*y[-2] + 3*y[-1]) / (x[-1] - x[-3])
		else:
			d[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])

	return d


def deriv5pt(y, x, periodic = False):
	n = y.size
	if(n < 6):
		raise RuntimeError('deriv: arrays must have at least 6 points')

	h = x[1] - x[0]	 # equidistant only
	d = np.zeros(n)
	for i in range(2, n-2):	 # inside the array
		d[i] = (-y[i+2] + 8*y[i+1] - 8*y[i-1] + y[i-2]) / (12.0*h)

	if periodic:
		d[1] = (-y[3] + 8*y[2] - 8*y[0] + y[-1]) / (12.0*h)
		d[0] = (-y[2] + 8*y[1] - 8*y[-1] + y[-2]) / (12.0*h)
		d[-1] = (-y[1] + 8*y[0] - 8*y[-2] + y[-3]) / (12.0*h)
		d[-2] = (-y[0] + 8*y[-1] - 8*y[-3] + y[-4]) / (12.0*h)
	else:
		d[1] = (-3*y[5] + 16*y[4] - 36*y[3] + 48*y[2] - 25*y[1]) / (12.0*h)
		d[0] = (-3*y[4] + 16*y[3] - 36*y[2] + 48*y[1] - 25*y[0]) / (12.0*h)
		d[-1] = (3*y[-5] - 16*y[-4] + 36*y[-3] - 48*y[-2] + 25*y[-1]) / (12.0*h)
		d[-2] = (3*y[-6] - 16*y[-5] + 36*y[-4] - 48*y[-3] + 25*y[-2]) / (12.0*h)
	return d


def deriv_halfGrid(y, x, periodic = False):
	import scipy.interpolate as scinter
	n = y.size
	if(n < 3):
		raise RuntimeError('deriv: arrays must have at least 3 points')

	d1 = np.zeros(n-1)
	xd1 = np.zeros(n-1)

	for i in range(0, n-1):
		xd1[i] = 0.5*(x[i+1] + x[i])
		d1[i] = (y[i+1] - y[i]) / (x[i+1] - x[i])

	f = scinter.UnivariateSpline(xd1, d1, s = 0)
	d = f(x)

	if periodic:
		d[0] = 0.5*(y[1] - y[-1]) / (x[1] - x[0])
		d[-1] = 0.5*(y[0] - y[-2]) / (x[-1] - x[-2])
	else:
		if(abs(x[2] - 2*x[1] + x[0]) < 1e-12):		# equidistant x only
			d[0] = (-y[2] + 4*y[1] - 3*y[0]) / (x[2] - x[0])
		else:
			d[0] = (y[1] - y[0]) / (x[1] - x[0])

		if(abs(x[-1] - 2*x[-2] + x[-3]) < 1e-12):	# equidistant x only
			d[-1] = (y[-3] - 4*y[-2] + 3*y[-1]) / (x[-1] - x[-3])
		else:
			d[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])
	return d
