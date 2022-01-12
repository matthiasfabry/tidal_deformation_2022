# %%
"""
Fitting of the single rotating star quantities, results computed here are written in app A
"""
import numpy as np
import scipy.integrate as spint
import scipy.optimize as spopt

import singlestarcheckingfuncs as single

pi = np.pi


def psi(x, om):
    """
    dimentionless potential of a rotating single star
    :param x: x coordinate
    :param om: fractional rotation rate
    :return: potential
    """
    return -1 / r(x, om) - x ** 2 / 2


def dpsidx(x, om):
    """
    dimensionless derivative of psi to x
    :param x: x coordinate
    :param om: fractonal rotation rate
    :return: d potential / dx
    """
    return x * (1 / r2(x, om) ** 1.5 - 1)


def dy2dx2(x, om):
    """
    second dimensionless derivative of y(x) of a rotating stellar shell
    :param x: x coordinate
    :param om: fractional rotation rate
    :return: dy/dx ^ 2
    """
    return x ** 2 / y2(x, om) * (alpha(x, om) ** -3 - 1) ** 2


def alpha(x, om):
    """
    auxiliary function
    :param x: x coordinate
    :param om: fractional rotation rate
    :return: alpha
    """
    return om ** (-2. / 3) + om ** (4. / 3) / 2 - x ** 2 / 2


def y(x, om):
    """
    y(x) curve of a rotating shell
    :param x: x coordinate
    :param om: fractional rotation rate
    :return: y(x)
    """
    return np.sqrt(y2(x, om))


def y2(x, om):
    """
    y^2(x) of a rotating shell
    :param x: x coordinate
    :param om: fractional rotation rate
    :return: y^2(x)
    """
    return alpha(x, om) ** -2 - x ** 2


def r(x, om):
    """
    radius length of a point on a rotating shell
    :param x: x coordinate of the point
    :param om: fractional rotation rate
    :return: r(x)
    """
    return np.sqrt(r2(x, om))


def r2(x, om):
    """
    r^2 of a point on a rotating shell
    :param x: x coordinate
    :param om: fractional rotation rate
    :return: r^2(x)
    """
    return x ** 2 + y2(x, om)


def dvoldx(x, om):
    """
    dimensionless derivative of the equipotential volume to x
    :param x: x coordinate
    :param om: fractional rotation rate
    :return: dV/dx
    """
    return 2 * pi * x * y(x, om)


def dsurfdx(x, om):
    """
    dimensionless derivative of the equipotential surface area to x
    :param x: x coordinate
    :param om: fractional rotation rate
    :return: dS/dx
    """
    return 2 * pi * x * np.sqrt(1 + dy2dx2(x, om))


def g(x, om):
    """
    dimensionless surface gravity at a point on a rotating shell
    :param x: x coordinate
    :param om: fractional rotation rate
    :return: g(x)
    """
    return 1 / r2(x, om) * np.sqrt(1 + x ** 2 * r(x, om) * (r2(x, om) ** 1.5 - 2))


def davgdx(x, om):
    """
    dimensionless derivative of the the average equipotential gravity to x
    :param x: x coordinate
    :param om: fractional rotation rate
    :return: d <g> / dx
    """
    return g(x, om) * dsurfdx(x, om)


def daveinvgdx(x, om):
    """
    dimensionless derivative of the the average equipotential inverse gravity to x
    :param x: x coordinate
    :param om: fractional rotation rate
    :return: d <g^-1> / dx
    """
    return 1. / g(x, om) * dsurfdx(x, om)


def dinerdx(x, om):
    """
    dimensionless derivative of the the average equipotential moment of inertia to x
    :param x: x coordinate
    :param om: fractional rotation rate
    :return: d i_rot / dx
    """
    return x ** 2 * daveinvgdx(x, om)


def volume(os):
    """
    computes numerically the dimensionless volumes of equipotentials with fractional rotation rates "os"
    for dimension full quantities scale this with the equatorial radius cubed.
    :param os: list, fractional rotation rates
    :return: V(os)
    """
    volm = np.zeros(len(os))
    for ii in range(len(os)):
        volm[ii] = 2 * spint.quad(dvoldx, 0, os[ii] ** (2 / 3), args=os[ii])[0] / (os[ii] ** 2)
    return volm


def surface(os):
    """
    computes numerically the dimensionless surface areas of equipotentials with fractional rotation rates "os"
    for dimension full quantities scale this with the equatorial radius squared.
    :param os: list, fractional rotation rates
    :return: S(os)
    """
    surf = np.zeros(len(os))
    for ii in range(len(os)):
        surf[ii] = 2 * spint.quad(dsurfdx, 0, os[ii] ** (2 / 3), args=os[ii])[0] / (os[ii] ** (4 / 3))
    return surf


def radius(os):
    """
    computes numerically the dimensionless volume equivalent radii of equipotentials with fractional rotation rates "os"
    for dimension full quantities scale this with the equatorial radius.
    :param os: list, fractional rotation rates
    :return: r_psi(os)
    """
    return np.cbrt(3 / (4 * pi) * volume(os))


def equatradius(os):
    """
    In this module we assume the equatorial radius of a shell with any rotation rate is 1.
    for dimension full quantities scale this with the equatorial radius.
    :param os: list, fractional rotation rates
    :return: list of ones
    """
    return np.ones(len(os))  # everything is scaled to this value


def surface_times_av_gravity(os):
    """
    computes numerically the dimensionless surface area times average gravity of equipotentials
    with fractional rotation rates "os"
    for dimension full quantities scale this with GM_Psi.
    :param os: list, fractional rotation rates
    :return: S(os) * <g>(os)
    """
    surftgrav = np.zeros(len(os))
    for ii in range(len(os)):
        surftgrav[ii] = 2 * spint.quad(davgdx, 0, os[ii] ** (2 / 3), args=os[ii])[0]
        # here should be times GM
    return surftgrav


def surface_times_av_inv_gravity(os):
    """
    computes numerically the dimensionless surface area times average inverse gravity of equipotentials
    with fractional rotation rates "os"
    for dimension full quantities scale this with 1/GM_Psi.
    :param os: list, fractional rotation rates
    :return: S(os) * <g^-1>(os)
    """
    surftinvgrav = np.zeros(len(os))
    for ii in range(len(os)):
        surftinvgrav[ii] = 2 * spint.quad(daveinvgdx, 0, os[ii] ** (2 / 3), args=os[ii])[0] / os[
            ii] ** (8 / 3)
        # here should be div by GM
    return surftinvgrav


def fP(os):
    """
    computes numerically the pressure correction factor of equipotentials
    with fractional rotation rates "os"
    this quantity is dimensionless by construction.
    :param os: list, fractional rotation rates
    :return: f_P(os)
    """
    return 4 * pi * radius(os) ** 4 / surface_times_av_inv_gravity(os)


def fT(os):
    """
    computes numerically the temperature correction factor of equipotentials
    with fractional rotation rates "os"
    this quantity is dimensionless by construction.
    :param os: list, fractional rotation rates
    :return: f_T(os)
    """
    return (4 * pi * radius(os) ** 2) ** 2 / (surface_times_av_inv_gravity(os) * surface_times_av_gravity(os))


def inertia(os):
    """
    computes numerically the dimensionless moment of inertia of equipotentials
    with fractional rotation rates "os"
    for dimension full quantities scale this with the equatorial radius squared.
    :param os: list, fractional rotation rates
    :return: i_rot(os)
    """
    ints = np.zeros(len(os))
    for ii in range(len(os)):
        ints[ii] = 2 * spint.quad(dinerdx, 0, os[ii] ** (2 / 3), args=os[ii])[0] / (os[ii] ** 4)
    ints /= surface_times_av_inv_gravity(os)
    return ints


def jrotgmr(os):
    """
    computes numerically the dimensionless fractional specific angular momentum of equipotentials
    with fractional rotation rates "os"
    :param os: list, fractional rotation rates
    :return: j_rot / sqrt(G M_psi r_psi)(os)
    """
    return os * inertia(os) * np.sqrt(1 / radius(os))


def polyvol(om, a):
    """
    describes the polynomial fit to V_Psi
    :param om: fractional rotation rate
    :param a: coefficient of om^4
    :return: fitted value of V_Psi
    """
    return 4 * pi / 3 * (1 - 0.5 * om ** 2 + a * om ** 4 + vol_boundary(a) * om ** 6)


def polysurf(om, a):
    """
    describes the polynomial fit to S_Psi
    :param om: fractional rotation rate
    :param a: coefficient of om^4
    :return: fitted value of S_Psi
    """
    return 4 * pi * (1 - om ** 2 / 3 + a * om ** 4 + surf_boundary(a) * om ** 6)


def polyrad(om, a):
    """
    describes the polynomial fit to r_Psi
    :param om: fractional rotation rate
    :param a: coefficient of om^4
    :return: fitted value of r_Psi
    """
    return 1 - om ** 2 / 6 + a * om ** 4 + rad_boundary(a) * om ** 6


def polyeqrad(om, a):
    """
    describes the polynomial fit to r_e
    since r_e is assumed 1, the fit of f in found by solving implicitly 1 = r_psi(om) * f(om)
    :param om: fractional rotation rate
    :param a: coefficient of om^4
    :return: fitted value of r_e
    """
    return radius(om) * (1 + om ** 2 / 6 + a * om ** 4 + eqrad_boundary(a) * om ** 6)


def polysurftimesgrav(om, a):
    """
    describes the polynomial fit to S_Psi <g>
    :param om: fractional rotation rate
    :param a: coefficient of om^4
    :return: fitted value of S_Psi <g>
    """
    return 4 * pi * (1 - 2 / 3 * om ** 2 + a * om ** 4 + surftimesgrav_boundary(a) * om ** 6)


def polysurftimesinvgrav(om, a, b, c, d):
    """
    describes the polynomial fit to S_Psi <g^-1>
    :param om: fractional rotation rate
    :param a: coefficient of om^4 in A
    :param b: coefficient of om^6 in A
    :param c: coefficient of the logarithmic term in A
    :param d: exponent of om within the logarithm in A
    :return: fitted value of S_Psi <g^-1>
    """
    return 4 * pi * A(om, a, b, c, d)


def A(om, a, b, c, d):
    """
    Auxiliary function in the polynomial fit to S_Psi <g^-1>
    :param om: fractional rotation rate
    :param a: coefficient of om^4
    :param b: coefficient of om^6
    :param c: coefficient of the logarithmic term
    :param d: exponent of om within the logarithm
    :return: A
    """
    return 1 + a * om ** 4 + b * om ** 6 + c * np.log1p(-om ** d)


def polyinertia(om, a):
    """
    describes the polynomial fit to i_rot
    requires fit of S <g^-1>
    :param om: fractional rotation rate
    :param a: coefficient of om^4 in B
    :param b: coefficient of om^6 in B
    :return: fitted value of i_rot
    """
    gminonefit = fits['surface*gravity-1'].get_fit_coeffs()
    return 2 / 3 * B(om, a, gminonefit[0], gminonefit[1], gminonefit[2], gminonefit[3]) / A(om, *gminonefit)


def B(om, a, gma, gmb, c, d):
    """
    Auxiliary function in the polynomial fit to i_rot
    :param om: fractional rotation rate
    :param a: coefficient of om^4
    :param gma: coefficient of om^4 in the S_Psi <g^-1> fit
    :param gmb: coefficient of om^6 in the S_Psi <g^-1> fit
    :param c: coefficient of the logarithmic term
    :param d: exponent of om within the logarithm
    :return: B
    """
    return 1 + 1 / 5 * om ** 2 + a * om ** 4 + inertia_boundary(a, gma, gmb) * om ** 6 + 3 / 2 * c * np.log1p(-om ** d)


def polyjrot(om, a, b):
    """
    describes the polynomial fit to j_rot / sqrt(G M_Psi r_psi)
    the fit to S <g^-1> has to be computed first in the dictionary "fits"
    :param om: fractional rotation rate
    :param a: coefficient of om^4
    :param b: coefficient of om^6
    :return: fitted value of j_rot / sqrt(G M_Psi r_psi)
    """
    gminonefit = fits['surface*gravity-1'].get_fit_coeffs()
    return 2 / 3 * om * C(om, a, b, gminonefit[2], gminonefit[3]) / A(om, *gminonefit)


def C(om, a, b, c, d):
    """
    Auxiliary function in the polynomial fit to j_rot
    :param om: fractional rotation rate
    :param a: coefficient of om^4
    :param b: coefficient of om^6
    :param c: coefficient of the logarithmic term of the S_Psi <g^-1> fit
    :param d: exponent of om within the logarithm of the S_Psi <g^-1> fit
    :return: C
    """
    return 1 + 17 / 60 * om ** 2 + a * om ** 4 + b * om ** 6 + jrot_boundary(c) * np.log1p(-om ** d)


def polyfP(om, a, b):
    """
    describes the polynomial fit to f_P
    the fit to S <g^-1> has to be computed first in the dictionary "fits"
    :param om: fractional rotation rate
    :param a: coefficient of om^4
    :param b: coefficient of om^6
    :return: fitted value of f_P
    """
    gminonefit = fits['surface*gravity-1'].get_fit_coeffs()
    return (1 - 2 / 3 * om ** 2 + a * om ** 4 + b * om ** 6) / A(om, *gminonefit)


def polyfT(om, a, b):
    """
    describes the polynomial fit to f_T
    the fit to S <g^-1> has to be computed first in the dictionary "fits"
    :param om: fractional rotation rate
    :param a: coefficient of om^4
    :param b: coefficient of om^6
    :return: fitted value of f_T
    """
    gminonefit = fits['surface*gravity-1'].get_fit_coeffs()
    return (1 + a * om ** 4 + b * om ** 6) / A(om, *gminonefit)


def vol_boundary(a):
    """
    provides the V_Psi(om=1) boundary condition for the om^6 coefficient as function of the om^4 coefficient
    :param a: coefficient om om^4 in the V_Psi fit
    :return: coefficient of om^6
    """
    return -a - .5 + volumeatone * 3 / (4 * pi)


def surf_boundary(a):
    """
    provides the S_Psi(om=1) boundary condition for the om^6 coefficient as function of the om^4 coefficient
    :param a: coefficient om om^4 in the S_Psi fit
    :return: coefficient of om^6
    """
    return -a - 2 / 3 + surfaceatone / (4 * pi)


def rad_boundary(a):
    """
    provides the r_psi(om=1) boundary condition for the om^6 coefficient as function of the om^4 coefficient
    :param a: coefficient om om^4 in the r_Psi fit
    :return: coefficient of om^6
    """
    return -a - 5 / 6 + radiusatone


def eqrad_boundary(a):
    """
    provides the r_e(om=1) boundary condition for the om^6 coefficient as function of the om^4 coefficient
    :param a: coefficient of om^4 in the r_e fit
    :return: coefficient of om^6
    """
    return -a - 7 / 6 + 1 / radiusatone


def surftimesgrav_boundary(a):
    """
    provides the S_Psi <g>(om=1) boundary condition for the om^6 coefficient as function of the om^4 coefficient
    :param a: coefficient of om^4 in the S_Psi <g> fit
    :return: coefficient of om^6
    """
    return -a - 1 / 3 + surftimesgravityatone / (4 * pi)


def inertia_boundary(a, gma, gmb):
    """
    provides the i_rot(om=1) boundary condition for the om^6 coefficient as function of the om^4 coefficient
    :param a: coefficient om om^4 in the i_rot fit
    :param gma: coefficient of om^4 in the S_Psi <g^-1> fit
    :param gmb: coefficient of om^6 in the S_Psi <g^-1> fit
    :return: coefficient of om^6
    """
    return -a - 6 / 5 + 3 / 2 * (1 + gma + gmb)


def jrot_boundary(gmc):
    """
    provides the j_rot(om=1) boundary condition for the om^6 coefficient as function of the om^4 coefficient
    :param gmc: coefficient of the logarithmic term the S_Psi fit
    :return: coefficient of om^6
    """
    return 1.5 * gmc / np.sqrt(radiusatone)


thresh = 0.001
# compute the dimentionless values of V_Psi, S_psi, r_psi, S_psi<g>, S_Psi<g^-1> at critical rotation
volumeatone = 2 * spint.quad(dvoldx, 0, 1, args=1.)[0]
surfaceatone = 2 * spint.quad(dsurfdx, 0, 1, args=1.)[0]
radiusatone = np.cbrt(3 / (4 * pi) * volumeatone)
surftimesgravityatone = 2 * spint.quad(davgdx, 0, 1, args=1.)[0]
surftimesinvgravityatone = np.inf

# protect the following against importing
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    
    
    class SinglePolynomializer:
        """
        Objectifies the curve-fit procedure to find the coefficients of the polynomial fits to the single
        rotating star quantities.
        """
        
        def __init__(self, name, os, fun_exact, fun_poly, powers, fittedfun=None):
            """
            constructor
            :param name: string to indentify this object
            :param os: range of fractional rotation rates you want to fit over
            :param fun_exact: python function that return the exact values you are trying to fit
            :param fun_poly: python function that parametrically computes the polynomial fit
            :param powers: list of exponents you want to use as error weightings in the fitting procedure (1/om^power is
            used)
            :param fittedfun: (optional) if you already have a fitted function, provide it here (eg for quick plotting)
            """
            self._name = name
            self._powers = powers
            self._mini = 0
            self._max_rel_distance = np.inf
            self._coeffs = []
            self._o = os
            self._exacts = fun_exact(self._o)
            # self._exacts = np.load(self._name+".npy")  # if you have the numerical results saved
            # np.save(self._name, self._exacts)  # save numerical results for quick loading later
            self._fun_poly = fun_poly
            self._fitted_fun = fittedfun
        
        def get_exacts(self):
            """
            returns the exact values
            :return: self.exacts
            """
            return self._exacts
        
        def get_fit_coeffs(self):
            """
            returns the best fit coefficients of the polynomial fit
            """
            return self._coeffs[self._mini]
        
        def fit_polynomial(self):
            """
            use scipy.curve_fit to to fit the function self._fun_poly to the exact values self._exacts.
            Do this for all self._powers weightings and save the maximal relative error and at which index in
            self._powers this happens.
            """
            for i in range(len(self._powers)):
                self._coeffs.append(
                    spopt.curve_fit(self._fun_poly, self._o, self._exacts, sigma=1. / self._o ** self._powers[i])[0])
                if max(np.abs(1 - self._fun_poly(self._o, *self._coeffs[i]) / self._exacts)) < self._max_rel_distance:
                    self._max_rel_distance = max(np.abs(1 - self._fun_poly(self._o, *self._coeffs[i]) / self._exacts))
                    self._mini = i
                    
            # this print to the console gives you the coefficients of the polynomials you want to save/write down
            print(self._name, self._powers[self._mini], self._coeffs[self._mini], self._max_rel_distance)
        
        def plot_rel_error(self):
            plt.plot(self._o, 1 - self._fun_poly(self._o, *self._coeffs[self._mini]) / self._exacts,
                     label=self._name + ' error')
        
        def plot_polyfit(self):
            plt.plot(o, self._fun_poly(o, *self._coeffs[self._mini]), label=self._name + ' poly')
        
        def plot_exacts(self):
            plt.plot(self._o, self._exacts, label=self._name + ' exact')
        
        def already_fitted_plots(self):
            """
            make a quick plot and compute the maximal relative error of fitted_fun
            :return:
            """
            diff = np.abs(1 - self._fitted_fun(self._o) / self._exacts)
            print(self._name, max(diff))
            plt.plot(self._o, self._fitted_fun(self._o), label=self._name + ' fit')
            plt.plot(self._o, self._exacts, label=self._name + " exacts")
    
    
    o = np.linspace(0.0001, 1, 10000)
    o_not_one_not_zero = np.linspace(0.02, 0.9999, 10000)
    
    print('critical star volume =', volumeatone)
    print('critical star surface =', surfaceatone)
    print('critical star equiv radius =', radiusatone)
    print('critical star surf*grav=', surftimesgravityatone)
    print('critical star average gravity=', surftimesgravityatone / surfaceatone)
    
    # noinspection PyDictCreation
    fits = {}
    
    # you can play with the powers to see if you can get ever better fits
    fits['volume'] = SinglePolynomializer('vol', o, volume, polyvol, [4], single.polyvol)
    fits['surface'] = SinglePolynomializer('surf', o, surface, polysurf, [2, 3, 4, 5, 6], single.polysurf)
    fits['radius'] = SinglePolynomializer('pot_rad', o, radius, polyrad, [3], single.polyequivrad)
    fits['equatradius'] = SinglePolynomializer('eq_rad', o, equatradius, polyeqrad, [2], single.polyequarad)
    fits['surface*gravity'] = SinglePolynomializer('surf*g', o, surface_times_av_gravity, polysurftimesgrav, [2],
                                                   single.polysurftimesgrav)
    fits['surface*gravity-1'] = SinglePolynomializer('surf*g-1', o_not_one_not_zero, surface_times_av_inv_gravity,
                                                     polysurftimesinvgrav, [7], single.polysurftimesinvgrav)
    fits['inertia'] = SinglePolynomializer('inertia', o_not_one_not_zero, inertia, polyinertia, [4], single.polyirot)
    fits['fP'] = SinglePolynomializer('fP', o_not_one_not_zero, fP, polyfP, [5], single.polyfP)
    fits['fT'] = SinglePolynomializer('fT', o_not_one_not_zero, fT, polyfT, [4], single.polyfT)
    fits['jrot'] = SinglePolynomializer('jrotovergmr', o_not_one_not_zero, jrotgmr, polyjrot, [2], single.polyjrot)
    
    plt.figure()
    for key, qqt in fits.items():
        qqt.fit_polynomial()  # heavy lifting here
        qqt.plot_rel_error()
    plt.legend()
    plt.figure()
    for key, qqt in fits.items():
        qqt.plot_exacts()
        qqt.plot_polyfit()
    plt.legend()
    plt.show()
