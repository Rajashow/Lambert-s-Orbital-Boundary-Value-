from numpy import linalg as LA
import numpy as np
from scipy.optimize import root_scalar


def get_init_guess(lambda_, T, number_rev=0):
    # single rev
    if number_rev == 0:
        T_0 = np.arccos(lambda_) + lambda_ * np.sqrt(1 - lambda_ ** 2) + \
            number_rev * np.pi
        T_1 = 2 * (1 - lambda_ ** 3) / 3
        if T >= T_0:
            x_0 = (T_0 / T) ** (2 / 3) - 1
        elif T < T_1:
            x_0 = 5 / 2 * T_1 / T * (T_1 - T) / (1 - lambda_ ** 5) + 1
        else:

            x_0 = (T_0 / T) ** (np.log2(T_1 / T_0)) - 1

        return (x_0)
    else:
        x_0l = (((number_rev * np.pi + np.pi) / (8 * T)) ** (2 / 3) - 1) / (
            ((number_rev * np.pi + np.pi) / (8 * T)) ** (2 / 3) + 1
        )
        x_0r = (((8 * T) / (number_rev * np.pi)) ** (2 / 3) - 1) / (
            ((8 * T) / (number_rev * np.pi)) ** (2 / 3) + 1
        )

        return (x_0l, x_0r)


def find_xy(lambda_, T):
    m_max = T//np.pi
    T_00 = np.arccos(lambda_ + lambda_*np.sqrt(1-(np.power(lambda_, 2))))
    if T < T_00+(m_max*np.pi) and m_max > 0:
        t_min = None
        def f(x): raise ValueError("not done")
        v = root_scalar(f, x0=0, args={"T": T_00})
        # start Halley iterations fromx= 0,T=T0 and findTmin(Mmax)
        if t_min > T:
            m_max = m_max-1
    T_1 = (1-np.power(lambda_, 3))*2/3
    # compute x0
    # start Householder iterations fromx0and findx,y
    while m_max > 0:
        # compute x0l and x0r from Eq.(31) with M=M_max

            # start Householder iterations from x0land find xr,yr
            # start Householder iterations from x0rand find xl,yl
        x_0s = get_init_guess(lambda_, T, 0)
        for x_0 in x_0s:
            #   x and y
            pass
        m_max = m_max-1


def solve_lambert(r1, r2, t, mu):
    if not (t > 0 and mu > 0):
        #     throw error
        pass
    if np.all(np.cross(r1, r2) == 0):

        pass
    c = r1-r2
    c = LA.norm(c)
    norm_r1 = LA.norm(r1)
    r2_ = LA.norm(r2)
    s = (norm_r1+r2_+c)/2
    i_r_1, i_r_2 = r1/norm_r1, r2/r2_
    i_h = np.cross(i_r_1, i_r_2)
    lambda_ = np.sqrt(1 - (c/s))
    if (np.dot(r1[0], r1[1]) - np.dot(r2[0], r2[1])) < 0:
        lambda_ = -lambda_
        i_t_1 = np.cross(i_r_1, i_h)
        i_t_2 = np.cross(i_r_2, i_r_2)
    else:
        i_t_1 = np.cross(i_r_1, i_h)
        i_t_2 = np.cross(i_r_2, i_r_2)
    T = np.sqrt(2*mu / np.power(s, 3))*t
    xy_s = find_xy(lambda_, T)
    gamma = np.sqrt(mu*s/2)
    rho = (norm_r1 - r2_)/c
    sigma = np.sqrt(1-rho**2)

    for x, y in xy_s:
        # Reconstruct solution velocity vectors
        num = gamma * (((lambda_*y) - x) - rho*(lambda_*y + x))
        v_r_1 = num/norm_r1
        v_r_2 = -num/r2_
        vt_num = (y + lambda_*x)*gamma*sigma
        v_t1 = vt_num/r1
        v_t2 = vt_num/r2

        v1 = v_r_1*i_r_1 + v_t1*i_t_1
        v2 = v_r_2*i_r_2 + v_t2*i_t_2

        yield v1, v2
