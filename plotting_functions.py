from numpy import *
from pylab import *
import os
import time
from shutil import copyfile
from scipy.optimize import leastsq

src_path = os.getcwd()
# useinj = False

needfix = ['EC235', 'EC236', 'EC237', 'EC238', 'EC239']


def ctim1(x):  # modified
    return int(float(x.split('\t')[1]))


def ctim(x):
    a = x.split('\t')[1].split(':')
    return float(a[0]) * 3600 + float(a[1]) * 60 + float(a[2])


def timearray(filen):
    Timefile = open(filen + '.txt', mode='r')
    b = Timefile.readlines()
    Timefile.close()
    ret = []
    for i in range(1, len(b)):
        ret.append(ctim1(b[i]))
    return array(ret)  # - ret[0]


def data(filen):
    candt = loadtxt(filen + '.txt', skiprows=1, usecols=(2, 3))
    Conductivity = candt[:, 0]
    Temperature = candt[:, 1]
    timme = timearray(filen)
    return [timme, Conductivity, Temperature]


def selecdat(dat, start, stop):
    return [dat[0][start:stop], dat[1][start:stop], dat[2][start:stop]]


def inj_time(onelist):
    injfile = open(onelist + '.txt', mode='r')
    c = injfile.readlines()
    injfile.close()
    ter = {}
    for j in range(len(c)):
        ter[c[j].split('\t')[0]] = ctim(c[j])
    return ter


def imp_data(src_path):
    src_objnames = os.listdir(src_path)
    store = {}
    for name in src_objnames:
        if name[-4:] == '.txt' and name[:3] != 'inj':
            store[name[:-4]] = data(name[:-4])
    inj = inj_time('injection_times')
    return [store, inj]


def idata_nooffset(src_path):
    IMPD = imp_data(src_path)
    for k in IMPD[0]:
        offset = IMPD[0][k][0][0]
        IMPD[0][k][0] = IMPD[0][k][0] - offset
        if k in IMPD[-1]:
            IMPD[1][k] = IMPD[1][k] - offset
    return IMPD


def allplot():
    impd = imp_data(src_path)
    store = impd[0]
    inj = impd[1]
    i = 1
    for k in store:
        figure(i)
        plot(store[k][0] - store[k][0][0], store[k][1], '.', markersize=2,
             label=k)
        if k in inj:
            injx = ones(10) * inj[k] - store[k][0][0]
            injy = linspace(min(store[k][1]), max(store[k][1]), 10)
            plot(injx, injy, 'rx')
        ylabel('Absorbance')
        xlabel('Time (s)')
        grid(which='major', linewidth=0.2, ls='-', alpha=0.8)
        legend()
        i += 1
    show()


def peval_lin(p, x):
    return p[0] + p[1] * x


def residuals_lin(p, y, x):
    err = y - peval_lin(p, x)
    return err


def peval_MM(p, x):  # p[0] = Vmax, p[1] = Km
    return p[0] * x / (p[1] + x)


def residuals_MM(p, y, x):
    err = y - peval_MM(p, x)
    return err


def fit_lin(dat, p0):
    p_final, cov_x, info, mesg, success = leastsq(residuals_lin, p0, args=
    (dat[1], dat[0]), full_output=True, maxfev=10000)
    avg = sum(dat[1]) / len(dat[1])
    sstot = sum((dat[1] - avg) ** 2)
    ssres = sum((dat[1] - peval_lin(p_final, dat[0])) ** 2)
    Rsquared = 1 - ssres / sstot
    return p_final, Rsquared


def fit_MM(dat, p0):
    p_final, cov_x, info, mesg, success = leastsq(residuals_MM, p0, args=
    (dat[1], dat[0]), full_output=True, maxfev=10000)
    avg = sum(dat[1]) / len(dat[1])
    sstot = sum((dat[1] - avg) ** 2)
    ssres = sum((dat[1] - peval_MM(p_final, dat[0])) ** 2)
    Rsquared = 1 - ssres / sstot
    return p_final, Rsquared


def looplinfit(uranges, impd):
    p0 = 80, 0.01
    linfits = {}
    for k in uranges:
        if k in needfix:
            dat = selecdat(impd[0][k], uranges[k][1] * 2 + 1, uranges[k][2] * 2 + 1)
        else:
            dat = selecdat(impd[0][k], uranges[k][1], uranges[k][2])
        param = fit_lin(dat, p0)
        linfits[k] = [param, uranges[k][0]]
    return linfits


def allplot_linfit(uranges, impd):
    linfs = looplinfit(uranges, impd)
    store = impd[0]
    inj = impd[1]
    i = 1
    for k in uranges:
        x = linspace(uranges[k][1], uranges[k][2], 200)
        param = linfs[k][0][0]
        figure(i)
        plot(store[k][0], store[k][1], '.', markersize=2, label=k
                                                                + ', ' + str(uranges[k][0]) + 'mM')  # uL if enz linear
        plot(x, peval_lin(param, x), label='fit: yint, slope, R**2 =' + '\n'
                                           + str(round(linfs[k][0][0][0], 4)) + ', '
                                           + str(round(linfs[k][0][0][1], 4)) + ', '
                                           + str(round(linfs[k][0][1], 4)))
        if k in inj:
            injx = ones(10) * inj[k]
            injy = linspace(min(store[k][1]), max(store[k][1]), 10)
            plot(injx, injy, 'rx')
        ylabel('Absorbance')
        xlabel('Time (s)')
        grid(which='major', linewidth=0.2, ls='-', alpha=0.8)
        legend()
        savefig('ratevconc ' + str(uranges[k][0]) + 'mM' + '.pdf')  # uL in enz linear
        i += 1
    show()


def LBfit(fits):
    invconc = []
    invrate = []
    for k in fits:
        invconc.append(fits[k][-1] ** (-1))
        invrate.append(fits[k][0][0][1] ** (-1))
    dats = [array(invconc), array(invrate)]
    p0 = -1, 20
    return fit_lin(dats, p0)


def MMfit(fits):
    conc = []
    rate = []
    for k in fits:
        conc.append(fits[k][-1])
        rate.append(-fits[k][0][0][1])
    dats = [array(conc), array(rate)]
    p0 = 0.05, 1.5
    return fit_MM(dats, p0)


def enzlin(fits):
    enzconc = []
    rate = []
    for k in fits:
        enzconc.append(fits[k][-1] * 17 / 42)
        rate.append(fits[k][0][0][1])
    dats = [array(enzconc), array(rate)]
    p0 = -1, 20
    return fit_lin(dats, p0)


def lbcustomfit(fits):
    x = linspace(-3, 6, 200)
    MMf = MMfit(fits)
    invparams = array([1 / MMf[0][0], MMf[0][1] / MMf[0][0]])
    plot(x, peval_lin(invparams, x), label='Vmax, Km, R**2 = ' + '\n'
                                           + str(round(MMf[0][0], 4)) + ' uS/cm*s, '
                                           + str(round(MMf[0][1], 4)) + ' mM, '
                                           + str(round(MMf[1], 3)))
    figure(1)
    for k in fits:
        plot(fits[k][-1] ** (-1), (-fits[k][0][0][1]) ** (-1), 'bo')
    ylabel('1 / d/dt Absorbance)')
    xlabel('1 / Concentration (1/mM)')
    title('Lineweaver-Burk Plot (Parameters from MM Fit)')
    grid(which='major', linewidth=0.2, ls='-', alpha=0.8)
    legend()
    savefig('L-B_fitfromMM' + '.pdf')
    show()


def lineweaverburk(fits):
    lbf = LBfit(fits)
    x = linspace(-3, 6, 200)
    plot(x, peval_lin(lbf[0], x), label='Vmax, Km, R**2 = ' + '\n'
                                        + str(round(lbf[0][0] ** (-1), 4)) + ' uS/cm*s, '
                                        + str(round(lbf[0][1] / lbf[0][0], 4)) + ' mM, '
                                        + str(round(lbf[1], 3)))
    figure(1)
    for k in fits:
        plot(fits[k][-1] ** (-1), fits[k][0][0][1] ** (-1), 'bo')
    ylabel('1 / d/dt Absorbance)')
    xlabel('1 / Concentration (1/mM)')
    title('Lineweaver-Burk Plot for M-M Kinetics')
    grid(which='major', linewidth=0.2, ls='-', alpha=0.8)
    legend()
    savefig('lineweaverburk' + '.pdf')
    show()


def ratevconc(fits, run=0):
    MMf = MMfit(fits)
    x = linspace(0, 11, 200)
    plot(x, peval_MM(MMf[0], x), label='Vmax, Km, R**2 = ' + '\n'
                                       + str(round(MMf[0][0], 4)) + ' uS/cm*s, '
                                       + str(round(MMf[0][1], 4)) + ' mM, '
                                       + str(round(MMf[1], 3)))
    figure(1)
    for k in fits:
        plot(fits[k][-1], -fits[k][0][0][1], 'bo')
    ylabel('d/dt Absorbance)')
    xlabel('Concentration (mM)')
    title('Rate vs Concentration with M-M Fit')
    grid(which='major', linewidth=0.2, ls='-', alpha=0.8)
    legend()
    savefig('run #' + str(run) + ' ratevconc.pdf')
    show()


def ratevconc_enzlin(fits):
    ezlinear = enzlin(fits)
    x = linspace(0, 35, 200)
    plot(x, peval_lin(ezlinear[0], x), label='Yint, slope, R**2 = ' + '\n'
                                             + str(round(ezlinear[0][0], 4)) + ' uS/cm*s, '
                                             + str(round(ezlinear[0][1], 4)) + ' (uS*/cm*s)/(U/mL), '
                                             + str(round(ezlinear[1], 3)))
    figure(1)
    for k in fits:
        plot(fits[k][-1] * 17 / 42, fits[k][0][0][1], 'bo')  # , label=k)
    ylabel('d/dt Absorbance)')
    xlabel('Concentration (U/mL)')
    title('Rate vs Concentration')
    grid(which='major', linewidth=0.2, ls='-', alpha=0.8)
    legend()
    savefig('enzyme_linearity.pdf')
    show()

