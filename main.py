#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 09:54:51 2017
@author: lihu
"""

import numpy as np
from numpy.ctypeslib import load_library, ndpointer
from ctypes import c_int

# params used by main program
pi    = 3.14159265359
omega = 2*pi/(24*3600)
theta = 35.0      # 区域所在纬度
nt    = 288000
nout  = 2400

txMAX = 0.2
tyMAX = 0.2

# params used by dll
nx    = 51
ny    = 51

# struct arg_params
dx    = 100.0
dy    = 100.0
dt    = 3.0
hmin  = 0.05
rho   = 1028.0
g     = 9.81

f     = 2*omega*np.sin(theta*pi/180)  # 科氏力参数
beta  = (0.5*dt*f)**2  # beta系数
r     = 0.001
taux  = 0.0
tauy  = 0.0
ah    = 1.0

mode  = 1
slip  = 0.0

arg_params = np.dtype([('dx',   'float64', 1), \
                       ('dy',   'float64', 1), \
                       ('dt',   'float64', 1), \
                       ('hmin', 'float64', 1), \
                       ('rho',  'float64', 1), \
                       ('g',    'float64', 1), \
                       ('f',    'float64', 1), \
                       ('beta', 'float64', 1), \
                       ('r',    'float64', 1), \
                       ('taux', 'float64', 1), \
                       ('tauy', 'float64', 1), \
                       ('ah',   'float64', 1), \
                       ('mode', 'int32',   1), \
                       ('slip', 'float64', 1)])
params = np.array([(dx, dy, dt, hmin, rho, g, f, beta, r, taux, tauy, ah, mode, slip)], dtype=arg_params)

#%% initialize var
hs      = np.zeros((nx+2, ny+2), dtype='float64', order='F')
eta_cur = np.zeros((nx+2, ny+2), dtype='float64', order='F')
eta_nex = np.zeros((nx+2, ny+2), dtype='float64', order='F')   
h       = np.zeros((nx+2, ny+2), dtype='float64', order='F')
mask    = np.zeros((nx+2, ny+2), dtype='float64', order='F')
u_cur   = np.zeros((nx+2, ny+2), dtype='float64', order='F')
u_nex   = np.zeros((nx+2, ny+2), dtype='float64', order='F')
v_cur   = np.zeros((nx+2, ny+2), dtype='float64', order='F')
v_nex   = np.zeros((nx+2, ny+2), dtype='float64', order='F')

# 读取净水深文件
inner_dep = np.loadtxt('./topo1.dat')
inner_dep = inner_dep.T
hs[1:-1, 1:-1] = inner_dep[:, :]
# 设置外围边界的净水深
hs[0, :]  = -0.5
hs[-1, :] = -0.5
hs[:, 0]  = -0.5
hs[:, -1] = -0.5

eta_cur[:, :] = -np.minimum(hs, 0.0)
eta_nex[:, :] = eta_cur[:, :]

h[:, :] = hs[:, :]+eta_cur[:, :]
mask[:, :] = 1
mask[np.where(h[:, :]<hmin)] = 0

u_cur[:, :] = 0.0
u_nex[:, :] = 0.0
v_cur[:, :] = 0.0
v_nex[:, :] = 0.0

#%% load dll
fdll = load_library('proc', './')

# cal_u args setup
fdll.cal_u.argtypes = [ndpointer(dtype='float64', ndim=2, flags='F'), \
                       ndpointer(dtype='float64', ndim=2, flags='F'), \
                       ndpointer(dtype='float64', ndim=2, flags='F'), \
                       ndpointer(dtype='float64', ndim=2, flags='F'), \
                       ndpointer(dtype='float64', ndim=2, flags='F'), \
                       ndpointer(dtype='float64', ndim=2, flags='F'), \
                       c_int,                                         \
                       c_int,                                         \
                       ndpointer(dtype=arg_params)]

# cal_v args setup
fdll.cal_v.argtypes = [ndpointer(dtype='float64', ndim=2, flags='F'), \
                       ndpointer(dtype='float64', ndim=2, flags='F'), \
                       ndpointer(dtype='float64', ndim=2, flags='F'), \
                       ndpointer(dtype='float64', ndim=2, flags='F'), \
                       ndpointer(dtype='float64', ndim=2, flags='F'), \
                       ndpointer(dtype='float64', ndim=2, flags='F'), \
                       c_int,                                         \
                       c_int,                                         \
                       ndpointer(dtype=arg_params)]

# cal_eta args setup
fdll.cal_eta.argtypes = [ndpointer(dtype='float64', ndim=2, flags='F'), \
                         ndpointer(dtype='float64', ndim=2, flags='F'), \
                         ndpointer(dtype='float64', ndim=2, flags='F'), \
                         ndpointer(dtype='float64', ndim=2, flags='F'), \
                         ndpointer(dtype='float64', ndim=2, flags='F'), \
                         c_int,                                         \
                         c_int,                                         \
                         ndpointer(dtype=arg_params)]

# open output file
eta_fid = open('./output_eta.txt', 'w+')
u_fid = open('./output_u.txt', 'w+')
v_fid = open('./output_v.txt', 'w+')

# calculate
for n in np.arange(1, nt+1):
    time = n*dt
    params['taux'] = txMAX*np.minimum(time/(1.0*24.0*3600.0), 1.0)
    params['tauy'] = tyMAX*np.minimum(time/(1.0*24.0*3600.0), 1.0)
    
    fdll.cal_u(mask, h, eta_cur, u_cur, v_cur, u_nex, nx, ny, params)
    fdll.cal_v(mask, h, eta_cur, u_cur, v_cur, v_nex, nx, ny, params)
    fdll.cal_eta(h, eta_cur, u_nex, v_nex, eta_nex, nx, ny, params)
    
    eta_cur[:, :] = eta_nex[:, :]
    h[:, :]       = hs[:, :]+eta_cur[:, :]
    mask[:, :]    = 1
    mask[np.where(h[:, :]<hmin)] = 0
    u_cur[:, :]   = u_nex[:, :]
    v_cur[:, :]   = v_nex[:, :]
    
    # output
    if n%nout == 0:
        np.savetxt(eta_fid, eta_cur[:, :], fmt='%10.6f', delimiter=' ', newline='\n')
        np.savetxt(u_fid, u_cur[:, :], fmt='%10.6f', delimiter=' ', newline='\n')
        np.savetxt(v_fid, v_cur[:, :], fmt='%10.6f', delimiter=' ', newline='\n')
        print "Output data at time = {}days".format(time/(24.0*3600))

# close
eta_fid.close()
u_fid.close()
v_fid.close()
