# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import math
######################################################################
def gauss_solve(arr_A,arr_B,arrx,n):#高斯消元法，用于求解矩阵M
    for k in range(n-1):
        if arr_A[k][k]==0:
            print('主元为0，消元失败')
            return -1
        for i in range(k+1,n):
            temp=arr_A[i][k] / arr_A[k][k]
            for j in range(k,n):
                arr_A[i][j] = arr_A[i][j] - temp * arr_A[k][j]
            arr_B[i] = arr_B[i] - temp * arr_B[k]

    arrx[n - 1] = arr_B[n - 1] / arr_A[n - 1][n - 1]
    for k in range(n-1,-1,-1):
        temp = arr_B[k]
        for j in range(k+1,n):
            temp = temp - arr_A[k][j] * arrx[j]
        arrx[k] = temp / arr_A[k][k]
###################################
def linear(x_0,y_0,x):#分段一次插值
    y=np.zeros(len(x))
    for j in range(len(x)):
        k=0
        for i in range(1,len(x_0)):#查找x值所属的区间段
            if x[j]<=x_0[i]:
                k=i
                break
        if k>0:
             y[j]=-((y_0[k] + (y_0[k] - y_0[k - 1]) * (x[j] - x_0[k]) / (x_0[k] - x_0[k - 1])))
        else:
            print("数据超出范围")
            return -1
    return y
###################################
def Newton(x_0,y_0,x):#分段二次牛顿插值
    y=np.zeros(len(x))
    for n in range(len(x)):
        k=0
        for i in range(1,len(x_0),2):#查找x值所属的区间段
            if x_0[i-1]<=x[n] and x[n]<x_0[i+1]:
                k=i
                break
    
        if k>0:
            x0=x_0[k-1:k+2]
            y0=y_0[k-1:k+2]
    
            for j in range(1,3):
                for i in range(2,j-1,-1):
                    y0[i] = ((y0[i] - y0[i-1]) / (x0[i] - x0[i-j]))

            y[n]=y0[2]
            for i in range(1,-1,-1):
                y[n] = y[n] * (x[n] - x0[i]) + y0[i]
            y[n]=-y[n]
        else:
            print("数据超出范围")
            return -1
    return y
###################################
def spline(x_0,y_0,x):#分段三次样条插值
    mat_D=np.zeros(len(x_0))
    for i in range(1,len(x_0)-1):
        x0_temp=x_0[i-1:i+2]
        y0_temp=y_0[i-1:i+2]
        for j in range(1,3):
            for q in range(2,j-1,-1):
                y0_temp[q] = (y0_temp[q]-y0_temp[q-1]) / (x0_temp[q]-x0_temp[q-j])
        mat_D[i]=6*y0_temp[2]
    mat_D[0]=0
    mat_D[-1]=0
    
    mat_h=np.zeros(len(x_0))
    for i in range(1,len(x_0)):
        mat_h[i]=x_0[i]-x_0[i-1]
    
    mat_A=np.identity(len(x_0))
    mat_p1=np.zeros(len(x_0))
    mat_p2=np.zeros(len(x_0))
    
    for i in range(1,len(x_0)-1):
        mat_p1[i] = mat_h[i + 1] / (mat_h[i] + mat_h[i + 1]) #λi
        mat_p2[i] = 1 - mat_p1[i] #μi
    
    for i in range(len(x_0)):
        mat_A[i][i] = 2
        if i-1>=0:
            mat_A[i][i - 1] = mat_p2[i]
        if i+1<len(x_0):
            mat_A[i][i + 1] = mat_p1[i]
    mat_M=np.zeros(len(x_0))
    gauss_solve(mat_A, mat_D, mat_M, len(x_0))

   
    y=np.zeros(len(x))
    sum=0
    dx=x[1]-x[0]
    for i in range(len(x)):
        k=0
        for j in range(1,len(x_0)):
            k=j
            if x[i]<=x_0[j]:
                break
        if k>0:
            a1 = (x_0[k] - x[i])
            a2 = (x[i] - x_0[k - 1])
            h = mat_h[k]
            y[i] =-( mat_M[k - 1] * a1 * a1 * a1 / 6
			+ mat_M[k] * a2 * a2 * a2 / 6
			+ (y_0[k - 1] - mat_M[k - 1] * h * h / 6) * a1
			+ (y_0[k] - mat_M[k] * h * h / 6) * a2 ) / h
            ###计算曲线长度
            df=(-mat_M[k - 1]*a1*a1/2+mat_M[k]*a2*a2/2+(y_0[k]-y_0[k-1])-(mat_M[k]-mat_M[k - 1])/6*h*h)/h
            sum+=math.sqrt(1+df*df)*dx
    print('三次样条曲线长度：',sum,'m')
            
    return y
######################################################################
#河床的初始数据
x0 = list(range(0,54,2))
y0 = [0,4.01,6.96,7.96,7.97,8.02,9.05,10.13,11.18,12.26,13.28,12.61,10.22,7.9,7.95,8.86,10.8,10.93,11.23,11.3,10.94,10.1,9.54,8.3,7.3,2.5,0.2]

x=np.linspace(0,51.9,500)
y1=linear(x0,y0,x)#分段一次
y2=Newton(x0,y0,x)#分段二次牛顿
y3=spline(x0,y0,x)#样条三次
#绘图
plt.subplot(311)
plt.plot(x,y1)
plt.title('sectional Linear interpolation')
plt.ylabel('depth(m)')

plt.subplot(312)
plt.plot(x,y2)
plt.title('sectional Newton interpolation')
plt.ylabel('depth(m)')

plt.subplot(313)
plt.plot(x,y3)
plt.title('sectional Spline interpolation')
plt.xlabel('width(m)')
plt.ylabel('depth(m)')
plt.show()
