# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
######################################################################
def gauss_solve(arr_A,arr_B,arrx,n,width):#高斯消元法，用于求解矩阵M
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
    #第一步，构造矩阵d，d1~dn-1是可以从y的三阶差商计算得来
    mat_D=np.zeros(len(x_0))
    for i in range(1,len(x_0)-1):
        x0_temp=x_0[i-1:i+2]
        y0_temp=y_0[i-1:i+2]
        for j in range(1,3):
            for q in range(2,j-1,-1):
                y0_temp[q] = (y0_temp[q]-y0_temp[q-1]) / (x0_temp[q]-x0_temp[q-j])
        mat_D[i]=6*y0_temp[2]
    #第二步，构造矩阵h存放hi，hi的计算从书上得来，hi的值会贯穿后面的程序
    mat_h=np.zeros(len(x_0))
    for i in range(1,len(x_0)):
        mat_h[i]=x_0[i]-x_0[i-1]
    #第三步，样条插值必须要有边界条件，我选的是第三种边界条件，需要计算d0和dn
    temp=y_0[0:4]
    for k in range(1,4):
        for i in range(3,k-1,-1):
            temp[i] = (temp[i]-temp[i-1]) / (x_0[i]-x_0[i-k])
    d0=-12*mat_h[1]*temp[3]#边界条件d0的计算

    temp=y_0[len(y_0)-4:]
    for k in range(1,4):
        for i in range(3,k-1,-1):
            temp[i] = (temp[i]-temp[i-1]) / (x_0[len(y_0)-4+i]-x_0[len(y_0)-4+i-k])
    dn=12*mat_h[len(mat_h)-1]*temp[3]#边界条件dn的计算
    mat_D[0]=d0
    mat_D[len(mat_D)-1]=dn#至此，矩阵dn计算完成
    #第四步，构造矩阵A。矩阵A是一个带状矩阵，对角线上是常量2，上带是λ0~λn-1，下带是μ1~μn
    mat_A=np.identity(len(x_0))
    mat_p1=np.zeros(len(x_0))
    mat_p2=np.zeros(len(x_0))
    for i in range(len(x_0)-1):
        mat_p1[i] = mat_h[i + 1] / (mat_h[i] + mat_h[i + 1]) #λi
        mat_p2[i] = 1 - mat_p1[i] #μi
    mat_p1[0] = -2
    mat_p2[len(x_0)-2] = -2
    for i in range(len(x_0)):
        mat_A[i][i] = 2
        if i-1>=0:
            mat_A[i][i - 1] = mat_p2[i - 1]
        if i+1<len(x_0):
            mat_A[i][i + 1] = mat_p1[i]
    #第五步，有了构造好的矩阵A，矩阵d，即可调用高斯消元求解矩阵M。A*M=d
    mat_M=np.zeros(len(x_0))
    gauss_solve(mat_A, mat_D, mat_M, len(x_0), 1)
    #第六步 利用计算求得的矩阵M，求解插值y
    y=np.zeros(len(x))
    for i in range(len(x)):
        k=0
        for j in range(1,len(x_0)):#查找x值所属的区间段
            if x[i]<=x_0[j]:
                k=j
                break
        if k>0:
            a1 = (x_0[k] - x[i])
            a2 = (x[i] - x_0[k - 1])
            h = mat_h[k]
            y[i] = ( mat_M[k - 1] * a1 * a1 * a1 / 6
			+ mat_M[k] * a2 * a2 * a2 / 6
			+ (y_0[k - 1] - mat_M[k - 1] * h * h / 6) * a1
			+ (y_0[k] - mat_M[k] * h * h / 6) * a2 ) / h
            #y[i]=-y[i]
        else:
            print("数据超出范围")
            return -1
    return y
######################################################################
#河床的初始数据
x0 = list(range(0,54,2))
y0 = [0,4.01,6.96,7.96,7.97,8.02,9.05,10.13,11.18,12.26,13.28,12.61,10.22,7.9,7.95,8.86,10.8,10.93,11.23,11.3,10.94,10.1,9.54,8.3,7.3,2.5,0.2]

x=np.linspace(0,51.9,500)
y1=linear(x0,y0,x)
y2=Newton(x0,y0,x)
y3=spline(x0,y0,x)

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
