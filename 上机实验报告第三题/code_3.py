import math
import numpy as np
import matplotlib.pyplot as plt
###################################################################
def sgn(x):#自定义一个sgn判断值正负的函数
    if x>0:
        return 1
    if x<0:
        return -1
    if x==0:
        return 0

def polyfit(x0,y0,p):#x0,y0是初始数据点，p表示拟合多项式的次数
    #第一步：构造矩阵G
    m=len(x0)
    n=p+1
    G=np.zeros((m,n+1))
    for i in range(m):
        for j in range(n):
            G[i][j]=pow(x0[i],j)
    G[0:m,n]=y0
    #第二步：构造矩阵Qk
    for k in range(0,n):
        a1=0
        for i in range(k,m):
            a1+=pow(G[i][k],2)
        a1=-sgn(G[k][k])*math.sqrt(a1)
        w=np.zeros(m)
        w[k]=G[k][k]-a1
        for j in range(k+1,m):
            w[j]=G[j][k]
        B=a1*w[k]
    #第三步：变换Gk-1到Gk
        G[k][k]=a1
        for j in range(k+1,n+1):
            t=0
            for i in range(k,m):
                t+=w[i]*G[i][j]
            t=t/B
            for i in range(k,m):
                G[i][j]+=t*w[i]
    #第四步：解三角矩阵Ra=h
    a=np.zeros(n)
    a[n-1]=G[n-1][n]/G[n-1][n-1]
    for i in range(n-2,-1,-1):
        temp=0
        for j in range(i+1,n):
            temp+=G[i][j]*a[j]
        a[i]=(G[i][n]-temp)/G[i][i]
    print('拟合多项式的各项系数为：')
    for i in range(len(a)):
        print('a',i,'= ',a[i],sep='')
    
    p=0
    for i in range(n+1,m):
        p+=pow(G[i][n],2)
    print('拟合多项式的误差E2=',p)
    
    x=np.linspace(x0[0],x0[-1],len(x0)*20)
    y=0
    for i in range(len(a)):
        y+=pow(x,i)*a[i]
    plt.plot(x0,y0,'*')
    plt.hold
    plt.plot(x,y)
    plt.show()
###################################################################
x0=list(range(1,32))
y0=[34,35,33,34,35,35,34,24,28,31,32,31,36,37,36,29,
    25,26,29,32,35,27,34,34,36,37,38,40,33,34,36]
polyfit(x0,y0,10)
    