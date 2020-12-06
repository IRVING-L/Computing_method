import numpy as np
import matplotlib.pyplot as plt
import time
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
def spline(x_0,y_0,x):#分段三次样条插值
    mat_D=np.zeros(len(x_0))
    for i in range(1,len(x_0)-1):
        x0_temp=x_0[i-1:i+2]
        y0_temp=y_0[i-1:i+2]
        for j in range(1,3):
            for q in range(2,j-1,-1):
                y0_temp[q] = (y0_temp[q]-y0_temp[q-1]) / (x0_temp[q]-x0_temp[q-j])
        mat_D[i]=6*y0_temp[2]

    mat_h=np.zeros(len(x_0))
    for i in range(1,len(x_0)):
        mat_h[i]=x_0[i]-x_0[i-1]

    temp=y_0[0:4]
    for k in range(1,4):
        for i in range(3,k-1,-1):
            temp[i] = (temp[i]-temp[i-1]) / (x_0[i]-x_0[i-k])
    d0=-12*mat_h[1]*temp[3]

    temp=y_0[len(y_0)-4:]
    for k in range(1,4):
        for i in range(3,k-1,-1):
            temp[i] = (temp[i]-temp[i-1]) / (x_0[len(y_0)-4+i]-x_0[len(y_0)-4+i-k])
    dn=12*mat_h[len(mat_h)-1]*temp[3]
    mat_D[0]=d0
    mat_D[len(mat_D)-1]=dn

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

    mat_M=np.zeros(len(x_0))
    gauss_solve(mat_A, mat_D, mat_M, len(x_0), 1)

    y=np.zeros(len(x))
    for i in range(len(x)):
        k=0
        for j in range(1,len(x_0)):
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
###################################
def TR(sy,a,b,error):
    #1
    h1=b-a
    #2
    Tn=h1*(sy[0]+sy[-1])/2
    #3
    flag0=1
    while flag0:
        h2=h1/2 
        S=0
        x=a+h2
        #4
        flag1=1
        while flag1:
            if x>b:
                flag1=0
            else:
                S+=sy[int((x-a)/(b-a)*len(sy))]
                x+=h1
        T2n=(Tn+h1*S)/2
        if abs(T2n-Tn)<error:
            I=T2n
            flag0=0 
        else:
            Tn=T2n
            h1=h2
    return I
###################################
def Romberg(sy,a,b,error,k=20):
    h1=b-a
    T=np.zeros(k)
    S=np.zeros(k)
    C=np.zeros(k)
    R=np.zeros(k)
    Q=0
    T[0]=(sy[0]+sy[-1])*h1/2
    for i in range(1,k):
        h2=h1/2
        sum=0
        x=a+h2
        while x<b:
            sum+=sy[int((x-a)/(b-a)*len(sy))]
            x+=h1
        T[i]=(T[i-1]+h1*sum)/2
        h1=h2
        if i==k-1:
            print('迭代次数过多，求解失败！')
            return -1
        if i-1>=0:
            S[i-1]=T[i]+(T[i]-T[i-1])/3
        if i-2>=0:
            C[i-2]=S[i-1]+(S[i-1]-S[i-2])/15
        if i-3>=0:
            R[i-3]=C[i-2]+(C[i-2]-C[i-3])/63
        if i-4>=0:
            if abs(R[i-3]-R[i-4])<error:
                Q=R[i-3]
                print('两者相减=',abs(R[i-3]-R[i-4]))
                break
    return Q       
######################################################################
x0=np.linspace(0,1,11)
y0=[]
for i in range(11):
    y0.append(np.exp(x0[i]))
beg=time.time()
x=np.linspace(0,1,10000)
s_y=spline(x0,y0,x)
Q=Romberg(s_y,0,1,0.0001)
I=np.exp(1)-1
print('积分准确值=',I,'数值积分的值=',Q,'error=',abs(Q-I))
