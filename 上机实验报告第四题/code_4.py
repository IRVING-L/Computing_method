import numpy as np
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt
font = FontProperties(fname=r"c:\windows\fonts\simsun.ttc", size=14) 
import time
######################################################################
def gauss_solve(arr_A,arr_B,arrx,n):#列主元高斯消元法，用于求解矩阵M
    for k in range(n-1):
        if arr_A[k+1][k]>arr_A[k][k]:
            uk=k+1
            for j in range(n):
                temp_a = arr_A[uk][j]
                arr_A[uk][j] = arr_A[k][j]
                arr_A[k][j] = temp_a
            temp_b = arr_B[k]
            arr_B[k] = arr_B[uk]
            arr_B[uk] = temp_b
    #print(arr_A)
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
def TR(sy,a,b,error):#自动求步长h的复化梯形数值积分函数
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

def Romberg(sy,a,b,error,k=20):#romberg函数。没有按照书上的步骤写，goto语句太恶心。
    #核心思想是，5个Tn->4个Cn->3个Sn->2个Rn->两个Rn相减与误差error相比较，如果不满足误差要求，再计算一个新Tn->新Cn->新Sn->新Rn->再相减与error比较
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
                #print('两者相减=',abs(R[i-3]-R[i-4]))
                break
    return Q  
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
            y[i] =( mat_M[k - 1] * a1 * a1 * a1 / 6
			+ mat_M[k] * a2 * a2 * a2 / 6
			+ (y_0[k - 1] - mat_M[k - 1] * h * h / 6) * a1
			+ (y_0[k] - mat_M[k] * h * h / 6) * a2 ) / h
    return y
     
######################################################################
x0=[0,1,2,3,4,5,6,7,8,9,10,11,12,13]
x1=[0,1,2,3,4,5,6,7,8,9,10,11,12]
# y0=[0,202074,177540,56644,17872,6617,2514,1100,462,289,127,150,41,47]
# y1=[0,99681,217126,140241,56448,26080,17194,5564,2289,1197,619,1029,210]
y0=[0,202117,379657,436301,454174,460791,463311,464414,464877,465167,465294,465445,465487,465534]
y1=[0,99691,316799,457025,513467,539545,556738,562302,564591,565788,566407,567437,567647]
x=np.linspace(0,13,1300)
xx=np.linspace(0,12,1300)
s_y0=spline(x0,y0,x)
s_y1=spline(x1,y1,xx)
plt.figure(figsize=(10,5), dpi=60)

plt.plot(x,s_y0, 'r',label="Wandering the earth")
plt.plot(xx,s_y1,'g',label="Warwolf 2")
plt.legend(loc='upper left')
plt.xlabel(u"时间/周", fontproperties=font)
plt.ylabel(u"票房数", fontproperties=font)
plt.title(u"两部电影的总票房变化图",fontproperties=font)
plt.show()

#Q=Romberg(s_y0,0,13,1)
#print('数值积分=',Q)
