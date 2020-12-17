import numpy as np
import time
import math
##########################################################
def judge(a,b):#该函数的作用：用于判断后面for循环里终点的位置，让for循环只处理矩阵里的非0数据
    if a>b:
        return int(b) 
    if a<b:
        return int(a) 
    if a==b:
        return int(a)
#############################
def exchange(atemp,a,rows,width):#该函数的作用是将读取的压缩矩阵（m*n）转变成（m*m）的方阵
    for i in range(rows):
        a[i][i]=atemp[i][width]
        for j in range(1,width+1):
            if atemp[i][width+j]!=0:
                a[i][i+j]=atemp[i][width+j]
            if atemp[i][width-j]!=0:
                a[i][i-j]=atemp[i][width-j]
#############################
def file_read():
    #如果要计算data20194.dat这个52100阶的矩阵，需要修改电脑的配置，开辟更大的虚拟内存（10G以上），否则会报错
    filepath='D:\\QQFiles\\2020-8班大作业\\data20194.dat' #文件的保存路径
    fileinfo=np.fromfile(filepath,dtype=np.int32,count=5) #读取矩阵的headinfo和fileinfo
    flag=fileinfo[1]   #flag是用于判断矩阵是否为压缩矩阵的标志
    rows=cols=fileinfo[2] #矩阵的阶数
    width=fileinfo[3] #矩阵的带宽
    
    mat_data=np.fromfile(filepath,dtype=np.float32)
    if flag==258:#根据矩阵的压缩判断编号确定读取矩阵的方式
        data_a=mat_data[5:rows*cols+5]
        mat_a=np.zeros((rows,cols))
        for i in range(rows):
            mat_a[i,:]=data_a[i*rows:(i+1)*rows]
        mat_b=mat_data[rows*cols+5:]
        return rows,cols,width,mat_a,mat_b
    if flag==514:
        data_a=mat_data[5:rows*(2*width+1)+5]
        mat_a=np.zeros((rows,cols))
        mat_atemp=np.zeros((rows,(2*width+1)))
        for i in range(rows):
            mat_atemp[i,:]=data_a[i*(2*width+1):(i+1)*(2*width+1)]
        exchange(mat_atemp,mat_a,rows,width)
        mat_b=mat_data[rows*(2*width+1)+5:]
        return rows,cols,width,mat_a,mat_b
    
#############################
def gauss_solve(arr_A,arr_B,n,width):
    arrx=np.zeros(n)
    for k in range(n-1):
        if arr_A[k][k]==0:
            print('主元为0，消元失败')
            return -1
        for i in range(k+1,judge(n,k+width+1)):
            temp=arr_A[i][k] / arr_A[k][k]
            for j in range(k,judge(n,k+2*width+2)):
                arr_A[i][j] = arr_A[i][j] - temp * arr_A[k][j]
            arr_B[i] = arr_B[i] - temp * arr_B[k]#消元过程，因为是带状矩阵，就不用列主元消元了

    arrx[n - 1] = arr_B[n - 1] / arr_A[n - 1][n - 1]
    for k in range(n-1,-1,-1):
        temp = arr_B[k]
        for j in range(k+1,judge(n,k+width+1)):
            temp = temp - arr_A[k][j] * arrx[j]
        arrx[k] = temp / arr_A[k][k]#回代过程结束
    return arrx
##########################################################
start=time.time()
rows,cols,width,mat_a,mat_b=file_read()
mid=time.time()
mat_x=gauss_solve(mat_a,mat_b,rows,width)
for i in range(5):
    print('第',i+1,'个解为：',mat_x[i],sep='')
print('......')
for i in range(5):
    print('第',rows-4+i,'个解为：',mat_x[rows-5+i],sep='')
end=time.time()
print('\n程序总执行时间=',end-start,'s')
print('读取矩阵时间=',mid-start,'s')
print('矩阵计算时间=',end-mid,'s')
