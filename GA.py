import math
import numpy.random as random
import numpy as np
import math
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
def Chorosome(lower1,uper1,lower2,uper2,i,n):#i精确度，染色体个数
    long1=int(math.log2((uper1-lower1)*pow(10,i)))+1#计算染色体长度
    long2=int(math.log2((uper2-lower2)*pow(10,i)))+1
    long=long1+long2
    chorosome=[[str(random.randint(0,2)) for i in range(long)] for j in range(n)]
    Chorosome=[]
    for j in chorosome:
        b=''.join(j)
        Chorosome.append(b)
    return Chorosome,long1
def trans(a,lower,uper,):#a为需转化的基因,lower为下界，uper为上界，两种基因的长度
    s=[]
    for j in a:
        i=int(j,base=2)
        s.append(round((lower)+(i*(uper-lower))/(pow(2,len(j))),5))#解码公式
    return s
def andlist(a,b):#合并列表
    L=[[a[i],b[i]] for i in range(len(a))]
    return L
def fiting(p):#计算适应度，p表现型
    value=[]
    for a,b in p:
        value1=a*math.sin(math.pi*4*a)
        value2=b*math.sin(math.pi*20*b)
        value.append(value1+value2+21.5)
    return value#list_fiting是染色体与适应度的字典
def chosing(value,population):#value适应度且返回为表现型,population种群
    average_P=[i/(sum(value)) for i in value]#计算平均概率
    saverage_P=sorted(average_P,reverse=False)#排序后的平均概率，升序
    list_averageP=andlist(population,average_P)#平均概率与染色体的列表
    sort_list_averageP=sorted(list_averageP,key=lambda x:x[1],reverse=False)#根据平均概率对字典排序，且获得对应的染色体，升序（列表）
    sum_P=[ sum(saverage_P[:i+1]) for i in range(len(average_P))]#根据排序后的平均概率计算累加概率
    newp=[]#新种群
#设计轮盘赌(选择）
    for i in range(len(population)):
        point=random.random()
        if  point<=sum_P[0]:
            newp.append(sort_list_averageP[0][0])
        else:
            for j in range(len(population)):
                try:
                    if point>=sum_P[j]and point<=sum_P[j+1]:
                        newp.append(sort_list_averageP[j+1][0])
                except:
                     continue
    return newp
def crossover(pc,JX):#pc为交叉概率，JX为所需交叉的基因组，（类型为列表）(一条染色体为字符串）
    for i in range(int(len(JX)/2)):
        if pc>=random.random():
            c1,c2=random.randint(1,len(JX),size=2)
            site=random.randint(1,len(JX[c1]))
            JX[c1]=JX[c1][:site]+JX[c2][site:]
            JX[c2]=JX[c2][:site]+JX[c1][site:]
            return JX
    return JX
def Mutation(pm,JX):#pm为变异概率概率，JX为所需变异的基因组，（类型为列表）(一条染色体为字符串）
    if pm>=random.random():
        Chorosome(-3.1,12,4.8,5.1,5,10)
        return Chorosome
    else:
        return JX#这个JX是没有发生变异的
def diving(Chorosome,long1):#将整条染色体分开成两个
    a1=[]
    a2=[]
    for i in range(len(Chorosome)):
        a1.append(Chorosome[i][:long1])
        a2.append(Chorosome[i][long1:])
    return a1,a2
Chorosome,long1=Chorosome(-3.1,12,4.8,5.1,5,10)
a1,a2=diving(Chorosome,long1)
L=[]#存放每代中的最优解
for i in range(100):
    p=andlist(trans(a1,-3.1,12),trans(a2,4.8,5.1))
    value=fiting(p)
    JX=chosing(value,Chorosome)
    JX=crossover(0.6,JX)
    JX=Mutation(0.01,JX)
    a1,a2=diving(JX,long1)
    L.append(max(value))
#画图
x,y=np.mgrid[-3.1:12:20j,4.1:5.8:20j]
z=x*np.sin(math.pi*4*x)+y*np.sin(math.pi*20*y)+21.5
ax=plt.subplot(211,projection='3d')
ax.plot_surface(x,y,z,rstride=1,cstride=1,cmap=plt.cm.coolwarm,alpha=0.8)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.subplot(212)
x=[i for i in range(100)]
y=L
plt.plot(x,y)
plt.show()
