from cmath import inf

class Heap:         #Heap class implemented using array
    def __init__(self,l):
        self.data=l
        self.access=[i[1] for i in l]       #initialised the access list to 0,1,2,...n-1
        t=len(self.data)-1
        while t>=0:
            self.heapdown(t)            #Used fast buildheap by heapdown on the initial list
            t-=1      
    def heapdown(self,u):
        t=u
        while (2*t+3<=len(self.data) and (self.data[2*t+1]<self.data[t] or self.data[2*t+2]<self.data[t])) or (2*t+2==len(self.data) and self.data[2*t+1]<self.data[t]) :
            if 2*t+3<=len(self.data):
                if self.data[2*t+1]<=self.data[2*t+2]:
                    v=2*t+1
                else:
                    v=2*t+2
            else:
                v=2*t+1
            self.data[v],self.data[t]=self.data[t],self.data[v]     #childeren swapped
            self.access[self.data[v][1]],self.access[self.data[t][1]]=self.access[self.data[t][1]],self.access[self.data[v][1]]
            #access list also changed due to swapping of the elements
            t=v
        
    def heapup(self,u):
        t=u
        while self.data[(t-1)//2]>self.data[t] and (t-1)//2>=0:
            v=(t-1)//2
            self.data[(t-1)//2],self.data[t]=self.data[t],self.data[(t-1)//2]       #parent swapped
            self.access[self.data[v][1]],self.access[self.data[t][1]]=self.access[self.data[t][1]],self.access[self.data[v][1]]
            #access list also changed due to swapping of the elements   
            t=v
    def extractmin(self):
        if len(self.data)>1:
            w=self.data[0]          #minimum extracted
            x=self.data[-1]
            self.data.pop()
            self.data[0]=x
            j=self.access[x[1]]
            self.access[w[1]]=j
            self.access[x[1]]=0
            self.heapdown(0)
            self.data.append((float('inf'),w[1]))
            return w            #minimum returned
        else:
            w=self.data[0]
            self.data[0]=(float('inf'),w[1])
            return self.data[0]

def sgn(x):
    if not x:
        return 0
    elif abs(x)==x:
        return -1
    else:
        return 1

def listCollisions(M,x,v,m,T):
    f=float('inf')
    l=[]
    for i in range(len(x)-1):       #storing time between successive particles to collide in the array
        if sgn((x[i+1]-x[i]))*sgn((v[i+1]-v[i])) != -1:
            l.append((f,i))
        else:
            l.append((-(x[i+1]-x[i])/(v[i+1]-v[i]),i))
    h=Heap(l)           #initial heap storing time between successive particles to collide
    p=0
    tarr=[0]*len(x)         #stores the time of last collision for every particle i
    tnew=0                  #stores the current time updated with each collision
    ans=[]
    while p<m and tnew<=T:          #p=collision number
        (tnew,i)=h.extractmin()
        if tnew==float('inf') or tnew>T:            #break condition
            break
        m1=M[i]
        m2=M[i+1]
        v1=v[i]
        v2=v[i+1]
        x[i]=x[i]+v[i]*(tnew-tarr[i])            #updated the position of i, x[i] since  the last collision
        x[i+1]=x[i+1]+v[i+1]*(tnew-tarr[i+1])            ##updated the position of i+1, x[i+1] since  the last collision
        tarr[i]=tnew                #updated the time of last collision for particle i
        tarr[i+1]=tnew           #updated the time of last collision for particle i+1
        v[i+1]=2*m1*v1/(m1+m2)-(m1-m2)*v2/(m1+m2)       #velocities change after collision
        v[i]=(m1-m2)*v1/(m1+m2)+2*m2*v2/(m1+m2)
        if 1<=i<=len(x)-3:          #change for i+1 and i-1
            key2=h.access[i+1]
            x1=x[i-1]+v[i-1]*(tnew-tarr[i-1])
            x2=x[i+2]+v[i+2]*(tnew-tarr[i+2])
            if sgn((x2-x[i+1]))*sgn((v[i+2]-v[i+1])) == 1 or sgn((x2-x[i+1]))*sgn((v[i+2]-v[i+1]))==0: 
                h.data[key2]=(f,i+1)        #changed time to infinite
            else:
                h.data[key2]=(tnew-(x2-x[i+1])/(v[i+2]-v[i+1]),i+1)
            h.heapup(key2)              #updated so should remain heap
            h.heapdown(key2)
            key1=h.access[i-1]
            if sgn((x[i]-x1))*sgn((v[i]-v[i-1])) == 1 or sgn((x[i]-x1))*sgn((v[i]-v[i-1])) == 0:
                h.data[key1]=(f,i-1)            #changed time to infinite
                
            else:
                h.data[key1]=(tnew-(x[i]-x1)/(v[i]-v[i-1]),i-1)
            h.heapup(key1)
            h.heapdown(key1)
        elif i==0:          #change for i+1 only
            key2=h.access[i+1]
            x2=x[i+2]+v[i+2]*(tnew-tarr[i+2])
            if sgn((x2-x[i+1]))*sgn((v[i+2]-v[i+1])) == 1 or sgn((x2-x[i+1]))*sgn((v[i+2]-v[i+1]))==0:
                h.data[key2]=(f,i+1)
            else:
                h.data[key2]=(tnew-(x2-x[i+1])/(v[i+2]-v[i+1]),i+1)
            h.heapup(key2)
            h.heapdown(key2)
        elif i==len(x)-2:       #change for i-1 only
            key1=h.access[i-1]      
            x1=x[i-1]+v[i-1]*(tnew-tarr[i-1])
            if sgn((x[i]-x1))*sgn((v[i]-v[i-1])) == 1 or sgn((x[i]-x1))*sgn((v[i]-v[i-1])) == 0 :
                h.data[key1]=(f,i-1)
            else:
                h.data[key1]=(tnew-(x[i]-x1)/(v[i]-v[i-1]),i-1)
            h.heapup(key1)
            h.heapdown(key1)
        ans.append((tnew,i,x[i]))           #appended the tuple tnew i.e the time for collision, i the particle which collided, x[i] the position of collision
        p+=1
    return ans


# Sample Test Case 1
# M = [1, 2, 8,1,2,3,4,5,6,7,8,9,10,11]
# x = [-20,-10,0, 1, 3,19,20,21,45,70,75,81,97,1002]
# v = [8, -4, -1,2,3,4,5,6,7,8,9,10,11,2]
# m = 100000000000000000000000000000
# T = 100000000000000000000000000000
# print((listCollisions(M,x,v,m,T)))
# [(0.8333333333333334, 0, -13.333333333333332), (3.3333333333333335, 1, -3.333333333333334), (100.55555555555556, 12, 1203.111111111111), (114.38418079096044, 11, 1224.8418079096045), (129.67653613498885, 10, 1242.0888252148998), (148.03671992764546, 9, 1254.2937594211637), (173.39281038691698, 8, 1258.7496727084188), (204.4821595128843, 7, 1247.8929570773057), (239.1544354941823, 6, 1215.7721774709116), (282.1550953654062, 5, 1147.620381461625), (337.54600468547943, 4, 1015.6380140564384), (399.74197088197275, 3, 800.4839417639455), (464.6256968006903, 2, 457.9590301340237), (526.731559162313, 3, 823.3719936418244), (796.8276046769422, 4, 1899.0437411723278), (820.9176610929139, 2, 317.1003977704987), (3229.5419936602257, 5, 10585.615899251803)]
