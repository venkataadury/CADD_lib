#!/usr/bin/python3
from rdkit import Chem,DataStructs
import matplotlib.pyplot as plt
import sys
from math import sqrt,cos,sin,acos,asin,exp,log,pi

sqr=lambda x: x*x
distance=lambda p1,p2: sqrt(sqr(p2[0]-p1[0])+sqr(p2[1]-p1[1]))
def printMatrix(mat):
    for v in mat:
        for el in v:
            print(round(el,4),end=", ")
        print()
def getplace(p1,p2,r1,r2):
    d=distance(p1,p2)
    ang=acos((sqr(d)+sqr(r1)-sqr(r2))/(2*d*r1))
    print("Using angle: ",ang*(180.0/pi))
    v=(p2[0]-p1[0],p2[1]-p1[1]); mag=distance((0,0),v)/r1; v=[el/mag for el in v]
    v1=(cos(ang)*v[0]-sin(ang)*v[1]+p1[0],sin(ang)*v[0]+cos(ang)*v[1]+p1[1])
    v2=(cos(-ang)*v[0]-sin(-ang)*v[1]+p1[0],sin(-ang)*v[0]+cos(-ang)*v[1]+p1[1])
    #print(distance(v1,p1))
    #print(distance(v2,p1))
    return (v1,v2)
'''def reducedim(bitvec):
    locs,n=0,0
    for i in range(len(bitvec)):
        if bitvec[i]:
            locs+=i
            n+=1
    return (locs/len(bitvec),n/len(bitvec))''' #Close
def reducedim(bitvec):
    locs,n=0,0
    for i in range(len(bitvec)):
        if bitvec[i]:
            locs+=i
            n+=1
    return (n/len(bitvec),locs/(n*len(bitvec)))
def mydistancemetric(d1,d2):
    nv=0
    for i in range(len(d1)): nv+=d1[i]*d2[i]
    return sqrt(sqr(distance(reducedim(d1),reducedim(d2)))+(nv/(sum(d1)+sum(d2))))
tani=lambda x: 1-x
def place(FDM):
    xs=[0,tani(FDM[0][1])]; ys=[0,0]
    p=getplace((xs[0],ys[0]),(xs[1],ys[1]),tani(FDM[0][2]),tani(FDM[1][2]))[0]
    xs.append(p[0]); ys.append(p[1])
    for i in range(3,len(FDM)):
        d2=tani(FDM[i-2][i])
        d1=tani(FDM[i-1][i])
        print("Placing ",i,"at",d1,"from",(xs[-1],ys[-1]),"and at",d2,"from",(xs[-2],ys[-2]),end="")
        p=getplace((xs[-1],ys[-1]),(xs[-2],ys[-2]),d1,d2)
        print(p)
        e1,e2=0,0
        for t in range(0,len(xs)-2):
            e1+=sqr((distance(p[0],(xs[t],ys[t])))-tani(FDM[t][i]))
            e2+=sqr((distance(p[1],(xs[t],ys[t])))-tani(FDM[t][i]))
        e1=sqrt(e1); e2=sqrt(e2)
        print(e1,e2)
        if e1<e2:
            xs.append(p[0][0])
            ys.append(p[0][1])
            print("Choice 0")
        else:
            xs.append(p[1][0])
            ys.append(p[1][1])
            print("Choice 1")
    return (xs,ys)

def clusterInternal(clusters,FDM,n,step=0.01,r=0.999):
    if(len(clusters)<=n): return clusters
    merge=[]
    for ind in range(len(clusters)):
        el=clusters[ind]
        for it in el:
            id=it[0];
            for ind2 in range(ind+1,len(clusters)):
                el2=clusters[ind2]
                for it2 in el2:
                    id2=it2[0]
                    if FDM[id][id2]>r:
                        merge.append((ind,ind2))
                        break
                if len(merge): break
            if len(merge): break
        if len(merge): break
    if not len(merge):
        r-=step
        print("Distance updated to:",r)
    switches=[i for i in range(len(clusters))];
    for i in range(len(clusters)):
        for mr in merge:
            if switches[mr[0]]==i or switches[mr[1]]==i:
                print("Merge on: ",mr[0],mr[1])
                clusters[mr[0]]+=clusters[mr[1]]
                switches[mr[1]]=mr[0]
                clusters[mr[1]]=[]
    clusters=[l for l in clusters if l]
    return clusterInternal(clusters,FDM,n,step,r)
def cluster(mols,n,distf,plotdistf=None):
    if plotdistf is None: plotdistf=distf;
    data=[Chem.RDKFingerprint(mol,useHs=False) for mol in mols]
    if n<=1: return data;
    FDM=[]
    pFDM=[]
    for i in range(len(data)):
        ci=[1]*i; #ci.append(1)
        di=[1]*i;
        for j in range(i,len(data)):
            ci.append(distf(data[i],data[j]))
            di.append(plotdistf(data[i],data[j]))
        FDM.append(ci)
        pFDM.append(di)
    for i in range(len(data)):
        for j in range(0,i):
            FDM[i][j]=FDM[j][i];
            pFDM[i][j]=pFDM[j][i];
    print("Distance matrix computed")
    printMatrix(pFDM)
    points=place(pFDM)
    nFDM=[]
    for i in range(len(points[0])):
        vi=[]
        for j in range(len(points[0])):
            vi.append(distance((points[0][i],points[1][i]),(points[0][j],points[1][j])))
        nFDM.append(vi);
    print("\n\nnFDM is here")
    printMatrix(nFDM)
    plt.scatter(points[0],points[1],marker='*')
    for i in range(len(mols)):
        plt.text(points[0][i],points[1][i],str(i)+": "+Chem.MolToSmiles(mols[i]),fontsize=9)
        #plt.text(points[0][i],points[1][i],str(i)+": "+str(FDM[i]),fontsize=8)
    #plt.text(points[0],points[1],[Chem.MolToSmiles(v) for v in mols],fontsize=9)
    plt.show()

    clusters=[[(i,data[i])] for i in range(len(data))]
    #while len(clusters)>n:
    return clusterInternal(clusters,FDM,n)

fn=sys.argv[1]
fl=open(fn,"r");
mols=[]
fnames=[]
for f in fl:
    f=f[:len(f)-1]; #Remove the '\n' in the end
    fnames.append(f)
    print(f)
    mols.append(Chem.MolFromPDBFile(f,proximityBonding=False,removeHs=True,sanitize=True))

fails=sum([1 for i in mols if i is None])
print(len(mols),"molecules loaded with",fails,"failures")
print("Cleaning ... ",end="")
mols=[mol for mol in mols if mol]
print("done")
print("Clustering: ")
clst=cluster(mols,int(sys.argv[2]),DataStructs.TanimotoSimilarity,mydistancemetric)
print("Clustered into ",len(clst),"clusters")
K=0
for el in clst:
    print(K)
    for rep in el:
        #print(rep[0]) #ID
        print(rep[0],Chem.MolToSmiles(mols[rep[0]]))
    K+=1
