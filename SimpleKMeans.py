'''
Created on Dec 4, 2012

@author: ralekar
'''
'''
INPUT THE DATA 
CREATE THE FIRST CLUSTER
SECOND DO THE ITERATIONS
MAX WILL GET THE CLUSTER NAME
SPLITTING CRITERIA IF THE NUMBER OF CLUSTERS < REQUIRED NUMBER OF CLUSTERS 
 THEN SPLIT
    check the dispersion 
produce a confidence interval and distribute all the points again to other 
'''


import sys,os,random,math
import re,pprint
import operator

class Ddict(dict):
    '''Create a 2 Dimensional Dictionary'''
    def __init__(self, default=None):

        self.default = default



    def __getitem__(self, key):

        if not self.has_key(key):

            self[key] = self.default()

        return dict.__getitem__(self, key)

def INITFeatureDictionary(dataPoint,feature):
    '''Create a Dictionary with List as the Value and Patient SCN as the Key'''
    global featureDictionary
    for feat in feature:
        featureDictionary.setdefault(dataPoint,[]).append(feat)

def INITClusterDictionary(centroid,centroidList):
    '''Clusters Centroid Coordinates'''
    global centroidPoints
    for cents in centroidList:
        centroidPoints.setdefault(centroid,[]).append(cents)



def breakDataPPV(ftrain):
    
    step=len(ftrain)/10
    if len(ftrain)%10!=0:
        remainder=len(ftrain)%10

    blockPart=[]     
    for testBlock in range(0,len(ftrain),step):
        if testBlock+step+remainder<=len(ftrain):
            blockPart.append(testBlock)
            blockPart.append(testBlock+step)
    blockPart[-1]+=remainder
            
    for i in range(2,len(blockPart),2):
        blockPart[i]+=1            
    return blockPart
            
def readData():
    
    ppvFlag=int(sys.argv[4])
    fread=open(sys.argv[1],"r").readlines()
    if ppvFlag==1:
       blockPart=breakDataPPV(fread)
       
       for block in range(0,len(blockPart),2):
           chunk=fread[blockPart[block]:blockPart[block+1]]
           initDataStructures()
           for row in chunk:
              line=re.split(",",str(row.strip()))
              fillFeatureDictionary(line)
           iterationTillConvergence()   
    else:
        initDataStructures()
        for read in fread:
          line=re.split(",",str(read.strip()))
          fillFeatureDictionary(line)
        iterationTillConvergence()
                    
def fillFeatureDictionary(line):
    global featureDictionary
    global idLabel
    
    if line[0] in featureDictionary:
        print "Not in Feature Dictionary"
    if line[0] not in featureDictionary:
        newList=[]
        for n in line[1:-1]:
            newList.append(float(n))
        INITFeatureDictionary(line[0],newList)            
        idLabel[line[0]]=int(line[-1])
          
def createClusterSeeds(flag):
    
    '''Creating Clusters Seeds - Centroids'''
    global featureDictionary
    global centroidPoints
    global prevCentroidPoints
    try:
       if flag==True:
         K=int(sys.argv[3])   
         for k in range(K):
             number=random.randint(1,len(featureDictionary)-1)
             ID=str(featureDictionary.keys()[number])
             INITClusterDictionary(str(k+1),featureDictionary[ID])
       if flag==False:
           calculateCentroidPoints()
                      
    except:
        raise


def clusterCalculations():
     
     global featureDictionary
     global centroidPoints
     global clusterSet
     global idLabel   
     initClusterPoints()
     for row in featureDictionary:
         rowFeature=featureDictionary[row]
         points={}
         for centroid in centroidPoints:
             list=centroidPoints[centroid]
             eDistance=0.0
             index=0
             for cell in rowFeature:
                 eDistance+=((list[index]-cell)*(list[index]-cell))
                 index+=1
             points[centroid]=math.sqrt(abs(eDistance))    
         pointsSorted= sorted(points.iteritems(), key=operator.itemgetter(1))
         clusterSet[str(pointsSorted[0][0])][str(row)]=idLabel[str(row)]
    

def numberOfFeatures():
    
    featureSize=int(sys.argv[2])
    lists=[]
    for f in range(featureSize):
        lists.append(0)
    return lists    


def initClusterPoints():
    global clusterSet
    global centroidPoints
    K=int(sys.argv[3])
    for k in centroidPoints.keys():
        clusterSet[str(k)]["-1"]=0.0


def calculateCentroidPoints():

    global centroidPoints
    global clusterSet
    global featureDictionary
    global prevCentroidPoints
    prevCentroidPoints=centroidPoints
    centroidPoints={}
    for cluster in clusterSet:
        avgCentroid=numberOfFeatures()
        for points in clusterSet[cluster]:
            if points!="-1":
                list=featureDictionary[points]
                index=0
                for l in list:
                    avgCentroid[index]+=float(l)
                    index+=1
        ind=0        
        for point in avgCentroid:
            temp=float(float(point)/(len(clusterSet[cluster])))
            avgCentroid[ind]=temp
            ind+=1
        INITClusterDictionary(cluster, avgCentroid)
    
              
def iterationTillConvergence():

    global prevCentroidPoints
    global centroidPoints
    mainFlag=False
    iterator=0
    createClusterSeeds(True)
    clusterCalculations()
    while mainFlag==False:
        iterator+=1
        createClusterSeeds(False)
        mainFlag=convergenceCheck(iterator)
        if mainFlag==False:
            clusterCalculations()
              

def convergenceCheck(iterator):

    global prevCentroidPoints
    global centroidPoints
    STOP=50
    threshold=0.000001
    if iterator==STOP:
        return True
    else:
        for point in centroidPoints:
            current=centroidPoints[point]
            previous=prevCentroidPoints[point]
            index=0
            for curr in current:
                if abs((curr-previous[index]))>threshold:
                    return False
                index+=1
    print "Iterations ",iterator
    eliminateNulls()
    distribution()
    maxLabel()        
    return True


def eliminateNulls():
    global clusterSet
    for cluster in clusterSet:
        if '-1' in clusterSet[cluster]:
            del clusterSet[cluster]["-1"]
    
                                
def distribution():
    global clusterSet
    distDict={}
    for point in clusterSet:
        for id in clusterSet[point]:
            val=clusterSet[point][id]
            if val in distDict:
                list=distDict[val]
                value=list[int(point)-1]
                value+=1
                list[int(point)-1]=value
                distDict[val]=list
            if val not in distDict:
                list=[]
                for l in range(0,len(clusterSet)):
                       list.append(0)
                value=list[int(point)-1]
                value+=1
                list[int(point)-1]=value        
                distDict[val]=list
    
    f=open("analyse.json","w")
    f.write(str(distDict))
    f.close()

def maxLabel():
    
    global clusterSet
    for cluster in clusterSet:
        dicts={}
        for point in clusterSet[cluster]:
            if clusterSet[cluster][point] in dicts:
               temp=int(dicts[clusterSet[cluster][point]])
               temp+=1
               dicts[clusterSet[cluster][point]]=temp
            if clusterSet[cluster][point] not in dicts:
                dicts[clusterSet[cluster][point]]=1
        maxLabel= sorted(dicts.iteritems(), key=operator.itemgetter(1))
        if int(sys.argv[5])==1 and int(sys.argv[3])==3:        
            dendogram(maxLabel)
        if int(sys.argv[5])==1 and int(sys.argv[3])!=3:
            print "Hierarchical Clustering works for 3 clusters only"  
            sys.exit(0)                   
        if int(sys.argv[5])!=1 and int(sys.argv[3])!=3:
            calculatePPV(maxLabel) 
        if int(sys.argv[5])!=1 and int(sys.argv[3])==3:
            calculatePPV(maxLabel) 
            

def dendogram(maxLabel):
    
    global clusterSet
    list1_8={1:0,2:0,3:0,4:0,5:0,6:0,7:0,8:0}
    list9_10={9:0,10:0}
    list11_29={11:0,12:0,13:0,14:0,15:0,16:0,17:0,18:0,19:0,20:0,21:0,22:0,23:0,24:0,25:0,26:0,27:0,28:0,29:0}
    dicts={1:0,2:0,3:0}
    for tple in maxLabel:
        if int(tple[0]) in list1_8:
            dicts[1]+=int(tple[1])
        if int(tple[0]) in list9_10:
            dicts[2]+=int(tple[1])
        if int(tple[0]) in list11_29:
            dicts[3]+=int(tple[1])
    maxLabel= sorted(dicts.iteritems(), key=operator.itemgetter(1))         
    calculatePPV(maxLabel)
    
def calculatePPV(maxLabel):
    
    global avgPPV
    print "*********For a New Cluster***********"
    false=0
    if maxLabel:
        true=maxLabel[-1][1]
        for label in range(0,len(maxLabel)-1):
            false+=int(maxLabel[label][1])
        print "Majority Label: ",maxLabel[-1][0],"PPV:" ,float(true)/float(false+true)        
        avgPPV+=(float(true)/float(false+true))
     
        
def initDataStructures():
    
    global featureDictionary
    global idLabel
    global centroidPoints
    global clusterSet
    featureDictionary={}
    centroidPoints={}
    clusterSet=Ddict(dict)
    idLabel={}
    
def main():
    global featureDictionary
    global idLabel
    global centroidPoints
    global clusterSet
    global avgPPV
    avgPPV=0.0
    readData()
    print "Average PPV ",avgPPV/int(sys.argv[3])
    
    
if __name__=="__main__":
    main()
    