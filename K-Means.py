import sys,os,random
import re,pprint
import operator

def INITFeatureDictionary(Patient):
    '''Create a Dictionary with List as the Value and Patient SCN as the Key'''
    global FeatureDict
    Feature=[0,0,0,0,0,0,0,0,0,0]
    for Feat in Feature:
        FeatureDict.setdefault(Patient,[]).append(Feat)


def INITClusterDictionary(Centroid):
    '''Clusters Centroid Coordinates'''
    global ClusterCentroids
    CentroidAttribute=[0,0,0,0,0,0,0,0,0,0]
    for cents in CentroidAttribute:
        ClusterCentroids.setdefault(Centroid,[]).append(cents)


    
class Ddict(dict):
    '''Create a 2 Dimensional Dictionary'''
    def __init__(self, default=None):

        self.default = default



    def __getitem__(self, key):

        if not self.has_key(key):

            self[key] = self.default()

        return dict.__getitem__(self, key)




def ReadBreastCancerFile():
    '''Read The Breast Cancer Data from File'''
    try:
        root=sys.argv[1]
        files=open(root,'r').readlines()
        return files 
    except :

        raise
    return None






def CreateFeatureDictonary(AttLength):

    '''Initialising the Feature Dictionary with the SCN as Key and Attributes as List as Values'''
    global FeatureDict
    global PatientLabel
    global AttributeLength
    #AttLength=AttributeLength
    TempHash={}
    files=ReadBreastCancerFile()
    PatientID=""
    try:
        
        for subject in files:
            row=re.split(r" ",subject)
            PatientID=str(row[1].strip())
            if PatientID in PatientLabel:
                PatientLabel[PatientID]=int(row[-1].strip())
            if PatientID not in PatientLabel:
                PatientLabel[PatientID]=int (row[-1].strip())
                
                
            if PatientID in FeatureDict:
                    #print PatientID
                    #temp=int(TempHash[str(PatientID)])
                    #temp+=1
                    #TempHash[PatientID]=int(temp)
                    index=0
                    for r in range(2,AttLength):
                        temp=float(FeatureDict[PatientID][index])
                        temp+=float(row[r].strip())
                        FeatureDict[PatientID][index]=temp
                        index+=1
                        
            if PatientID not in FeatureDict:
                    INITFeatureDictionary(PatientID)
                    #TempHash[PatientID]=1
                    index=0
                    for r in range(2,AttLength):
                        FeatureDict[PatientID][index]=float(row[r].strip())
                        index+=1
            
                    
        NormalizeRepeatedValues(TempHash)                  
                    
    except:
            raise


def NormalizeRepeatedValues(TempHash):
    '''Not in Use'''
    
    global FeatureDict
    
    for ids in TempHash:
        if int(TempHash[ids])>1:
            divisor=int(TempHash[ids])
            for attributes in range(0,9):
                temp=float(FeatureDict[ids][attributes])
                temp=float(temp/divisor)
                FeatureDict[ids][attributes]=float(temp);
            
            
    CreatingClusterSeeds()

def CreatingClusterSeeds():
    '''Creating Clusters Seeds - Centroids'''
    global FeatureDict
    global ClusterIDs

    ClusterIDs={}
    
    try:
         ClusNumber=int(sys.argv[2])   
         while(len(ClusterIDs)<ClusNumber):
             number=random.randint(1,len(FeatureDict)-1)
             pId=str(FeatureDict.keys()[number])
             if pId not in ClusterIDs:
                 ClusterIDs[pId]=1
         InitialiseClusters()         
    except:
        raise

def InitialiseClusters():
   '''Initialise the Global 2 Dim Cluster Dictionary'''

   global ClusterIDs
   global Clusters
   Clusters=Ddict(dict)
   
   if ClusterIDs:
       for clus in ClusterIDs:
           if clus not in Clusters:
               Clusters[clus]["null"]=0
   
   

def CalculatingInitialCluster():
    '''Calculating the First Iteration Clusters'''

    global FeatureDict
    global Clusters
    global ClusterIDs
    for r in FeatureDict:
        NearestDistance={}
        count=len(ClusterIDs)
        for clus in ClusterIDs:
            EuclideanDistance=0.0
            for i in range(len(FeatureDict[r])-1):
                AttributeValue=float(FeatureDict[r][i])
                ClusterAttribute=float(FeatureDict[clus][i])
                EuclideanDistance+=(float(ClusterAttribute-AttributeValue))
            NearestDistance[clus]=abs(EuclideanDistance)
            count-=1
            if count==0:
                NearestDistanceSorted= sorted(NearestDistance.iteritems(), key=operator.itemgetter(1))
                if NearestDistanceSorted[0][0] in Clusters:
                    Clusters[str(NearestDistanceSorted[0][0])][str(r)]=float(NearestDistanceSorted[0][1])


    
            
            
def RecreateClusters(flag):
    '''Create Centroid Points for the new Iteration'''

    global FeatureDict
    global Clusters
    global ClusterCentroids
    global ClusterIDs
    
    NumberClusters=len(Clusters.keys())
    temp=FeatureDict.keys()[1]
    NumberofAttributes = len(FeatureDict[temp])
    lists=[]
    PointsCount=0
    for i in range(NumberofAttributes):
       lists.append(0.0)
    ClusID=1
    for clus in Clusters:
        PointsCount=0
        for points in Clusters[clus]:
            for i in range(NumberofAttributes):
              if points!='null':
                attvalue=float(FeatureDict[points][i])
                temp=float(lists[i])
                temp+=attvalue
                lists[i]=temp
                
            PointsCount+=1
        PointsCount-=1
        for i in range(NumberofAttributes):
           try:
                lists[i]/=PointsCount
           except:
               pass
                
        
        if ClusID<=NumberClusters:
               INITClusterDictionary(str(ClusID))
               for i in range(NumberofAttributes):
                   ClusterCentroids[str(ClusID)][i]=lists[i]
               ClusID+=1
               
        
def SecondIterationInitialise(number):
     '''Initialise the 2dim Cluster for second iteration'''
     global Clusters
     global ClusterIDs
     Clusters=Ddict(dict)
     ClusterIDs=[]
     for i in range(1,number+1):
         if i not in Clusters:
             Clusters[str(i)]["null"]=0
             ClusterIDs.append(str(i))

def PositivePredictiveValue():
    '''Calculate the PPV value'''
    
    global PatientLabel
    global Clusters
    Labels=Ddict(dict)
    
    for key1 in Clusters:
      for key2 in Clusters[key1]:
        if key2 in PatientLabel:
           if key1 in Labels:     
               if int (PatientLabel[key2]) in Labels[key1]:
                   temp=int(Labels[key1][PatientLabel[key2]])
                   temp+=1
                   Labels[key1][PatientLabel[key2]]=temp
               if int(PatientLabel[key2]) not in Labels[key1]:
                   Labels[key1][PatientLabel[key2]]=1
           if key1 not in Labels:
              Labels[key1][PatientLabel[key2]]=1
    
    #print Labels.items(),'\n'
    for clus in Labels:
         maximum,LabelID=FindMax(Labels,clus)
         Falsepos=FalsePositive(Labels,clus,LabelID)
         PPV=float(float(maximum)/float(maximum+Falsepos))
         print "LabelID: ",LabelID,"True Positive: ",maximum,"False Positives: ",Falsepos,"  PPV: ",PPV      
    

def FalsePositive(Labels,clus,LabelID):

     Fpositives=0
     try:
         for key1 in Labels:
          if key1!=clus:
                if LabelID in Labels[key1]:
                      Fpositives+=int(Labels[key1][LabelID])
     except:
         raise

     return Fpositives                       
              
                               

def FindMax(Labels,Ind):
   maximum=-1
   LabelID=0                             
   if Ind in Labels:
       for item in Labels[Ind]:
             if int(Labels[Ind][item])>maximum:

                    maximum=int(Labels[Ind][item])
                    LabelID=item
   return maximum,LabelID                             
    
        
    
def ClusterIterations():
    '''Loop for the Iterations after the First Iterations'''
    global ClusterIDs

    global ClusterCentroids
    global PrevCluster
    global CurrentCluster

    PrevCluster={}
    CurrentCluster={}
    EndFlag=False
    #PositivePredictiveValue()
    RecreateClusters(True)
    numClus=len(ClusterIDs.keys())
    Iteration=1

    while EndFlag!=True:
        PrevCluster=ClusterRatios()
        SecondIterationInitialise(numClus)
        SubsequentCalculations()
        CurrentCluster=ClusterRatios()
        ClusterCentroids={}
        RecreateClusters(False)
        count=numClus
        Iteration+=1
        for prev in PrevCluster:
          for curr in CurrentCluster:
              if PrevCluster[prev]==CurrentCluster[curr]:
                    if count>0:
                        
                        count-=1
                    if count==0:
                        PositivePredictiveValue()

                        print "\nIterations: ",Iteration 
                        EndFlag=True
                        


def ClusterRatios():
    '''Check for Cluster ratio for Convergence'''
    global Clusters
    global PrevCluster
    global CurrentCluster
    count=0
    Dicts={}
    
    for key1 in Clusters:
        for key2 in Clusters[key1]:
            count+=1
        Dicts[key1]=count    
        count=0
    return Dicts   

def SubsequentCalculations():
    '''Calculations of the Min Distance after First Iteration'''
    global ClusterIDs
    global ClusterCentroids
    global FeatureDict
    global Clusters
    global AttributeLength
    for r in FeatureDict:
        NearestDistance={}
        count=len(ClusterIDs)
        for clus in Clusters:
            EuclideanDistance=0.0
            for i in range(len(FeatureDict[r])-1):
                AttributeValue=float(FeatureDict[r][i])
                ClusterAttribute=float(ClusterCentroids[str(clus)][i])
                EuclideanDistance+=(float(ClusterAttribute-AttributeValue))
                 
            NearestDistance[str(clus)]=abs(EuclideanDistance)
            count-=1
            if count==0:
                NearestDistanceSorted= sorted(NearestDistance.iteritems(), key=operator.itemgetter(1))
                if str(NearestDistanceSorted[0][0]) in Clusters:
                    Clusters[str(NearestDistanceSorted[0][0])][str(r)]=float(NearestDistanceSorted[0][1])
    for c in Clusters:
        if len(Clusters[c])<2:
                 print "Restart the Process\n"
                 CreateFeatureDictonary(AttributeLength)
                 CalculatingInitialCluster()
                 ClusterIterations()
def main():

   global FeatureDict
   global ClusterIDs
   global Clusters
   global ClusterCentroids
   global PatientLabel
   global AttributeLength
   PatientLabel={}
   FeatureDict={}
   ClusterCentroids={}
   Clusters=Ddict(dict)
   AttributeLength=9

   while AttributeLength>=1:
       print "\n--------------New Set of Attributes---------------"
       
       if AttributeLength<3:
           AttributeLength=2
       
       CreateFeatureDictonary(AttributeLength)
       CalculatingInitialCluster()
       print "Attributes: ",AttributeLength
       ClusterIterations()
       AttributeLength-=2

if __name__=="__main__":
    main()

    
    


 

