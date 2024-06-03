########################### Calculate Absolute Local G #######################################
import clusterpy
import numpy as np
import shapefile
import pandas as pd

# Function to calculate the local Getis-Ord G* statistic
def calculateGetisG(keyList, dataMean, dataStd, dataDictionary, dataLength):
    """
    This function returns the local G statistic for a given region.
    keyList is the list of keys of neighbors
    dataLength is the total number of input data units
    dataMean = global mean
    dataStd = global standard deviation
    """
    sum_val = 0
    for i in keyList:
        sum_val += dataDictionary[i]
    neighborNumber = len(keyList)
    numerator = sum_val - (dataMean * neighborNumber)
    denominator = dataStd * ((float(dataLength * neighborNumber - (neighborNumber ** 2)) / (dataLength - 1)) ** 0.5)
    G = numerator / denominator
    return G

# Load data
# Import shapefile data for the regions
Location = clusterpy.importArcData(".../Mongolian/Mongolia_district_attribute")
wt = Location.Wqueen

# Create new weight matrix that includes each region itself
new_wt = {}
list1 = []
list2 = []
for x in wt:
    wtlist = wt[x]
    wtlist.append(x)
    list1.append(wtlist)
    list2.append(x)
for a, b in enumerate(list2):
    b = int(b)
    new_wt[b] = list1[a]

# Load population data and create data dictionaries
shape = shapefile.Reader(".../Mongolian/Mongolia_district_attribute")
data = pd.read_csv(".../Mongolian/livestock1992_2006.csv")
ori_data = pd.read_csv(".../Mongolian/livestock1992_2006.csv")
id = data.ADM2_CODE
pop = data.Summary99
oripop = ori_data.Summary92 # Set the 1st year population data
ori_Dic = {}
dataDictionary = {}
for a, b in enumerate(id):
    b = int(b)
    dataDictionary[b] = pop[a]
    ori_Dic[b] = oripop[a]

areaKeys = dataDictionary.keys()

# Calculate data mean and standard deviation
dataMean = np.mean(np.double(list(dataDictionary.values())))
dataStd = np.std(np.double(list(dataDictionary.values())))

# Get data length
dataLength = len(shape)

# Calculate first year data mean and standard deviation
fst_dataMean = np.mean(np.double(list(ori_Dic.values())))
fst_dataStd = np.std(np.double(list(ori_Dic.values())))

# Calculate Local G* values for each region
clusterGstrValues = {}
resultstr = []
for x in range(dataLength):
    keyList = new_wt[x]
    currentG = calculateGetisG(keyList, fst_dataMean, fst_dataStd, dataDictionary, dataLength)
    resultstr.append(currentG)
    clusterGstrValues[x] = currentG

# Monte Carlo permutation test for significance
plist = []
for x in areaKeys:
    Nlist = list(range(0, dataLength))
    betterClusters = 0
    number = len(new_wt[x][:-1])
    Nlist.pop(x)
    for j in range(999):  # Monte Carlo permutation test
        permKey = np.random.choice(Nlist, number, False)
        permKey = permKey.tolist()
        permKey.append(x)
        randomG = calculateGetisG(permKey, fst_dataMean, fst_dataStd, ori_Dic, dataLength)
        if clusterGstrValues[x] >= 0:
            if clusterGstrValues[x] < randomG:
                betterClusters += 1
        else:
            if clusterGstrValues[x] > randomG:
                betterClusters += 1
    pValue = (betterClusters + 1) / 1000.00
    plist.append(pValue)

# Save Local G* results to a CSV file
df_G = pd.DataFrame()
df_G['G_str'] = resultstr
df_G['P_sim'] = plist
df_G.to_csv(".../Mongolian/Mongolian_99_92.csv", index=False)

########################### Calculate Absolute Local Moran's I #######################################
def calculateMoranI(ikey, keyList, dataMean, dataStd, dataDictionary, dataLength):
    """
    This function returns the local Moran's I statistic for a given region.
    keyList is the list of the keys of i's neighbors
    dataLength is the total number of input data units
    """
    sum_val = 0
    for j in keyList:
        sum_val += np.double(dataDictionary[j] - dataMean)
    neighborNumber = len(keyList)

    numerator = dataLength * (dataDictionary[ikey] - dataMean) * sum_val
    denominator = dataStd ** 2
    denominator *= neighborNumber

    I = np.double(numerator) / np.double(denominator)

    return I

# Load data for Moran's I
ca = clusterpy.importArcData(".../1KM_2012/Rwan_poly_2012")
wt = ca.Wqueen

# Load population data and create data dictionaries
shape = shapefile.Reader(".../1KM_2012/Rwan_poly_2012")
data = pd.read_csv(".../Rwanda/Rwd_2012.csv")
ori_data = pd.read_csv(".../Rwanda/Rwd_2011.csv")
id = data.FID
pop = data.r__2012
oripop = ori_data.r__2011 # Set the 1st year population data
ori_Dic = {}
dataDictionary = {}
for a, b in enumerate(id):
    b = int(b)
    dataDictionary[b] = pop[a]
    ori_Dic[b] = oripop[a]

areaKeys = dataDictionary.keys()

# Calculate data mean and standard deviation
dataMean = np.mean(np.double(list(dataDictionary.values())))
dataStd = np.std(np.double(list(dataDictionary.values())))
fst_dataMean = np.mean(np.double(list(ori_Dic.values())))
fst_dataStd = np.std(np.double(list(ori_Dic.values())))

# Get data length
dataLength = len(shape)

# Calculate Local Moran's I values for each region
MoranValues = {}
result = []
for x in range(dataLength):
    ikey = x
    keylist = wt[x]
    currentI = calculateMoranI(ikey, keylist, fst_dataMean, fst_dataStd, dataDictionary, dataLength)
    result.append(currentI)
    MoranValues[x] = currentI

# Monte Carlo permutation test for significance
plist = []
for x in areaKeys:
    Nlist = list(range(0, dataLength))
    betterClusters = 0
    number = len(wt[x])
    Nlist.pop(x)
    for j in range(999):  # Monte Carlo permutation test
        permKey = np.random.choice(Nlist, number, False)
        permKey = permKey.tolist()
        randomI = calculateMoranI(x, permKey, fst_dataMean, fst_dataStd, ori_Dic, dataLength)
        if MoranValues[x] >= 0:
            if MoranValues[x] < randomI:
                betterClusters += 1
        else:
            if MoranValues[x] > randomI:
                betterClusters += 1
    pValue = (betterClusters + 1) / 1000.00
    plist.append(pValue)

# Categorize regions based on Moran's I and p-values
idx = []
for x in areaKeys:
    if plist[x] <= 0.05 and MoranValues[x] >= 0:
        if dataDictionary[x] < fst_dataMean:
            idx.append('LL')
        else:
            idx.append('HH')
    elif plist[x] <= 0.05 and MoranValues[x] < 0:
        if dataDictionary[x] < fst_dataMean:
            idx.append('LH')
        else:
            idx.append('HL')
    else:
        idx.append('NS')

print("LL count:", idx.count('LL'))
print("HH count:", idx.count('HH'))
print("LH count:", idx.count('LH'))
print("NS count:", idx.count('NS'))

# Save Local Moran's I results to a CSV file
df_I = pd.DataFrame()
df_I['Sq'] = idx
df_I.to_csv(".../I_new12_11.csv", index=False)
