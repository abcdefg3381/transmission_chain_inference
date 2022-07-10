from scipy import stats
import xlrd
import datetime
import pandas as pd
import numpy as np
import os, sys
import random
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gamma
from statistics import mean

import os
if os.name == "nt": # Windows Systems
    plt.rcParams['font.sans-serif'] = ['SimHei']  
else: # MacOS Systems
    plt.rcParams["font.family"] = 'Arial Unicode MS'  
plt.rcParams['axes.unicode_minus'] = False  
os.chdir(os.path.dirname(sys.argv[0]))

# Read Excel data

def readExcel(file_path, sheet_index):
    # open the excel file
    table = xlrd.open_workbook(file_path)
    # choose the sheet to open
    data = table.sheet_by_name(table.sheet_names()[sheet_index])
    rows = data.nrows
    cols = data.ncols
    headers = data.row_values(0)

    table_matrix = []
    for i in range(rows):
        row_data = []
        for j in range(cols):
            cell_value = data.cell(i, j).value
            if cell_value == "":
                cell_value = "NAN"
            row_data.append(cell_value)
        table_matrix.append(row_data)
    return table_matrix, headers

# Get the set of all cluster id
def getClusterIdSet(table_data):
    cluster_set = []
    for data in table_data:
        cluster_set.append(data[3])
    return set(cluster_set)

# Get the set of all patient id
def getPatientSet(table_data):
    patient_set = []
    for data in table_data:
        patient_set.append(data[0])
    return patient_set

# classify the dictionary of each cluster’s data according to cluster id
def classifiedByClusterId(cluster_set, table_data):
    cluster_data_dict = {}
    for cluster_id in cluster_set:
        for i in range(len(table_data)):
            if table_data[i][3] == cluster_id:
                if cluster_id not in cluster_data_dict.keys():
                    cluster_data_dict[cluster_id]=[table_data[i]]                  
                else:
                    cluster_data_dict[cluster_id].append(table_data[i])
    return cluster_data_dict

# classify the dictionary of each patient‘s data according to patient id
def classifiedByPatient(patient_set, table_data):
    patient_data_dict = {}
    for patient_id in patient_set:
        for i in range(len(table_data)):
            if table_data[i][0] == patient_id:
                if patient_id not in patient_data_dict.keys():
                    patient_data_dict[patient_id]=[table_data[i]]                  
                else:
                    patient_data_dict[patient_id].append(table_data[i])
    return patient_data_dict

# Convert date str to datetime
def strToDate(str):
    date = datetime.datetime.strptime(str, '%m/%d/%Y').date()
    return date

# Sorting patients by exposure possible date
def sortPatient(group_data):
    sortdata = sorted(group_data, key = lambda group_data: strToDate(group_data[4]))
    return sortdata

# Calculate the exposure time difference for each cluster's patients two by two
def calculateTimeDifference(group_data, cluster_id):
    data_records = []
    for i in range(len(group_data)):
        for j in range(len(group_data)):
            if i == j or group_data[i][4] == "NA" or group_data[j][4] == "NA":
                continue
            exposure_time_difference = (strToDate(group_data[i][4]) - strToDate(group_data[j][4])).days
            if exposure_time_difference < 0:
                continue
            data_records.append([group_data[i][0], group_data[i][2], group_data[j][0], cluster_id, group_data[i][1], exposure_time_difference])
    return data_records

# Define the segmentation function for human viral load based on the peak (4.31, 8.1) and the slope of decline -0.168
def viralRoad(x, x1, x2):
    if x <= x1:
        return  0
    elif x1 < x < x2:
        return 8.1/4.31 * x 
    else:
        return -0.168 * x + 0.168 * 4.31 + 8.1

# Set transmission risk based on time difference and virus shedding vector
def setTransmissionRisk(data, virus_shedding_vector):
    # There are 3 choices for virus shedding vector
    # 1. Use cynomolgus macaques virus shedding as the virus shedding  vector
    if virus_shedding_vector == 1:
        transmission_rate ={1:0.00065, 2:0.00095, 3:0.00105, 4:0.00125, 5:0.00120, 6:0.00115, 7:0.00105, 8:0.00095, 9:0.00075, 10:0.00050, 11:0.00045, 12:0.00035, 13:0.00030, 14:0}
        for i in range(len(data)):
            if data[i][5] in transmission_rate:
                data[i].append(transmission_rate[data[i][5]])
            else:
                data[i].append(0)
    # 2. Use human virus load time series as the virus shedding  vector
    if virus_shedding_vector == 2:
        for i in range(len(data)):
            y = viralRoad(data[i][5], 0, 4.31)
            data[i].append(y)
    # 3. Use the Gamma probability density function as the virus shedding  vector
    if virus_shedding_vector == 3:        
        for i in range(len(data)):
            if data[i][5] == 0:
                data[i].append(0)
            else:
                data[i].append(gamma.pdf(data[i][5], a=3, loc = 0, scale = 2))
    return data
    
# Calculate the correct rate of matching with manually inferred transmission pair
# Accuracy = number of successful matches / number of known transmission pairs
def calculateAccuarcy(data):
    t, p =0, 0
    for i in range(len(data)): 
        if data[i][2] == 'NA':
            continue
        # p is used to record the number of known infectors
        p+=1
        # When inferred infector = known manually inferred infector:
        if data[i][5] == data[i][2]:
            t+=1
    accuracy = t/p
    return accuracy

# Calculate R 
# R = average number of people an infectee can transmit to = number of infectee/number of infector
def calculateR(data):
    infectee_set = []
    infector_set = []
    for i in range(len(data)): 
        if data[i][5] == 'NA':
            continue
        infectee_set.append(data[i][0])
        infector_set.append(data[i][5])
    R =len(set(infectee_set))/len(set(infector_set))
    return R

# Get the root in the inferred transmission chain
def getRoot(data):
    root = []
    for i in range(len(data)): 
        # If a patient's possible infector is NA, it means that the patient is the root of the cluster
        if data[i][5] == 'NA':
            tmp = [data[i][3], data[i][0]]
            root.append(tmp)
    return root

# Get the infector in the inferred transmission chain
def getInfector(data):
    infector_set = []
    for i in range(len(data)):
        if data[i][5] == 'NA':
            continue
        if data[i][5] not in infector_set:
            infector_set.append(data[i][5])
    return infector_set

# Search all infectee for each infector
def searchPair(queue, data):
    pair = {}
    for i in range(len(queue)):
        infector = queue[i]
        for j in range(len(data)):
            if data[j][5] == infector: 
                if infector not in pair:
                    pair[infector] = [data[j][0]]
                else:   
                    pair[infector].append(data[j][0])
    return pair

# Calculate the generation of each cluster 
# for each cluster, how many generations are infected, and how many people are infected in each generation
def searchGeneration(root,data):
    generation_dict = {}
    for i in range(len(root)):
        # root means there is the first generation, so generation=1, size=1
        generation = 1
        cluster_root = [root[i][1]]
        cluster_id = root[i][0]      
        generation_dict[cluster_id] = {1:1}   
        for j in range(len(data)):
            if data[j][3] != cluster_id:
                continue
            if data[j][0] == cluster_root[0]:
                pair = searchPair(cluster_root, data)
                if not pair:
                    break
                generation+=1
                infectee_number = 0
                for key in pair.keys():
                    infectee = pair[key]
                    tmp = len(infectee)
                    infectee_number += tmp
                generation_dict[cluster_id][generation] = infectee_number              
                while pair:
                    pair = searchPair(infectee, data)
                    if not pair:
                        break
                    generation+=1
                    infectee_number = 0
                    for key in pair.keys():
                        infectee = pair[key]
                        tmp = len(infectee)
                        infectee_number += tmp
                        generation_dict[cluster_id][generation] = infectee_number   
    return generation_dict

if __name__ == '__main__':
    # Read excel data
    table_matrix, headers = readExcel('./hongkong cluster.xls', 0)
    
    # Get the set of all cluster id
    cluster_id_set = getClusterIdSet(table_matrix[1:])

    # Get the set of all patient id
    patient_id_set = getPatientSet(table_matrix[1:])

    # Get the dictionary of each cluster’s data 
    data_by_cluster = classifiedByClusterId(cluster_id_set, table_matrix[1:])

    # Sorting patients by exposure possible date
    for key in data_by_cluster.keys():
        data_by_cluster[key] = sortPatient(data_by_cluster[key])
    
    '''
    Repeat the inference 100 times
    '''  
    # new list to store the results of 100 inferences
    pair_for100 = []
    for j in range(100):  
        # Calculate time difference for each cluster's patients two by two
        all_time_difference_data = []
        for key in data_by_cluster.keys():
            data_record = calculateTimeDifference(data_by_cluster[key], key)
            all_time_difference_data += data_record

        # Set transmission risk based on time difference and choiced virus shedding vector  
        all_time_difference_data = setTransmissionRisk(all_time_difference_data, 3)
         
        # Special case: 
        # When there are multiple cases infected on day 0, one patient is randomly selected as the primary case and all remaining cases automatically become second cases to that primary case.
        new_data_by_cluster = classifiedByClusterId(cluster_id_set, all_time_difference_data) 
        # Find the special clusters with multiple primary cases and randomly select the root
        for key in data_by_cluster.keys():
            cluster_id = key
            earliest_date = data_by_cluster[key][0][4]
            primary_case = []
            for i in range(len(data_by_cluster[key])):
                exposure_time = data_by_cluster[key][i][4]
                patient = data_by_cluster[key][i][0]
                if exposure_time == earliest_date:
                    primary_case.append(patient)
            if len(primary_case) == 1:
                continue
            root = random.choice(primary_case)
            # Removal of root data in the all_time_difference_data list
            for key in new_data_by_cluster.keys():
                if key != cluster_id:
                    continue
                for i in range(len(new_data_by_cluster[key])):
                    patient = new_data_by_cluster[key][i][0]
                    infector = new_data_by_cluster[key][i][2]
                    if patient not in primary_case:
                        continue
                    if patient == root or infector != root:
                        all_time_difference_data.remove(new_data_by_cluster[key][i])

        # Calculate weights for transmission risk
        # Get the dictionary of each patient’s data 
        data_by_patient = classifiedByPatient(patient_id_set, all_time_difference_data)
        # Calculate the weight of transmission risk for all possible infectors for each patient
        all_probability = []
        for key in data_by_patient.keys():
            header = ["id", "infector", "inferred_infector", "cluster_id", "age", "exposure_time_difference", "transmission_rate"]
            data = pd.DataFrame(data_by_patient[key], columns=header)
            data['%'] = 100 * (data['transmission_rate'] / data['transmission_rate'].sum())
            data = data.fillna(0)
            probability = data.values.tolist()
            all_probability += probability
        
        # Randomly select one infector for each patient 
        probability_by_patient = classifiedByPatient(patient_id_set, all_probability)
        transmission_pair = []
        for key in probability_by_patient.keys():
            flag = False
            # When the patient has only one possible infector, select that infector
            if len(probability_by_patient[key])==1:
                flag = True
                transmission_pair.append(probability_by_patient[key][0])
            else:
                # When the patient has more than one possible infector, select one by comparison with a random number
                s = 0 
                rnd = random.uniform(1,100) 
                for i in range(len(probability_by_patient[key])):
                    probability = probability_by_patient[key][i][7]
                    s += probability         
                    if s > rnd or s == 100:              
                        transmission_pair.append(probability_by_patient[key][i])
                        flag = True
                        break  
                if not flag: # Test for no random selection
                    print(key, s, rnd)   

        # Read the initial excel data again
        table_matrix, headers = readExcel('./hongkong cluster.xls', 0)
        # Add randomly selected infector to the initial data
        for i in range(1, len(table_matrix)):
            for k in range(len(transmission_pair)):
                most_probable_infector =  transmission_pair[k][2]
                if transmission_pair[k][0] == table_matrix[i][0]: 
                    table_matrix[i].append(most_probable_infector)
        # Missing case fill in "NA" ('NA' is the root of the cluster)
        for i in range(1, len(table_matrix)):
            if len(table_matrix[i]) < 6:
                table_matrix[i].append('NA')
        pair_for100.append(table_matrix[1:])

    '''
    Calculate epidemiological indicators based on 100 inference results
    '''
    # Used to sum 100 times accuracy, R, and 90 percentile of individual infector R, respectively
    accuracy, inferR, percentile = 0, 0, 0 
    # Store R calculated for each loop (100 values in total) to calculate the 95% CI of 100 R
    infR_set = [] 
    # Used to sum 100 times the generation (how many generations, how many people per generation) and the depth
    size1, size2, size3, size4, size5, size6, size7, size8, size9, size10, size11, depth = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    num1, num2, num3, num4, num5, num6, num7, num8, num9, num10, num11= 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  
    # Store the individual infector R obtained from each loop to calculate the frequency
    infector_R_set = []
    # Store the list including infectors, age and R from each loop
    age_R_set = [] 
    # Store the dictionary including age groups and R mean from each loop
    age_Rmean_set= {} 

    for j in range(len(pair_for100)): 
        # 1. Calculate the accuracy of each inference
        if headers[2] == 'infector_id':
            acc = calculateAccuarcy(pair_for100[j])
            print('(1)accuracy of each inference:', acc)
            accuracy += acc
        if headers[2] == 'contact_id':
            print("(1)This dataset can't calculate accuracy")

        # 2. Calculate R of each inference
        infR = calculateR(pair_for100[j])
        print('(2)R of each inference:', infR)
        inferR += infR
        infR_set.append(infR) 

        # 3. Get generation
        # get root set
        root = getRoot(pair_for100[j])
        # get generation dict for each cluster（cluster id: generation: number of patient）
        generation_dict = searchGeneration(root,pair_for100[j])
        # Calculate average depth = Total number of generations per cluster / total number of clusters
        s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 # frequency of generation
        n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 # number of patient for each generation
        for cluster_id in generation_dict.keys():
            flag = False
            size = len(generation_dict) # number of cluster
            generation = len(generation_dict[cluster_id]) # cluser’s depth
            if generation == 1:
               flag = True
               s1 += 1
               print(cluster_id) 
            if generation == 2:
               flag = True
               s2 += 1
               n1 += generation_dict[cluster_id][1]
               n2 += generation_dict[cluster_id][2]
            if generation == 3:
                flag = True
                s3 += 1
                n1 += generation_dict[cluster_id][1]
                n2 += generation_dict[cluster_id][2]
                n3 += generation_dict[cluster_id][3]
            if generation == 4:
                flag = True
                s4 += 1
                n1 += generation_dict[cluster_id][1]
                n2 += generation_dict[cluster_id][2]
                n3 += generation_dict[cluster_id][3]
                n4 += generation_dict[cluster_id][4]
            if generation == 5:
                flag = True
                s5 += 1
                n1 += generation_dict[cluster_id][1]
                n2 += generation_dict[cluster_id][2]
                n3 += generation_dict[cluster_id][3]
                n4 += generation_dict[cluster_id][4]
                n5 += generation_dict[cluster_id][5]
            if generation == 6:
                flag = True
                s6 += 1
                n1 += generation_dict[cluster_id][1]
                n2 += generation_dict[cluster_id][2]
                n3 += generation_dict[cluster_id][3]
                n4 += generation_dict[cluster_id][4]
                n5 += generation_dict[cluster_id][5]
                n6 += generation_dict[cluster_id][6]
            if generation == 7:
                flag = True
                s7 += 1
                n1 += generation_dict[cluster_id][1]
                n2 += generation_dict[cluster_id][2]
                n3 += generation_dict[cluster_id][3]
                n4 += generation_dict[cluster_id][4]
                n5 += generation_dict[cluster_id][5]
                n6 += generation_dict[cluster_id][6]
                n7 += generation_dict[cluster_id][7]
            if generation == 8:
                flag = True
                s8 += 1
                n1 += generation_dict[cluster_id][1]
                n2 += generation_dict[cluster_id][2]
                n3 += generation_dict[cluster_id][3]
                n4 += generation_dict[cluster_id][4]
                n5 += generation_dict[cluster_id][5]
                n6 += generation_dict[cluster_id][6]
                n7 += generation_dict[cluster_id][7]
                n8 += generation_dict[cluster_id][8]
            if generation == 9:
                flag = True
                s9 += 1
                n1 += generation_dict[cluster_id][1]
                n2 += generation_dict[cluster_id][2]
                n3 += generation_dict[cluster_id][3]
                n4 += generation_dict[cluster_id][4]
                n5 += generation_dict[cluster_id][5]
                n6 += generation_dict[cluster_id][6]
                n7 += generation_dict[cluster_id][7]
                n8 += generation_dict[cluster_id][8]
                n9 += generation_dict[cluster_id][9]
            if generation == 10:
                flag = True
                s10 += 1
                n1 += generation_dict[cluster_id][1]
                n2 += generation_dict[cluster_id][2]
                n3 += generation_dict[cluster_id][3]
                n4 += generation_dict[cluster_id][4]
                n5 += generation_dict[cluster_id][5]
                n6 += generation_dict[cluster_id][6]
                n7 += generation_dict[cluster_id][7]
                n8 += generation_dict[cluster_id][8]
                n9 += generation_dict[cluster_id][9]
                n10 += generation_dict[cluster_id][10]
            if generation == 11:
                flag = True
                s11 += 1
                n1 += generation_dict[cluster_id][1]
                n2 += generation_dict[cluster_id][2]
                n3 += generation_dict[cluster_id][3]
                n4 += generation_dict[cluster_id][4]
                n5 += generation_dict[cluster_id][5]
                n6 += generation_dict[cluster_id][6]
                n7 += generation_dict[cluster_id][7]
                n8 += generation_dict[cluster_id][8]
                n9 += generation_dict[cluster_id][9]
                n10 += generation_dict[cluster_id][10]
                n11 += generation_dict[cluster_id][11]
            if not flag:
                print(cluster_id,size) 
        dep = (2*s2+3*s3+4*s4+5*s5+6*s6+7*s7+8*s8+9*s9+10*s10 +11*s11)/size
        print('(3)generation of each inference: frequency of generation:',s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,'number of patient:',n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,'number of cluster:',size,'average depth:',dep)
        num1 += n1
        num2 += n2
        num3 += n3
        num4 += n4
        num5 += n5
        num6 += n6
        num7 += n7
        num8 += n8
        num9 += n9
        num10 += n10
        num11 += n11
        size1 += s1  
        size2 += s2
        size3 += s3
        size4 += s4
        size5 += s5
        size6 += s6
        size7 += s7
        size8 += s8
        size9 += s9
        size10 += s10
        size11 += s11
        depth += dep

        # 4.Calculate infector R  
        # The x-axis is the number of secondary cases/person（R), and the y-axis is the frequency
        # Get all infectors in inferred tranmission chain
        infectors = getInfector(pair_for100[j])
        # Get infectee of each infector, key is inferred infector id, value is all infectees
        transmission_pair = searchPair(infectors, pair_for100[j])
        # New list to store all infector‘s R
        r_set = []
        for key in transmission_pair.keys():
            # value‘s lengh is the number of infectee for each infector
            r = len(transmission_pair[key])
            r_set.append(r) 
            # Find the infector with the biggest R
            if r > 10: 
                print(key, len(transmission_pair[key]))
        # Sort and calculate 90th percentile
        sorted_r_set = sorted(r_set)          
        per = np.percentile(sorted_r_set,90)  
        print('(4)infector R 90 percentile:', per)  
        percentile += per 
        # Store the individual R obtained for each inference into list of infector_R_set
        infector_R_set.extend([x for x in r_set])

        # 5. Calculate infector R for each age group
        # According to the dictionary of pairs obtained in the previous step into a list of all infector id and R
        infector_R = []
        for key in transmission_pair.keys():
            r = len(transmission_pair[key])
            tmp = [key, r]
            infector_R .append(tmp)
        # Define age group dictionary
        age_group_dict = {1:'0-4岁',2:'5-9岁',3:'10-14岁',4:'15-19岁',5:'20-24岁',6:'25-29岁',7:'30-34岁',8:'35-39岁',9:'40-44岁',10:'45-49岁',11:'50-54岁',12:'55-59岁',13:'60-64岁',14:'65-69岁',15:'70-74岁',16:'75-79岁',17:'80-84岁',18:'85岁以上'}
        # Add the age information of each infector
        for i in range(len(infector_R)):
            infector = infector_R [i][0]
            for k in range(len(pair_for100[j])):
                age = pair_for100[j][k][1]
                if pair_for100[j][k][0] == infector:
                    infector_R[i].append(age)
        infector_age_R = []
        for i in range(len(infector_R)):
            age = infector_R[i][2]
            if age == 'NA':
                continue
            for key in age_group_dict.keys():
                if key == age:
                    tmp = [infector_R[i][0],infector_R[i][1],infector_R[i][2],age_group_dict[key]]
                    infector_age_R.append(tmp) 
        # Use groupby to get the mean value of each age group‘s R
        age_Rmean = pd.DataFrame(infector_age_R, columns=['inferred_infector', 'R','age_group','age'])  
        age_Rmean = age_Rmean.groupby(['age'])
        age_Rmean = age_Rmean['R'].mean()
        age_Rmean = age_Rmean.to_dict() 
        print('(5)R mean of each age group:', age_Rmean) 
        age_R_set.extend(infector_age_R) 
        for key in age_Rmean.keys():
            if key not in age_Rmean_set:
                age_Rmean_set[key] = [age_Rmean[key]]
            else:
                age_Rmean_set[key].append(age_Rmean[key]) 
        print('---------------------------------------------------------------')
    
    print('!!!100 loops end!!!')    
    print('---------------------------------------------------------------')
    
    # End of 100 loops, calculate the mean value of each index, draw distribution
    # 1. Mean value of inferred correct rate
    Accuracy = accuracy/100
    if Accuracy == 0:      
        print("1.This dataset can't calculate accuracy") 
    else:
        print('1.Accuracy is:', Accuracy)

    # 2. R mean and 95% Confidence interval for 100 R
    InferR = inferR /100 
    R_interval = stats.t.interval(0.95, len(infR_set)-1,loc=np.mean(infR_set),scale=stats.sem(infR_set))
    print('2.R is:',InferR,'95% CI is:',R_interval)

    # 3. generation
    print('3.generation(number of patient/cluster): 1st:', num1/100, size1/100,'2nd:',num2/100, size2/100, '3rd:',num3/100, size3/100,'4th:',num4/100, size4/100,'5th:',num5/100, size5/100,'6th:',num6/100, size6/100,'7th:',num7/100, size7/100,'8th:',num8/100, size8/100,'9th:',num9/100, size9/100,'10th:',num10/100, size10/100, 'average depth is:', depth/100)
    
    # 4. frequency mean for infector R and 90 percentile
    R_frequency = {}
    for key in infector_R_set:
            R_frequency[key] = R_frequency.get(key,0) +1 
    for key in R_frequency:
        R_frequency[key] = R_frequency[key]/100 
    Percentile = percentile/100 
    print('4.infector R frequency mean:',R_frequency, '90 percentile:', Percentile)

    # 5. infector R distribution of age group
    # 1.1. draw boxplot(R mean for each age group)
    age_Rmean_distribution = []
    for key in age_Rmean_set.keys():
        for i in range(len(age_Rmean_set[key])):         
            tmp = [key, age_Rmean_set[key][i]]
            age_Rmean_distribution.append(tmp)
    age_Rmean_distribution = sorted(age_Rmean_distribution, key = lambda age_distribution: age_distribution[0])
    age_Rmean_distribution = pd.DataFrame(age_Rmean_distribution, columns=['age', 'R mean'])
    age_Rmean_distribution.head()
    arm = sns.boxplot(x=age_Rmean_distribution['age'],y=age_Rmean_distribution['R mean'])
    arm.set_title('Hongkong')
    y_ticks = np.arange(1,10,1)
    plt.yticks(y_ticks)
    plt.xticks(rotation=45,fontsize=8)
    plt.show()
    # 1.2. draw violinplot (R for each age group) 
    age_R_set = sorted(age_R_set, key = lambda age_R_set: age_R_set[2])
    age_R_set = pd.DataFrame(age_R_set, columns=['inferred infector', 'secondary cases','age group','age'])
    age_R_set.head()
    ar = sns.violinplot(x=age_R_set['age'], y=age_R_set['secondary cases'])
    ar.set_title('Hongkong')
    y_ticks = np.arange(1,12,1)
    plt.yticks(y_ticks)
    plt.xticks(rotation=45,fontsize=8)
    plt.show()
    print(age_R_set)

    # 6. distribution of cluster size(the number of each cluster' patients)
    patient_number = []  
    for key in data_by_cluster.keys():   
        patient_number.append(len(data_by_cluster[key]))
    # Sort and calculate 90th percentile
    sorted_patient_number = sorted(patient_number)
    cluster_size_percentile = np.percentile(sorted_patient_number,90) 
    # calculate the frequency of each cluster size
    cluster_size_frequency = {}
    for key in sorted_patient_number:
        cluster_size_frequency[key] = cluster_size_frequency.get(key,0) +1
    print('5.cluster size frequency:', cluster_size_frequency, '90 percentile:',cluster_size_percentile)

        

        
       





    



  




   