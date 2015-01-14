from collections import Counter
import math

# first, read in data, split it by linebreak and remove the last line
with open('iris.data', 'r') as f:
	data_rows = f.read().split("\n")[0:-1]

# transform data so that we have a 5-tuple
data_tuples = [d.split(',') for d in data_rows]

# cast first four tuple values to floats
data_tuples = [(float(d[0]), float(d[1]), float(d[2]), float(d[3]), d[4]) for d in data_tuples]

# print out the first data point to take a look
print data_tuples[0]	

# split data into training and test data (last 9 are test elements)
data_tuple_length = len(data_tuples)
training_data, test_data = data_tuples[0:data_tuple_length-9], data_tuples[data_tuple_length-9:]

# define a distance function
def calculate_euclidian_distance(e1, e2):
	val = (e1[0] - e2[0])**2 \
		+ (e1[1] - e2[1])**2 \
		+ (e1[2] - e2[2])**2 \
		+ (e1[3] - e2[3])**2 
	return math.sqrt(val)
	
# define 
k = 7
	
# use all test data to predict a label for them
for test_element in test_data:
	
	distances = []
	for training_element in training_data:
		distance = calculate_euclidian_distance(test_element, training_element)
		distances.append([distance, test_element[4]])
		
	distances = sorted(distances, key=lambda e: e[0])
	
	class_labels_of_top_k = []
	
	for i in range(k):
		class_labels_of_top_k.append(distances[i][1])
		
	# Sums up the count for the elemts in the list 
	counted = Counter(class_labels_of_top_k)
	# Returns the highest occurring item
	predicted_label = counted.most_common(1)[0][0]  
	
	print "correctly classified:", predicted_label == test_element[4]

 	
