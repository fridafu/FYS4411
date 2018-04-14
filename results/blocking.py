import numpy as np
import matplotlib.pyplot as plt
import glob

def block_mean(vec):
    return sum(vec)/len(vec)

def meanAndVariance(vec):
    mean = sum(vec)/len(vec)
    var = sum([i ** 2 for i in vec])/len(vec) - mean*mean
    return mean, var

outfile =  open('blocking_taskc.txt', 'w')
outfile.write("#Filename     Energy     Variance \n")

import os
for filename in os.listdir("/Users/frida/Happyday/FYS4411/build-Project1-Desktop_Qt_5_9_1_clang_64bit-Release/c_dt_50vals"):
    if filename.endswith(".dat"):

        data = [float( line.rstrip('\n')) for line in open(filename)]
        n_blocks = 200
        block_size_min = 100
        block_size_max = len(data)/100
        block_step = int ((block_size_max - block_size_min + 1) / n_blocks)
        mean_vec = []
        var_vec = []
        block_sizes = []
        for i in range(0, n_blocks):
            mean_temp_vec = []
            start_point = 0
            end_point = block_size_min + block_step*i
            block_size = end_point
            block_sizes.append(block_size)

        mean_temp_vec.append(block_mean(data[start_point:end_point]))
        start_point = end_point
        end_point += block_size_min + block_step*i
        mean, var = meanAndVariance(mean_temp_vec)
        mean_vec.append(mean)
        var_vec.append(np.sqrt(  var/(len(data)/float(block_size) - 1.0)     ))
        line = []
        mean, var = meanAndVariance(data)
        line.append('{} {} {}\n'.format(filename, mean, var))
        outfile.writelines(line)

outfile.close()

#numerical variance is shit because it scales with the size of h squared, the number of particles and the number of cycles
