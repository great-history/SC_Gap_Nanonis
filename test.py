import numpy as np
import matplotlib.pyplot as plt

# datalist = []
# with open('C:/Users/yg/Desktop/MATLAB/CsV3Sb5/CsV3Sb5-80mV-110-0 Oe-mapping-001.3ds', 'r') as of:
#     N = -1   # 用N记录当前读取的行号(空行并不算在内)
#     of.seek(0,2)  # 指针移动到文件末尾
#     size = of.tell()  # 此时指针的位置即为文件末尾的位置
#     of.seek(0,0)  # 把指针重新移回文件开头
#     # 开始读取dat文件中的数据
#     print(size)
#     read_line = of.readline()
#     read_line = read_line.strip()
#     print(size)
#     while read_line != '[DATA]':
#         if of.tell() >= size:
#             print("End Of File")  # 在GUI中出现信息框,并返回
#             break
#         read_line = of.readline()
#         print(read_line)
#         read_line = read_line.strip()
#         if read_line.strip() == "":  # 读到的是空行
#             continue
#         N += 1
#     # 读取频道
#     label_line = of.readline()  # 字符串
#     label_line = label_line.strip()  # 仍为字符串
    # label_list = label_line.split('	')  # 列表list
    # # print(label_list)
    # AVG_index_list = []
    # for i in range(len(label_list)):
    #     if label_list[i].find("LI Demod 1 X [AVG]") != -1:
    #         AVG_index_list.append(i)
    # # 读取数据
    # data_lines = of.readlines()  # a list of strings
    # raw_data_matrix = []
    # for line in data_lines:
    #     line_new = line.strip('\n')
    #     line_new = line_new.split('\t')
    #     raw_data_matrix.append(line_new) # a list of lists of strings
    # raw_data_matrix = np.array(raw_data_matrix,dtype="float64")  # ndarray
    # dim = raw_data_matrix.shape  # 例如(301,45)
    # # 获取偏压
    # bias_array = np.zeros((1,dim[0]),dtype='float64')
    # bias_array += (raw_data_matrix.T)[0]  # shape例如:(1,301)
    # # 获取原始DIDV数据的平均
    # raw_avg_DIDV_data = np.zeros((1,dim[0]),dtype='float64')  # shape例如:(1,301)
    # for index in AVG_index_list:
    #     raw_avg_DIDV_data = np.add(raw_avg_DIDV_data,(raw_data_matrix.T)[index])
    # raw_avg_DIDV_data = np.true_divide(raw_avg_DIDV_data, np.mean(raw_avg_DIDV_data[0][0:10])) # 高能量处的DIDV进行归一化,为什么不是全部??
    # bias_array = bias_array[0]
    # raw_avg_DIDV_data = raw_avg_DIDV_data[0]



# # nothing valuable
# bias_array_new = bias_array*10 - 0.02
# bias_array_new = bias_array_new[::-1]
# fit_Num = 30
# tot_Num = 301
# bias_fit = np.hstack((bias_array_new[0:fit_Num+1],
#                                   bias_array_new[tot_Num-fit_Num:tot_Num]))
# data_fit = np.hstack((raw_avg_DIDV_data[0:fit_Num+1],
#                       raw_avg_DIDV_data[tot_Num-fit_Num:tot_Num]))
# from scipy import optimize
# def f_lin(x: float, k: float, b: float):
#     return k * x + b
# k,b = optimize.curve_fit(f_lin, xdata=bias_fit,ydata=data_fit,bounds=([-np.inf, -np.inf], [np.inf, np.inf]))[0]

# a = (1,2,3)
# print(a[1])

import re
a = re.findall("[1-9]{1}", "s+2")
print(a)