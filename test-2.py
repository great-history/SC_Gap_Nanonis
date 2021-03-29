import os,sys,struct
import binascii
import traceback
import numpy as np
import re

def find_para(str):
    if re.findall("Experiment=", a) != []:
        para_dict["Experiment"] = a.split("=")[1]
    if re.findall("Grid dim", a) != []:
        para_dict["Grid_dim"] = a.split("=")[1]
    if re.findall("Experiment parameters", a) != []:
        list = (a.split("=")[1]).split(";")
    if re.findall("Channels", a) != []:
        list = (a.split("=")[1]).split(";")

para_dict = {"Experiment":"", "Grid_dim":"","Sweep_Start":"","Sweep_End":"", "Experiment_parameters":{}, "Channel":{}}
with open('C:/Users/yg/Desktop/MATLAB/CsV3Sb5/CsV3Sb5-80mV-110-0 Oe-mapping-001.3ds', 'rb') as of:
    N = -1   # 用N记录当前读取的行号(空行并不算在内)
    of.seek(0,2)  # 指针移动到文件末尾
    size = of.tell()  # 此时指针的位置即为文件末尾的位置
    of.seek(0,0)  # 把指针重新移回文件开头
    # 开始读取dat文件中的数据
    read_line = of.readline()
    a = str(read_line, encoding="utf-8")
    
    # print(read_line)
    i = 0
    while 'END:' not in a:
        if of.tell() >= size:  # 讀到文件末尾了
            break
        read_line = of.readline()
        try:
            a = str(read_line, encoding="utf-8")

            a.split("=")
            print(a)
        except Exception as e:
            print(e)
        if read_line.strip() == "":  # 读到的是空行
            continue
        i += 1
    byte_offset = of.tell()
    of.seek(0,0)
    of.seek(byte_offset)
    griddata = np.fromfile(of, dtype='>f4')

    # membuf = of.read()
    # tag, length = struct.unpack("<HL", membuf[:6])
    # # tag, length = struct.unpack("<HL", membuf[7:13])
    # ta = binascii.hexlify(membuf)
    # a = membuf[1:6]
    # print(a)
    # a = membuf[3:4]
    # print(membuf)
    # ss = str(b'\x18\xfd\xab\x06\xe6\x0b\xadA\x07\x99\xad\xe1\x17j\xaeKH\xd4', 'utf8')
    # decimal = int.from_bytes(membuf[7:13], 'big')
    # print(hex(decimal))  # 0x100
    # i = 0
    # while True:
    #     i += 1
    #     print(i)
    #     line = of.readline()
    #     if not line:
    #         break
    #     else:
    #         try:
    #             #             print(line)
    #             #             print(line.decode('utf8'))
    #             line.decode('utf8')
    #             # 为了暴露出错误，最好此处不print
    #         except:
    #             print(str(line))
