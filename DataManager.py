# python version 3.8

from PyQt5.QtWidgets import *

class HEADER_of_3DS:
    def __init__(self,
                 name=None,
                 poltype=None,
                 thin=None,
                 thout=None,
                 phi=None,
                 ominc=None,
                 eloss=None,
                 gamma_c=None,
                 gamma_f=None,
                 scattering_axis=None,
                 eval_i=None,
                 eval_n=None,
                 trans_op=None,
                 gs_list=None,
                 temperature=None,
                 spectra=None):
        self.name = name if name is not None else ""
        self.poltype = poltype if poltype is not None else ""  # (str,str)
        self.thin = thin  # float
        self.thout = thout if thout is not None else ""  # float
        self.phi = phi  # float
        self.ominc = ominc if ominc is not None else []  # list of float
        self.eloss = eloss if eloss is not None else []  # list of float
        self.gamma_c = gamma_c if gamma_c is not None else []  # list of float
        self.gamma_f = gamma_f if gamma_f is not None else []  # list of float
        self.scattering_axis = scattering_axis if scattering_axis is not None else [[]]  # list of list
        self.eval_i = eval_i if eval_i is not None else []  # list of list
        self.eval_n = eval_n if eval_n is not None else []  # list of list
        self.trans_op = trans_op if trans_op is not None else [[]]  # list of list
        self.gs_list = gs_list if gs_list is not None else []  # list
        self.temperature = temperature if temperature is not None else "" # float
        self.spectra = spectra if spectra is not None else "{}"

class DataManager_3DS:
    def __init__(self):
        self.spectraBasicDataList = {}
        self.currentSpectraBasicData = {}

    def getNameFromSpectraData(self) -> str:
        return ""

    def addSpectraData(self, spectraData) -> bool:
        return False

    def getSpectraDataByName(self, name: str):
        if name in self.spectraBasicDataList.keys():
            return self.spectraBasicDataList[name]
        else:
            return None

def Read3DSinWindow(mainWindow=None) -> HEADER_of_3DS or None:
    fileName, fileType = QFileDialog.getOpenFileName(mainWindow, r'Load json',
                                                     r'D:\Users\yg\PycharmProjects\spectra_data',
                                                     r'json Files(*.json)')
    if fileName == "":
        return None



if __name__ == "__main__":
    pass