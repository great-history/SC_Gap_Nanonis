import json
import re
import numpy as np
import matplotlib
from scipy import optimize
from scipy import interpolate as interp
from scipy.integrate import quad
from sympy import *
from scipy.special import legendre,lpmv

k_B = 0.0862  # Boltzman常数,单位:meV/K # delta:meV eff_T1:单位K gamma:meV
Im = 1j  # 虚数单位
class DOS():
    def Extend_Energy(self, bias):
        if len(bias) % 2 == 1:
            N = 2 * len(bias) - 1
        else:
            N = 2 * len(bias)
        extend_E = np.linspace(2 * min(bias), 2 * max(bias), N,
                               endpoint=True)  # 把能量范围扩大,因为理论上DOS做的是[-\infty,+\infty]上的积分
        return extend_E

    def DIRAC_SMEAR(self, extend_E, eff_T1):
        Dirac_smear = np.zeros((len(extend_E), 1))
        for i in range(len(extend_E)):
            if extend_E[i] < 0:
                Dirac_smear[i] = (1 / (k_B * eff_T1)) * (np.exp(extend_E[i] / (k_B * eff_T1))) / \
                                 ((1 + np.exp(extend_E[i] / (k_B * eff_T1))) ** 2)
            else:
                dirac_numerator = np.exp(-extend_E[i] / (k_B * eff_T1))
                dirac_denominator = k_B * eff_T1 * (np.exp(-extend_E[i] / (k_B * eff_T1)) + 1) ** 2
                Dirac_smear[i] = dirac_numerator / dirac_denominator
        return Dirac_smear.flatten()  # 列向量之间不能卷积,只有转换为行向量,虽然他们都是(603,1),但实际上是二维的矩阵,因此要转换为一维

    def S_Wave_DOS(self, extend_E, delta1, gamma1):
        dos_numerator = extend_E + Im * gamma1
        dos_denominator = np.sqrt((extend_E + Im * gamma1) ** 2 - delta1 ** 2)
        DOS = abs((dos_numerator / dos_denominator).real)  # 一维数组
        DOS = np.array(DOS)
        return DOS.flatten()

    def Low_Symmetry_Wave_DOS_theta(self, theta, symmetry, E, delta2, gamma2, delta1=0):
        # delta = delta2 * abs(np.cos(symmetry * theta))
        delta_theta = delta2 * np.cos(symmetry * theta) ** 2 + delta1
        DOS_theta = (E - Im * gamma2) / np.sqrt((E - Im * gamma2) ** 2 - delta_theta ** 2)
        return DOS_theta.real

    def Low_Symmetry_Wave_DOS(self, symmetry, extend_E, delta2, gamma2, delta1=0):
        DOS = np.zeros(len(extend_E), dtype=float)
        for i in range(len(extend_E)):
            En = extend_E[i]
            integ = quad(self.Low_Symmetry_Wave_DOS_theta, 0, pi / symmetry, args=(symmetry, En, delta2, gamma2, delta1))[0]
            DOS[i] = abs(integ) * symmetry * 2
        return DOS.flatten()   # 積分一般不用，因爲發現跟離散化求和差不多

    def dos_by_k_dependence_delta(self, extend_E, delta2, gamma2, delta1=0):
        dos_numerator = extend_E + Im * gamma2
        DOS = np.zeros(dos_numerator.shape, dtype=float) # 一定要加complex，否则会默认是实数
        delta = delta2 + delta1
        i = 0
        try:
            for E_n in extend_E:
                dos_denominator = np.sqrt((E_n + Im * gamma2) ** 2 - delta ** 2)
                DOS[i] = sum(abs((dos_numerator[i] / dos_denominator).real))
                i += 1
            DOS = np.array(DOS)
            return DOS.flatten()
        except Exception as e:
            print(str(e) + ":dos_by_k_dependence_delta计算出错")


class DIDV_one_band():
    def __init__(self):
        self.dos = DOS()

    def Return_DIDV(self, extend_E, DOS, Dirac_smear, bias):
        extend_DIDV = np.convolve(DOS, Dirac_smear, 'same')  # 用same确定出与extend_E长度相同的extend_DIDV
        f = interp.interp1d(extend_E, extend_DIDV, kind="quadratic")
        dIdV = f(bias)
        dIdV = dIdV / np.mean(dIdV[0:10])  # 归一化,主要为了便于实验数据和拟合结果的比较,实验数据也经过类似的归一化
        return dIdV

    def isotropic_s_didv(self, bias, delta1, gamma1, eff_T1):  # to do list:其实不一定是扩大两倍,也可以是更高倍数,精度也越高,运算也更慢
        extend_E = self.dos.Extend_Energy(bias)
        DOS = self.dos.S_Wave_DOS(extend_E, delta1, gamma1) * 2 * np.pi
        Dirac_smear = self.dos.DIRAC_SMEAR(extend_E, eff_T1)
        dIdV = self.Return_DIDV(extend_E, DOS, Dirac_smear, bias)
        return dIdV

    def symmetry_didv(self, input, delta2, gamma2, eff_T2):
        bias = input[0]
        n = input[1][0]
        r = np.linspace(0, 1, 2500)
        theta = np.pi * r / n
        delta_k = delta2 * abs(np.cos(n * theta))
        extend_E = self.dos.Extend_Energy(bias)
        DOS = self.dos.dos_by_k_dependence_delta(extend_E, delta_k, gamma2) * (theta[1] - theta[0]) * 2 * n
        Dirac_smear = self.dos.DIRAC_SMEAR(extend_E, eff_T2)
        dIdV = self.Return_DIDV(extend_E, DOS, Dirac_smear, bias)
        return dIdV

    def s_symmetry_one_band_didv(self, input, x, delta0, delta1, gamma1, eff_T1):  # to do list:其实不一定是扩大两倍,也可以是更高倍数,精度也越高,运算也更慢
        bias = input[0]
        extend_E = self.dos.Extend_Energy(bias)
        n = input[1][0]
        r = np.linspace(0, 1, 2500)
        theta = 2 * np.pi * r / n
        delta = abs(np.cos(n*theta) * delta0 * (1-x) + delta1*x)
        DOS = self.dos.dos_by_k_dependence_delta(extend_E, delta, gamma1) * (theta[1] - theta[0]) * n
        Dirac_smear = self.dos.DIRAC_SMEAR(extend_E, eff_T1)
        dIdV = self.Return_DIDV(extend_E, DOS, Dirac_smear, bias)
        return dIdV

    def symmetry_symmetry_one_band_didv(self, input, x, phi, delta0, delta1, gamma1, eff_T1):  # phi 是指相位差
        bias = input[0]
        extend_E = self.dos.Extend_Energy(bias)
        n1 = input[1][0]
        n2 = input[2][0]
        r = np.linspace(0, 1, 2500)
        theta = 2 * np.pi * r
        delta = abs(np.cos(n1 * theta + phi) * delta0 * x + np.cos(n2 * theta) * delta1 * (1-x))
        DOS = self.dos.dos_by_k_dependence_delta(extend_E, delta, gamma1) * (theta[1] - theta[0])
        Dirac_smear = self.dos.DIRAC_SMEAR(extend_E, eff_T1)
        dIdV = self.Return_DIDV(extend_E, DOS, Dirac_smear, bias)
        return dIdV


class DIDV_two_band():
    def __init__(self):
        self.dos = DOS()
        self.didv_one_band = DIDV_one_band()

    def isotropic_2s_didv(self, bias, x, delta1, gamma1, eff_T1, delta2, gamma2, eff_T2):
        dIdV_1 = self.didv_one_band.isotropic_s_didv(bias=bias, delta1=delta1, gamma1=gamma1, eff_T1=eff_T1)
        dIdV_2 = self.didv_one_band.isotropic_s_didv(bias=bias, delta1=delta2, gamma1=gamma2, eff_T1=eff_T2)
        return dIdV_1 * x + dIdV_2 * (1 - x)

    def symmetry_symmetry_two_band_didv(self, input, x, delta1, delta2, gamma1, gamma2, eff_T1, eff_T2):
        bias = input[0]
        n1 = input[1][0]
        n2 = input[2][0]
        input1 = np.vstack((bias, np.ones(len(bias)) * n1))
        input2 = np.vstack((bias, np.ones(len(bias)) * n2))
        dIdV_1 = self.didv_one_band.symmetry_didv(input=input1,delta2=delta1,gamma2=gamma1,eff_T2=eff_T1)
        dIdV_2 = self.didv_one_band.symmetry_didv(input=input2,delta2=delta2,gamma2=gamma2,eff_T2=eff_T2)
        return dIdV_1*x + dIdV_2*(1-x)

    def s_symmetry_two_band_didv(self, input, x, delta1, delta2, gamma1, gamma2, eff_T1, eff_T2):
        bias = input[0]
        dIdV_1 = self.didv_one_band.isotropic_s_didv(bias=bias, delta1=delta1, gamma1=gamma1, eff_T1=eff_T1)
        dIdV_2 = self.didv_one_band.symmetry_didv(input=input, delta2=delta2, gamma2=gamma2, eff_T2=eff_T2)
        return dIdV_1 * x + dIdV_2 * (1 - x)


class Fitting_Methods():
    def __init__(self):
        self.didv_one_band = DIDV_one_band()
        self.didv_two_band = DIDV_two_band()
        self.one_band_functions = {"isotropic_s":self.didv_one_band.isotropic_s_didv, "symmetry":self.didv_one_band.symmetry_didv,
                                   "s_symmetry":self.didv_one_band.s_symmetry_one_band_didv, "symmetry_symmetry":self.didv_one_band.symmetry_symmetry_one_band_didv}
        self.two_band_functions = {"isotropic_s": self.didv_two_band.isotropic_2s_didv, }
        self.info_two_band = {
            "2s": "双带S波拟合结果:\n 成分1:%f \n delta1 = %f \n gamma1 = %f\n eff_T1 = %f \n delta2 = %f \n gamma2 = %f \n eff_T2 = %f",
            "s+2-two-band": "S波+二度对称拟合结果:\n 成分1:%f \n delta1 = %f \n gamma1 = %f\n eff_T1 = %f \n delta2 = %f \n gamma2 = %f \n eff_T2 = %f",
            "s+4-two-band": "S波+四度对称拟合结果:\n 成分1:%f \n delta1 = %f \n gamma1 = %f\n eff_T1 = %f \n delta2 = %f \n gamma2 = %f \n eff_T2 = %f",
            "s+6-two-band": "S波+六度对称拟合结果:\n 成分1:%f \n delta1 = %f \n gamma1 = %f\n eff_T1 = %f \n delta2 = %f \n gamma2 = %f \n eff_T2 = %f",
            "2+4-two-band": "二度+四度对称拟合结果:\n 成分1:%f \n delta1 = %f \n gamma1 = %f\n eff_T1 = %f \n delta2 = %f \n gamma2 = %f \n eff_T2 = %f",
            "2+6-two-band": "二度+六度对称拟合结果:\n 成分1:%f \n delta1 = %f \n gamma1 = %f\n eff_T1 = %f \n delta2 = %f \n gamma2 = %f \n eff_T2 = %f",
            "4+6-two-band": "四度+六度对称拟合结果:\n 成分1:%f \n delta1 = %f \n gamma1 = %f\n eff_T1 = %f \n delta2 = %f \n gamma2 = %f \n eff_T2 = %f",
            }
        self.info_one_band = {"s": "各项同性S波拟合结果:\n delta1 = %f \n gamma1 = %f \n eff_T1 = %f",
                              "2": "二度对称性波拟合结果:\n delta1 = %f \n gamma1 = %f \n eff_T1 = %f",
                              "4": "四度对称性波拟合结果:\n delta1 = %f \n gamma1 = %f \n eff_T1 = %f",
                              "6": "六度对称性波拟合结果:\n delta1 = %f \n gamma1 = %f \n eff_T1 = %f",
                              "s+2": "s+二度对称单带拟合结果:\n x = %f\n delta0 = %f\n delta1 = %f \n gamma1 = %f \n eff_T1 = %f",
                              "s+4": "s+四度对称单带拟合结果:\n x = %f\n delta0 = %f\n delta1 = %f \n gamma1 = %f \n eff_T1 = %f",
                              "s+6": "s+六度对称单带拟合结果:\n x = %f\n delta0 = %f\n delta1 = %f \n gamma1 = %f \n eff_T1 = %f",
                              "2+4": "二度+六度对称性波拟合结果:\n x = %f\n phi = %f\n delta0 = %f\n delta1 = %f \n gamma1 = %f \n eff_T1 = %f",
                              "2+6": "二度+六度对称性波拟合结果:\n x = %f\n phi = %f\n delta0 = %f\n delta1 = %f \n gamma1 = %f \n eff_T1 = %f",
                              "4+6": "四度+六度对称性波拟合结果:\n x = %f\n phi = %f\n delta0 = %f\n delta1 = %f \n gamma1 = %f \n eff_T1 = %f", }

    def fitting_one_band_methods_manager(self, fit_type:str, kwargs:dict):
        bias = kwargs["bias"]
        didv = kwargs["didv"]
        if "para":
            x_bound = kwargs["x"]
            delta0_bound = kwargs["delta0"]
            delta1_bound = kwargs["delta1"]
            delta2_bound = kwargs["delta2"]
            gamma1_bound = kwargs["gamma1"]
            gamma2_bound = kwargs["gamma2"]
            eff_T1_bound = kwargs["eff_T1"]
            eff_T2_bound = kwargs["eff_T2"]
            phi_bound = kwargs["phi"]

        if "fitting":
            if fit_type == "s":
                delta1, gamma1, eff_T1 = \
                    optimize.curve_fit(self.didv_one_band.isotropic_s_didv, bias, didv,
                                       bounds=([delta1_bound[0], gamma1_bound[0], eff_T1_bound[0]],
                                               [delta1_bound[1], gamma1_bound[1], eff_T1_bound[1]]))[0]
                fit_data = self.didv_one_band.isotropic_s_didv(bias=bias, delta1=delta1, gamma1=gamma1, eff_T1=eff_T1)
                para_dict = {"delta1":delta1, "gamma1":gamma1, "eff_T1":eff_T1}
                return (delta1,gamma1,eff_T1), para_dict, fit_data

            if fit_type == "2" or fit_type == "4" or fit_type == "6":
                input = np.vstack((bias, np.ones(len(bias)) * int(float(fit_type)/2)))
                try:
                    delta1, gamma1, eff_T1 = \
                        optimize.curve_fit(self.didv_one_band.symmetry_didv, input, didv,
                                           bounds=([delta1_bound[0], gamma1_bound[0], eff_T1_bound[0]],
                                                   [delta1_bound[1], gamma1_bound[1], eff_T1_bound[1]]))[0]
                    fitdata = self.didv_one_band.symmetry_didv(input=input, delta2=delta1, gamma2=gamma1, eff_T2=eff_T1)
                    para_dict = {"delta1": delta1, "gamma1": gamma1, "eff_T1": eff_T1}
                    return (delta1, gamma1, eff_T1), para_dict, fitdata
                except Exception as e:
                    print(str(e))

            if fit_type == "s+2" or fit_type == "s+4" or fit_type == "s+6":
                try:
                    import  re  # 一定要在这里重新import re，否则会报错，不知道为什么
                    n = float(re.findall("[1-9]{1}", fit_type)[0])
                    input = np.vstack((bias, np.ones(len(bias)) * (n/2)))
                    try:
                        x, delta0, delta1, gamma1, eff_T1 = \
                            optimize.curve_fit(self.didv_one_band.s_symmetry_one_band_didv, input, didv,
                                               bounds=([x_bound[0], delta0_bound[0], delta1_bound[0], gamma1_bound[0], eff_T1_bound[0]],
                                                       [x_bound[1], delta0_bound[1], delta1_bound[1], gamma1_bound[1], eff_T1_bound[1]]))[0]
                        fit_data = self.didv_one_band.s_symmetry_one_band_didv(input=input, x=x, delta0=delta0, delta1=delta1, gamma1=gamma1, eff_T1=eff_T1)
                        para_dict = {"x":x,"delta0":delta0,"delta1": delta1, "gamma1": gamma1, "eff_T1": eff_T1}
                        return (x, delta0, delta1, gamma1, eff_T1), para_dict, fit_data
                    except Exception as e:
                        print(str(e))
                except Exception as e:
                    print(str(e))

            if fit_type == "2+4" or fit_type == "2+6" or fit_type == "4+6":
                import re
                n1 = float((re.findall("[1-9]{1}", fit_type))[0])
                n2 = float((re.findall("[1-9]{1}", fit_type))[1])
                input = np.vstack((bias, np.ones(len(bias)) * int(n1 / 2), np.ones(len(bias)) * int(n2 / 2)))
                try:
                    x, phi, delta0, delta1, gamma1, eff_T1 = \
                        optimize.curve_fit(self.didv_one_band.symmetry_symmetry_one_band_didv, input, didv,
                                           bounds=([x_bound[0], phi_bound[0]*np.pi, delta0_bound[0], delta1_bound[0], gamma1_bound[0], eff_T1_bound[0]],
                                                   [x_bound[1], phi_bound[1]*np.pi, delta0_bound[1], delta1_bound[1], gamma1_bound[1], eff_T1_bound[1]]))[0]
                    fit_data = self.didv_one_band.symmetry_symmetry_one_band_didv(input=input, x=x, phi=phi, delta0=delta0, delta1=delta1, gamma1=gamma1, eff_T1=eff_T1)
                    para_dict = {"x":x, "phi":phi/np.pi,"delta0":delta0, "delta1": delta1, "gamma1": gamma1, "eff_T1": eff_T1}
                    return (x, phi/np.pi, delta0, delta1, gamma1, eff_T1), para_dict, fit_data
                except Exception as e:
                    print(str(e))

    def fitting_two_band_methods_manager(self, fit_type:str, kwargs:dict):
        bias = kwargs["bias"]
        didv = kwargs["didv"]

        if "fitting":
            if fit_type == "2s":
                try:
                    x, delta1, gamma1, eff_T1, delta2, gamma2, eff_T2 = \
                        optimize.curve_fit(self.didv_two_band.isotropic_2s_didv, bias, didv,
                                           bounds=([kwargs["x"][0], kwargs["delta1"][0], kwargs["gamma1"][0], kwargs["eff_T1"][0], kwargs["delta2"][0], kwargs["gamma2"][0], kwargs["eff_T2"][0]],
                                                   [kwargs["x"][1], kwargs["delta1"][1], kwargs["gamma1"][1], kwargs["eff_T1"][1], kwargs["delta2"][1], kwargs["gamma2"][1], kwargs["eff_T2"][1]]))[0]
                    fit_data = self.didv_two_band.isotropic_2s_didv(bias=bias, x=x, delta1=delta1, gamma1=gamma1, eff_T1=eff_T1, delta2=delta2, gamma2=gamma2, eff_T2=eff_T2)
                    para_dict = {"x":x,"delta1":delta1,"gamma1":gamma1,"eff_T1":eff_T1,"delta2":delta2,"gamma2":gamma2,"eff_T2":eff_T2}
                    return (x, delta1, gamma1, eff_T1, delta2, gamma2, eff_T2), para_dict, fit_data
                except Exception as e:
                    print(str(e))
                return

            if fit_type == "s+2-two-band" or fit_type == "s+4-two-band" or fit_type == "s+6-two-band":
                try:
                    import  re  # 一定要在这里重新import re，否则会报错，不知道为什么
                    n = float(re.findall("[1-9]{1}", fit_type)[0])
                    input = np.vstack((bias, np.ones(len(bias)) * int(n/2)))
                    try:
                        x, delta1, delta2, gamma1,gamma2, eff_T1, eff_T2 = \
                            optimize.curve_fit(self.didv_two_band.s_symmetry_two_band_didv, input, didv,
                                               bounds=([0, kwargs["delta1"][0], kwargs["delta2"][0], kwargs["gamma1"][0], kwargs["gamma2"][0],kwargs["eff_T1"][0],kwargs["eff_T2"][0]],
                                                       [1, kwargs["delta1"][1], kwargs["delta2"][1], kwargs["gamma1"][1], kwargs["gamma2"][1],kwargs["eff_T1"][1],kwargs["eff_T2"][1]]))[0]
                        fit_data = self.didv_two_band.s_symmetry_two_band_didv(input=input, x=x, delta1=delta1, gamma1=gamma1, eff_T1=eff_T1,
                                                                               delta2=delta2, gamma2=gamma2, eff_T2=eff_T2)
                        para_dict = {"x": x, "delta1": delta1, "delta2": delta2, "gamma1": gamma1, "gamma2": gamma2, "eff_T1": eff_T1, "eff_T2": eff_T2}
                        return (x, delta1, delta2, gamma1,gamma2, eff_T1, eff_T2), para_dict, fit_data
                    except Exception as e:
                        print(str(e))
                except Exception as e:
                    print(str(e))

            if fit_type == "2+4-two-band" or fit_type == "2+6-two-band" or fit_type == "4+6-two-band":
                import re
                n1 = float((re.findall("[1-9]{1}", fit_type))[0])
                n2 = float((re.findall("[1-9]{1}", fit_type))[1])
                input = np.vstack((bias, np.ones(len(bias)) * int(n1/2), np.ones(len(bias)) * int(n2/2)))
                try:
                    x, delta1, delta2, gamma1, gamma2, eff_T1, eff_T2 = \
                        optimize.curve_fit(self.didv_two_band.symmetry_symmetry_two_band_didv, input, didv,
                                           bounds=([0, kwargs["delta1"][0], kwargs["delta2"][0], kwargs["gamma1"][0], kwargs["gamma2"][0], kwargs["eff_T1"][0], kwargs["eff_T2"][0]],
                                                   [1, kwargs["delta1"][1], kwargs["delta2"][1], kwargs["gamma1"][1], kwargs["gamma2"][1], kwargs["eff_T1"][1], kwargs["eff_T2"][1]]))[0]
                    fit_data = self.didv_two_band.symmetry_symmetry_two_band_didv(input=input, x=x,
                                                                                  delta1=delta1, gamma1=gamma1, eff_T1=eff_T1,
                                                                                  delta2=delta2, gamma2=gamma2, eff_T2=eff_T2)
                    para_dict = {"x": x, "delta1": delta1, "delta2": delta2, "gamma1": gamma1, "gamma2": gamma2, "eff_T1": eff_T1, "eff_T2": eff_T2}
                    return (x, delta1, delta2, gamma1, gamma2, eff_T1, eff_T2), para_dict, fit_data
                except Exception as e:
                    print(str(e))

    def cal_one_band_with_given_parameters(self,fit_type:str, kwargs:dict):
        bias = kwargs["bias"]
        if fit_type == "s":
            dIdV = self.didv_one_band.isotropic_s_didv(bias=kwargs["bias"], delta1=kwargs["delta1"], gamma1=kwargs["gamma1"], eff_T1=kwargs["eff_T1"])
            return dIdV

        if fit_type == "2" or fit_type == "4" or fit_type == "6":
            input = np.vstack((bias, np.ones(len(bias)) * (float(fit_type) / 2)))
            dIdV = self.didv_one_band.symmetry_didv(input=input, delta2=kwargs["delta1"], gamma2=kwargs["gamma1"], eff_T2=kwargs["eff_T1"])
            return dIdV

        if fit_type == "s+2" or fit_type == "s+4" or fit_type == "s+6":
            import re
            n = float(re.findall("[1-9]{1}", fit_type)[0])
            input = np.vstack((bias, np.ones(len(bias)) * (n/2)))
            dIdV = self.didv_one_band.s_symmetry_one_band_didv(input=input, x=kwargs["x"], delta0=kwargs["delta0"], delta1=kwargs["delta1"], gamma1=kwargs["gamma1"], eff_T1=kwargs["eff_T1"])
            return dIdV

        if fit_type == "2+4" or fit_type == "2+6" or fit_type == "4+6":
            import re
            n1 = float(re.findall("[1-9]{1}", fit_type)[0])
            n2 = float(re.findall("[1-9]{1}", fit_type)[1])
            input = np.vstack((bias, np.ones(len(bias)) * (n1 / 2), np.ones(len(bias)) * (n2 / 2)))
            dIdV = self.didv_one_band.symmetry_symmetry_one_band_didv(input=input, x=kwargs["x"], phi=kwargs["phi"], delta0=kwargs["delta0"],
                                                               delta1=kwargs["delta1"], gamma1=kwargs["gamma1"], eff_T1=kwargs["eff_T1"])
            return dIdV

    def cal_two_band_with_given_parameters(self,fit_type:str, kwargs:dict):
        bias = kwargs["bias"]
        if fit_type == "2s":
            dIdV = self.didv_two_band.isotropic_2s_didv(bias=kwargs["bias"], x=kwargs["x"],
                                                        delta1=kwargs["delta1"], gamma1=kwargs["gamma1"], eff_T1=kwargs["eff_T1"],
                                                        delta2=kwargs["delta2"], gamma2=kwargs["gamma2"], eff_T2=kwargs["eff_T2"])
            return dIdV

        if fit_type == "s+2-two-band" or fit_type == "s+4-two-band" or fit_type == "s+6-two-band":
            import re
            n = float(re.findall("[1-9]{1}", fit_type)[0])
            input = np.vstack((bias, np.ones(len(bias)) * int(n / 2)))
            dIdV = self.didv_two_band.s_symmetry_two_band_didv(input=input, x=kwargs["x"],
                                                               delta1=kwargs["delta1"], gamma1=kwargs["gamma1"], eff_T1=kwargs["eff_T1"],
                                                               delta2=kwargs["delta2"], gamma2=kwargs["gamma2"], eff_T2=kwargs["eff_T2"])
            return dIdV

        if fit_type == "2+4-two-band" or fit_type == "2+6-two-band" or fit_type == "4+6-two-band":
            import re
            n1 = float(re.findall("[1-9]{1}", fit_type)[0])
            n2 = float(re.findall("[1-9]{1}", fit_type)[1])
            input = np.vstack((bias, np.ones(len(bias)) * int(n1 / 2), np.ones(len(bias)) * int(n2 / 2)))
            dIdV = self.didv_two_band.symmetry_symmetry_two_band_didv(input=input, x=kwargs["x"],
                                                                      delta1=kwargs["delta1"], gamma1=kwargs["gamma1"], eff_T1=kwargs["eff_T1"],
                                                                      delta2=kwargs["delta2"], gamma2=kwargs["gamma2"], eff_T2=kwargs["eff_T2"])
            return dIdV

    def generate_one_band_gap_func(self,fit_type:str, kwargs:dict):
        r = np.linspace(0, 1, 2500)
        theta = np.pi * r * 2
        if fit_type == "s":
            gap = np.ones(len(theta)) * kwargs["delta1"]
        elif fit_type == "2" or fit_type == "4" or fit_type == "6":
            gap = kwargs["delta1"] * abs(np.cos(float(fit_type) / 2 * theta))
        elif fit_type == "s+2" or fit_type == "s+4" or fit_type == "s+6":
            import re
            n = float(re.findall("[1-9]{1}", fit_type)[0])
            gap = abs(np.cos(n/2 * theta) * kwargs["delta0"] * (1 - kwargs["x"]) + kwargs["delta1"] * kwargs["x"])
        elif fit_type == "2+4" or fit_type == "2+6" or fit_type == "4+6":
            import re
            n1 = float(re.findall("[1-9]{1}", fit_type)[0])
            n2 = float(re.findall("[1-9]{1}", fit_type)[1])
            gap = abs(np.cos(n1/2 * theta + kwargs["phi"]) * kwargs["delta0"] * kwargs["x"] + np.cos(n2/2 * theta) * kwargs["delta1"] * (1 - kwargs["x"]))
        else:
            return False
        return theta, gap

    def generate_two_band_gap_func(self,fit_type:str, kwargs:dict):
        r = np.linspace(0, 1, 2500)
        theta = np.pi * r * 2
        if fit_type == "2s":
            gap1 = np.ones(len(theta)) * kwargs["delta1"]
            gap2 = np.ones(len(theta)) * kwargs["delta2"]

        elif fit_type == "s+2-two-band" or fit_type == "s+4-two-band" or fit_type == "s+6-two-band":
            import re
            n = float(re.findall("[1-9]{1}", fit_type)[0])
            gap1 = np.ones(len(theta)) * kwargs["delta1"]
            gap2 = abs(np.cos(n/2 * theta)) * kwargs["delta2"]

        elif fit_type == "2+4-two-band" or fit_type == "2+6-two-band" or fit_type == "4+6-two-band":
            import re
            n1 = float(re.findall("[1-9]{1}", fit_type)[0])
            n2 = float(re.findall("[1-9]{1}", fit_type)[1])
            gap1 = abs(np.cos(n1 / 2 * theta)) * kwargs["delta1"]
            gap2 = abs(np.cos(n2 / 2 * theta)) * kwargs["delta2"]
        else:
            return False
        return theta, gap1, gap2

