import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

density = 850


class tank:
    # volume:初始油量,limit:输送油量上限，i, j, k:油箱坐标,x, y, z:油箱尺寸
    def __init__(self,  i, j, k, x, y, z, volume, limit):
        self.coordinate = np.array([i, j, k])
        self.xyz = np.array([x, y, z])
        self.limit = limit
        self.volume = volume*density
        self.v = x*y*z
        # 每个时刻的质量
        self.m_array = np.zeros(7200)
        self.m_array[0] = self.volume
        # 以油箱左下角为坐标原点的质心坐标
        self.centroid = np.zeros((7200, 3))


class plane:
    # m:plane质量
    def __init__(self, angle_array, m):
        self.angle_array = angle_array
        self.m = m
        self.centroid = np.zeros((7200, 3))


tank1 = tank(8.91304348, 1.20652174, 0.61669004, 1.5, 0.9, 0.3, 0.3, 1.1)
tank2 = tank(6.91304348, -1.39347826, 0.21669004, 2.2, 0.8, 1.1, 1.5, 1.8)
tank3 = tank(-1.68695652, 1.20652174, -0.28330996, 2.4, 1.1, 0.9, 2.1, 1.7)
tank4 = tank(3.11304348, 0.60652174, -0.18330996, 1.7, 1.3, 1.2, 1.9, 1.5)
tank5 = tank(-5.28695652, -0.29347826, 0.41669004, 2.4, 1.2, 1, 2.6, 1.6)
tank6 = tank(-2.08695652, -1.49347826, 0.21669004, 2.4, 1, 0.5, 0.8, 1.1)
tank_list = [tank1, tank2, tank3, tank4, tank5, tank6]


gas_supply = np.array(pd.read_excel('附件2-问题1数据.xlsx'))[:, 1:]
plane_angle = pd.read_excel('附件2-问题1数据.xlsx', sheet_name=1)

plane1 = plane(np.array(plane_angle.iloc[:, 1]), 3000)

M_list = []
for i in range(6):
    for j in range(1, 7200):
        tank_list[i].m_array[j] = tank_list[i].m_array[j-1]-gas_supply[j, i]
        if i == 1:
            tank_list[i].m_array[j] += gas_supply[j, 0]
        if i == 4:
            tank_list[i].m_array[j] += gas_supply[j, 5]

for i in range(7200):
    for j in range(6):
        s = tank_list[j].m_array[i]/density/tank_list[j].xyz[1]
        tanthe = np.tan(plane1.angle_array[i]*np.pi/180)
        x = tank_list[j].xyz[0]
        y = tank_list[j].xyz[1]
        z = tank_list[j].xyz[2]

        if tanthe < z/x:
            if s < 0.5*(x**2)*tanthe:
                a = 2*s/(2*s*tanthe)**0.5
                b = (2*s*tanthe)**0.5
                tank_list[j].centroid[i, 0] = a/3
                tank_list[j].centroid[i, 2] = b/3
            elif s > (x*z-0.5*(x**2)*tanthe):
                b = ((x*z-2*s)*2*tanthe)**0.5
                a = x*z/b
                s1 = 0.5*(x-a)*0.5*z
                rx1, ry1 = 0.5*(x-a), 0.5*z
                s2 = s-s1
                rx2, ry2 = x-a*(3*z-b)/(6*z-3*b), (z**2 +
                                                   (z-b)**2+z*(z-b))/(6*z-3*b)
                tank_list[j].centroid[i, 0] = (s1*rx1+s2*rx2)/s
                tank_list[j].centroid[i, 2] = (s2*ry1+s2*ry2)/s
            else:
                a = s/x-0.5*x*tanthe
                b = s/x+0.5*x*tanthe
                tank_list[j].centroid[i, 0] = x-x*(a+2*b)/(3*a+3*b)
                tank_list[j].centroid[i, 2] = (a**2+b**2+a*b)/(3*a+3*b)
        else:
            if s < 0.5*(z**2)/tanthe:
                a = 2*s/(2*s*tanthe)**0.5
                b = (2*s*tanthe)**0.5
                tank_list[j].centroid[i, 0] = a/3
                tank_list[j].centroid[i, 2] = b/3
            elif s > (x*z-0.5*(z**2)/tanthe):
                b = ((x*z-2*s)*2*tanthe)**0.5
                a = x*z/b
                s1 = 0.5*(z-a)*0.5*x
                rx1, ry1 = 0.5*(z-a), 0.5*x
                s2 = s-s1
                rx2, ry2 = z-a*(3*x-b)/(6*x-3*b), (x**2 +
                                                   (x-b)**2+x*(x-b))/(6*x-3*b)
                tank_list[j].centroid[i, 0] = (s1*rx1+s2*rx2)/s
                tank_list[j].centroid[i, 2] = (s2*ry1+s2*ry2)/s
            else:
                a = s/x-0.5*z/tanthe
                b = s/x+0.5*z/tanthe
                tank_list[j].centroid[i, 0] = z-z*(a+2*b)/(3*a+3*b)
                tank_list[j].centroid[i, 2] = (a**2+b**2+a*b)/(3*a+3*b)

        tank_list[j].centroid[i, 0] += tank_list[j].coordinate[0]
        tank_list[j].centroid[i, 0] -= x/2
        tank_list[j].centroid[i, 1] += tank_list[j].coordinate[1]
        tank_list[j].centroid[i, 2] += tank_list[j].coordinate[2]
        tank_list[j].centroid[i, 2] -= z/2

    moment = 0
    M = plane1.m
    for k in range(6):
        moment += tank_list[k].centroid[i, 0]*tank_list[k].m_array[i]
        M += tank_list[k].m_array[i]
    plane1.centroid[i, 0] = moment/M

    M_list.append(M)

    moment = 0
    for k in range(6):
        moment += tank_list[k].centroid[i, 2]*tank_list[k].m_array[i]
    plane1.centroid[i, 2] = moment/M

    moment = 0
    for k in range(6):
        moment += tank_list[k].centroid[i, 1]*tank_list[k].m_array[i]
    plane1.centroid[i, 1] = moment/M
