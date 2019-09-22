
import math

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  

from DubinsAirplaneFunctions import * 
from PlottingTools import plot3
import numpy as np
import time
import sys
from Dubins_test import get_dubin_L
from pconf import * 



if DATA == 1:
    ALPHA1 = 25
    ALPHA2 = 15
    BETA1 = 20
    BETA2 = 25
    THETA = 30
    DATASET = "dataset.csv"
    DELTA = 0.001

else:
    ALPHA1 = 20
    ALPHA2 = 10
    BETA1 = 15
    BETA2 = 20
    THETA = 20
    DATASET = "dataset2.csv"
    DELTA = 0.001


def f1(n):
    return n
def f2(n):
    return n*n
def f3(n):
    return n*n*n
def f4(n):
    return math.exp(n)
def f5(n):
    math.log(n)

funcs = [f1,f2,f3,f4,f5]

def get_rotation_angle(x,y,z):
    his_cos = (x * 0 + y * 0 + z * 1)/((math.sqrt(x*x+y*y+z*z))*(1))
    angle_fake = math.acos(his_cos)

    if x>0 and y<0 and z<0:
        return 2*pi - angle_fake
    elif x<0 and y<0 and z<0:
        return 2*pi - angle_fake
    elif x>0 and y<0 and z>0:
        return 2*pi - angle_fake
    elif x<0 and y<0 and z>0:
        return 2*pi - angle_fake
    
    return angle_fake


def read_dataset_dubin(file):
    f = open(file,"r")
    points = []
    for line in f:
        # pid, x, y, z, ptype, ptype2 = line.split(",")
        points.append(line.split(","))
    graph = []
    graph_line_points = []
    error_cnt = 0
    angles = []
    for i in range(len(points)):
        graph_row = []
        angle_row = []
        # print(i,'/',len(points))
        sys.stdout.write('\r dubin processing :' + str(i) + '/' + str(len(points)))
        sys.stdout.flush()
        for j in range(len(points)):

            if i == j:
                graph_row.append(0)
                angle_row.append([])
                continue
            x1 = float(points[i][1])
            y1 = float(points[i][2])
            z1 = float(points[i][3])

            x2 = float(points[j][1])
            y2 = float(points[j][2])
            z2 = float(points[j][3])

            xn = float(points[-1][1])
            yn = float(points[-1][-2])
            zn = float(points[-1][3])

            qishi_angle = get_rotation_angle(xn-x1,yn-y1,zn-z1)
            end_angle = get_rotation_angle(xn-x2, yn-y2, zn-z2)
            angle_row.append([qishi_angle, end_angle])
            
            try:
                # print("a")
                L12, LSolution = get_dubin_L(x1,y1,z1,qishi_angle,x2,y2,z2,end_angle)

                if GEN_POINTS:
                    path_dubins_airplane = ExtractDubinsAirplanePath( LSolution)

                    np.save("./dubins_routes/d"+str(DATASET)+"/" + str(i)+"to"+str(j)+".npy", path_dubins_airplane.T ) 

            except:
                # print(x1,y1,z1,qishi_angle,x2,y2,z2,end_angle)
                L12 = MAX_value
                error_cnt += 1 


            # graph_row.append(math.sqrt(dx*dx + dy*dy + dz*dz))
            graph_row.append(L12)
        angles.append(angle_row)
        graph.append(graph_row)
    
    pointType = []
    for p in points:
        pointType.append(int(p[4]))
    
    jiaozhengType = []
    for p in points:
        jiaozhengType.append(int(p[5]))
    print(error_cnt)
    

    # np.save("A.npy",A)
    return graph, pointType, points, jiaozhengType, graph_line_points, angles


def read_dataset(file):
    f = open(file,"r")
    points = []
    for line in f:
        # pid, x, y, z, ptype, ptype2 = line.split(",")
        points.append(line.split(","))
    graph = []
    for i in range(len(points)):
        graph_row = []
        for j in range(len(points)):
            x1 = float(points[i][1])
            y1 = float(points[i][2])
            z1 = float(points[i][3])

            x2 = float(points[j][1])
            y2 = float(points[j][2])
            z2 = float(points[j][3])

            dx = x2 - x1
            dy = y2 - y1
            dz = z2 - z1
            graph_row.append(math.sqrt(dx*dx + dy*dy + dz*dz))
        graph.append(graph_row)
    
    pointType = []
    for p in points:
        pointType.append(int(p[4]))
    
    jiaozhengType = []
    for p in points:
        jiaozhengType.append(int(p[5]))

    return graph, pointType, points, jiaozhengType


def dijkstra(graph, pointType, s):
    # 判断图是否为空，如果为空直接退出
    # pointType: -1:起始点,0:水平点,1:垂直点,2:终点
    if graph is None:
        return None
    dist = [MAX_value]*len(graph)
    dist[s] = 0

    shuiping = [MAX_value] * len(graph)
    chuizhi = [MAX_value] * len(graph)

    shuiping[s] = 0
    chuizhi[s] = 0

    qianqu = [MAX_value] * len(graph)
    tujing_num = [MAX_value] * len(graph)
    tujing_num[s] = 0

    S = []
    Q = [i for i in range(len(graph))]
    dist_init = [i for i in graph[s]]
    while Q:
        u_dist = min([d for v, d in enumerate(dist_init) if v in Q])
        u = dist_init.index(u_dist)
 

        S.append(u)
        Q.remove(u)
 
        for v, d in enumerate(graph[u]):
            if 0 < d < MAX_value:
                if dist[v] > dist[u]+d:
                # if ALPHA * dist[v] + BETA * tujing_num[v] > ALPHA * (dist[u] + d) + BETA * (tujing_num[u] + 1):
                    if pointType[v] == 0 and shuiping[u] + d * DELTA < BETA2 \
                        and chuizhi[u] + d * DELTA < BETA1: 
                        dist[v] = dist[u]+d
                        dist_init[v] = dist[v]
                        shuiping[v] = 0
                        chuizhi[v] = chuizhi[u] + d*DELTA
                        qianqu[v] = u 
                        tujing_num[v] = tujing_num[u] + 1

                    if pointType[v] == 1 and shuiping[u] + d * DELTA < ALPHA2 \
                        and chuizhi[u] + d * DELTA < ALPHA1: 
                        dist[v] = dist[u]+d
                        dist_init[v] = dist[v]
                        chuizhi[v] = 0
                        shuiping[v] = shuiping[u] + d * DELTA
                        qianqu[v] = u
                        tujing_num[v] = tujing_num[u] + 1

                    if pointType[v] == 2 and shuiping[u] + d* DELTA < THETA \
                        and chuizhi[u] + d * DELTA < THETA:
                        dist[v] = dist[u] + d
                        dist_init[v] = dist[v]
                        qianqu[v] = u
                        tujing_num[v] = tujing_num[u] + 1
    
    return dist, qianqu, tujing_num


def dijkstra_with_stepnum(graph, pointType, s, max_stepnum, func_index, jiaozhengType):
    # 判断图是否为空，如果为空直接退出
    # pointType: -1:起始点,0:水平点,1:垂直点,2:终点
    if graph is None:
        return None
    dist = [MAX_value]*len(graph)
    dist[s] = 0

    shuiping = [MAX_value] * len(graph)
    chuizhi = [MAX_value] * len(graph)

    shuiping[s] = 0
    chuizhi[s] = 0

    qianqu = [MAX_value] * len(graph)
    tujing_num = [MAX_value] * len(graph)
    tujing_num[s] = 0

    p_safe = [0] * len(graph)
    p_safe[s] = 1

    S = []
    Q = [i for i in range(len(graph))]
    dist_init = [i for i in graph[s]]

    dist_weighted = [MAX_value] * len(graph)
    dist_weighted[s] = 0

    while Q:
        u_dist = min([d for v, d in enumerate(dist_init) if v in Q])
        u = dist_init.index(u_dist)
 

        S.append(u)
        Q.remove(u)
        # print(u)
        sys.stdout.write('\r dijkstra running:' + str(len(Q)) + '/' + str(len(graph)))
        sys.stdout.flush()
        # print(len(Q),'/',len(graph))
        for v, d in enumerate(graph[u]):
            if 0 < d < MAX_value and v in Q:
                # if dist[v] > dist[u]+d:
                # print(u,v)
                route1, route2 = get_two_routes(u,v,qianqu)
                p1,p2 = get_safe_P(route1,0,0,graph,pointType,jiaozhengType),get_safe_P(route2,0,0,graph,pointType,jiaozhengType)
                # print(route1,route2,p1,p2)
                if ALPHA * (dist[v]/graph[0][-1]) + BETA * funcs[FUNC_INDEX](((tujing_num[v] + 0.0)/max_stepnum)) + GAMA * p2 > \
                     ALPHA * ((dist[u] + d)/graph[0][-1]) + BETA * funcs[FUNC_INDEX](((tujing_num[u] + 1.0)/max_stepnum)) + GAMA * p1:
                    if pointType[v] == 0 and shuiping[u] + d * DELTA < BETA2 \
                        and chuizhi[u] + d * DELTA < BETA1:
                        # if jiaozhengType[v] == 1:
                        dist[v] = dist[u]+d
                        # dist_init[v] = dist[v]
                        dist_init[v] = ALPHA * (dist[v]/graph[0][-1]) + BETA * funcs[FUNC_INDEX](((tujing_num[v]+0.0)/max_stepnum)) + GAMA * p2
                        shuiping[v] = 0
                        chuizhi[v] = chuizhi[u] + d*DELTA
                        qianqu[v] = u 
                        tujing_num[v] = tujing_num[u] + 1

                    if pointType[v] == 1 and shuiping[u] + d * DELTA < ALPHA2 \
                        and chuizhi[u] + d * DELTA < ALPHA1: 
                        dist[v] = dist[u]+d
                        # dist_init[v] = dist[v]
                        dist_init[v] = ALPHA * (dist[v]/graph[0][-1]) + BETA * funcs[FUNC_INDEX](((tujing_num[v]+0.0)/max_stepnum)) + GAMA * p2

                        chuizhi[v] = 0
                        shuiping[v] = shuiping[u] + d * DELTA
                        qianqu[v] = u
                        tujing_num[v] = tujing_num[u] + 1

                    if pointType[v] == 2 and shuiping[u] + d* DELTA < THETA \
                        and chuizhi[u] + d * DELTA < THETA:
                        dist[v] = dist[u] + d
                        # dist_init[v] = dist[v]
                        dist_init[v] = ALPHA * (dist[v]/graph[0][-1]) + BETA * funcs[FUNC_INDEX](((tujing_num[v]+0.0)/max_stepnum)) + GAMA * p2

                        qianqu[v] = u
                        tujing_num[v] = tujing_num[u] + 1
    
    return dist, qianqu, tujing_num


def astar(graph, pointsType, s, jiaozhengType, max_time):
    open_list = []
    close_list = []

    parent = [MAX_value] * len(graph)
    parent[s] = 0
    
    G = [MAX_value] * len(graph)
    H = [MAX_value] * len(graph)

    shuiping = [MAX_value] * len(graph)
    chuizhi = [MAX_value] * len(graph)
    shuiping[s] = 0
    chuizhi[s] = 0

    open_list.append(s)
    G[s] = 0
    H[s] = 0
    
    Dist_ = [MAX_value] * len(graph)
    Cishu_ = [MAX_value] * len(graph)

    Dist_[s] = 0
    Cishu_[s] = 0

    while checklist(graph, open_list, close_list):
        
        min_d = open_list[0]
        min_f = G_P *  G[min_d] + H_P * H[min_d]
        for v,d in enumerate(open_list):
            if G_P * G[d] + H_P * H[d] < min_f:
                min_f = G_P * G[d] + H_P * H[d]
                min_d = d
        
        open_list.remove(min_d)
        close_list.append(min_d)
        node_current = min_d

        for i,d in enumerate(graph[node_current]):
            
            if i in close_list:
                continue

            sp = shuiping[node_current] + d * DELTA
            cz = chuizhi[node_current] + d * DELTA
            if pointsType[i] == 0 and sp  <= BETA2 and cz <= BETA1:
                #可达
                current_route = get_route(parent, points, pointsType, node_current)
                current_route.append(i)
                sp_ = get_safe_P(current_route, 0, 0, graph, pointsType, jiaozhengType)
                # print(sp_)
                if i not in open_list:
                    open_list.append(i)
                    parent[i] = node_current
                    shuiping[i] = 0
                    chuizhi[i] = cz
                    Dist_[i] = Dist_[node_current] + d
                    Cishu_[i] = Cishu_[node_current] + 1
                    G[i] = ALPHA * (Dist_[i]/graph[0][-1]) + BETA * funcs[FUNC_INDEX](((Cishu_[i]+0.0)/max_time)) + GAMA * sp_
                    H[i] = graph[i][-1]
                else:
                    G_current = ALPHA * ((Dist_[node_current] + d)/graph[0][-1]) + BETA * funcs[FUNC_INDEX](((Cishu_[node_current] + 1.0)/max_time)) + GAMA * sp_
                    if G_current < G[i]:
                        parent[i] = node_current
                        G[i] = G_current
                        shuiping[i] = 0
                        chuizhi[i] = cz
                        Dist_[i] = Dist_[node_current] + d
                        Cishu_[i] = Cishu_[node_current] + 1
                        H[i] = graph[i][-1]
            elif pointsType[i] == 1 and sp <= ALPHA2 and cz <= ALPHA1:
                #可达
                current_route = get_route(parent, points, pointsType, node_current)
                current_route.append(i)
                sp_ = get_safe_P(current_route, 0, 0, graph, pointsType, jiaozhengType)
                if i not in open_list:
                    open_list.append(i)
                    parent[i] = node_current
                    shuiping[i] = sp
                    chuizhi[i] = 0
                    Dist_[i] = Dist_[node_current] + d
                    Cishu_[i] = Cishu_[node_current] + 1
                    G[i] = ALPHA * (Dist_[i]/graph[0][-1]) + BETA * funcs[FUNC_INDEX](((Cishu_[i]+0.0)/max_time)) + GAMA * sp_
                    H[i] = graph[i][-1]
                else:
                    G_current = ALPHA * ((Dist_[node_current] + d)/graph[0][-1]) + BETA * funcs[FUNC_INDEX](((Cishu_[node_current] + 1.0)/max_time)) + GAMA * sp_
                    if G_current < G[i]:
                        parent[i] = node_current
                        G[i] = G_current
                        shuiping[i] = sp
                        chuizhi[i] = 0
                        Dist_[i] = Dist_[node_current] + d
                        Cishu_[i] = Cishu_[node_current] + 1
                        H[i] = graph[i][-1]
                
            elif pointsType[i] == 2 and sp <= THETA and cz <= THETA:
                #可达
                current_route = get_route(parent, points, pointsType, node_current)
                current_route.append(i)
                sp_ = get_safe_P(current_route, 0, 0, graph, pointsType, jiaozhengType)
                if i not in open_list:
                    open_list.append(i)
                    parent[i] = node_current
                    shuiping[i] = 0
                    chuizhi[i] = 0
                    Dist_[i] = Dist_[node_current] + d
                    Cishu_[i] = Cishu_[node_current] + 1
                    G[i] = ALPHA * (Dist_[i]/graph[0][-1]) + BETA * funcs[FUNC_INDEX](((Cishu_[i]+0.0)/max_time)) + GAMA * sp_
                    H[i] = graph[i][-1]
                else:
                    G_current = ALPHA * ((Dist_[node_current] + d)/graph[0][-1]) + BETA * funcs[FUNC_INDEX](((Cishu_[node_current] + 1.0)/max_time)) + GAMA * sp_
                    if G_current < G[i]:
                        parent[i] = node_current
                        G[i] = G_current
                        shuiping[i] = 0
                        chuizhi[i] = 0
                        Dist_[i] = Dist_[node_current] + d+ GAMA * sp_
                        Cishu_[i] = Cishu_[node_current] + 1    
                        H[i] = graph[i][-1]   
            else:
                #不可达
                continue
    
    return parent, Dist_
            




def checklist(graph, open_list, close_list):
    if len(open_list) == 0:
        return False
    if len(graph) - 1 in open_list:
        return False
    return True


def get_safe_P(routes, chushi_shuiping, chushi_chuizhi, graph, pointType, jiaozhengType):
    if len(routes) <= 1:
        return 1
    
    # if chushi_shuiping + d * DELTA
    d = graph[routes[0]][routes[1]]
    shuiping = chushi_shuiping + d * DELTA
    chuizhi = chushi_chuizhi + d * DELTA
    if pointType[routes[1]] == 0 and shuiping < BETA2 \
        and chuizhi < BETA1:
        if jiaozhengType[routes[1]] == 0: #shuiping
            return 1 * get_safe_P(routes[1:],0.0, chuizhi, graph, pointType, jiaozhengType)
        else:
            return 0.8* get_safe_P(routes[1:],0.0, chuizhi, graph, pointType, jiaozhengType) + 0.2 *\
                 get_safe_P(routes[1:],min(shuiping,5), chuizhi, graph, pointType, jiaozhengType)
    elif pointType[routes[1]] == 1 and shuiping< ALPHA2 and chuizhi < ALPHA1:
        if jiaozhengType[routes[1]] == 0:
            return 1 * get_safe_P(routes[1:],shuiping, 0.0, graph, pointType, jiaozhengType)
        else:
            return 0.8*get_safe_P(routes[1:],shuiping, 0.0, graph, pointType, jiaozhengType) +\
                0.2 * get_safe_P(routes[1:],shuiping, min(chuizhi,5), graph, pointType, jiaozhengType)
    elif pointType[routes[1]] == 2 and shuiping < THETA and chuizhi < THETA:
        return 1
    else:
        return 0


def get_two_routes(u,v,qianqu):
    #routes1: 要比较的
    #routes2: 现有的
    route1 = [v,u]
    route2 = [v]
    
    tmp = qianqu[u]
    while tmp != 0 and tmp != MAX_value:
        route1.append(tmp)
        tmp = qianqu[tmp]
    route1.append(0)
    route1_r = []
    for i in range(len(route1)):
        route1_r.append(route1[-1-i])
    
    tmp = qianqu[v]
    while tmp!=0 and tmp != MAX_value:
        route2.append(tmp)
        tmp = qianqu[tmp]
    route2.append(0)
    route2_r = []
    for i in range(len(route2)):
        route2_r.append(route2[-1-i])
    
    return route1_r, route2_r

def get_route(qianqu, points, pointsType, qidian ):
    s = ""
    
    q = qianqu[qidian]
    qs = []
    if qidian == -1:
        qs.append(len(points) - 1)
    else:
        qs.append(qidian)
    while q!=0:
        qs.append(q)
        # typ = ""
        
        # s += str(q) + "("+ str(pointsType[q]) + ") <- "
        # s += str(q) + "-"
        q = qianqu[q]
    qs.append(0)

    qss = []
    for i in range(len(qs)):
        qss.append(qs[-1-i])
    return qss
    # return qs

def plot_route(distance, qianqu, points, pointsType, tujing_num):
    # print(distance)
    # print("--------")
    # print(qianqu)
    # print("++++++++")
    # print(tujing_num)
    s = ""
    
    q = qianqu[-1]
    qs = []
    qs.append(len(points)-1)
    while q!=0:
        qs.append(q)
        typ = ""
        
        s += str(q) + "("+ str(pointsType[q]) + ") <- "
        # s += str(q) + "-"
        q = qianqu[q]
    qs.append(0)
    # print(qianqu[485])
    # print(s)
    # print(qs)
    x1, y1, z1 = [], [], []
    x2, y2, z2 = [], [], []
    x3, y3, z3 = [], [], []
    for p in points:
        if int(p[4]) == 0:
            x2.append(float(p[1]))
            y2.append(float(p[2]))
            z2.append(float(p[3]))
        else:
            x3.append(float(p[1]))
            y3.append(float(p[2]))
            z3.append(float(p[3]))            
    
    for qi in qs:
        x1.append(float(points[qi][1]))
        y1.append(float(points[qi][2]))
        z1.append(float(points[qi][3]))

    fig = plt.figure()
    # print(z1)
    # ax = Axes3D(fig)
    
    ax = fig.gca(projection='3d')
    # # 添加坐标轴(顺序是Z, Y, X)
    ax.set_zlabel('Z', fontdict={'size': 15, 'color': 'red'})
    ax.set_ylabel('Y', fontdict={'size': 15, 'color': 'red', 'label':'??'})
    ax.set_xlabel('X', fontdict={'size': 15, 'color': 'red'})
    ax.scatter(x2, y2, z2, c = 'g')
    ax.scatter(x3, y3, z3, c = 'b')

    
    ax.plot(x1, y1, z1, label='route', marker='o', mec='r', mfc='w')

    ax.legend()
    plt.show()


def plot_route_DUBIN(distance, qianqu, points, pointsType, tujing_num, graph_line_points, angles):
    # print(distance)
    # print("--------")
    # print(qianqu)
    # print("++++++++")
    # print(tujing_num)
    s = ""
    
    q = qianqu[-1]
    qs = []
    qs.append(len(points)-1)
    while q!=0:
        qs.append(q)
        typ = ""
        
        s += str(q) + "("+ str(pointsType[q]) + ") <- "
        # s += str(q) + "-"
        q = qianqu[q]
    qs.append(0)

    qs_new = []
    for i in range(len(qs)):
        qs_new.append(qs[-1-i])
    qs = qs_new
    # print(qianqu[485])
    # print(s)
    # print(qs)
    x1, y1, z1 = [], [], []
    x2, y2, z2 = [], [], []
    x3, y3, z3 = [], [], []
    x4, y4, z4 = [], [], []
    for p in points:
        if int(p[4]) == 0:
            x2.append(float(p[1]))
            y2.append(float(p[2]))
            z2.append(float(p[3]))
        else:
            x3.append(float(p[1]))
            y3.append(float(p[2]))
            z3.append(float(p[3]))            
    
    # for qi in qs:
    #     x1.append(float(points[qi][1]))
    #     y1.append(float(points[qi][2]))
    #     z1.append(float(points[qi][3]))
    


    fig = plt.figure()
    # print(z1)
    # ax = Axes3D(fig)
    
    ax = fig.gca(projection='3d')
    # # 添加坐标轴(顺序是Z, Y, X)
    ax.set_zlabel('Z', fontdict={'size': 15, 'color': 'red'})
    ax.set_ylabel('Y', fontdict={'size': 15, 'color': 'red', 'label':'??'})
    ax.set_xlabel('X', fontdict={'size': 15, 'color': 'red'})
    # ax.scatter(x2, y2, z2, c = 'g')
    # ax.scatter(x3, y3, z3, c = 'b')


    # route_arr = np.load("./dubins_routes/d"+str(DATA)+"/"+str(qi_dian)+"to"+str(zhong_dian)+".npy")
    for i in range(len(qs)-1):
        
        qi_dian = qs[i]
        zhong_dian = qs[i+1]
        print(qi_dian)
        # route_arr = np.load("./dubins_routes/d"+str(DATA)+"/"+str(qi_dian)+"to"+str(zhong_dian)+".npy")
        # print(route_arr[0])
        # y1.append(route_arr[1])
        # z1.append(route_arr[2])
        L12, LSolution = get_dubin_L(float(points[qi_dian][1]),float(points[qi_dian][2]),float(points[qi_dian][3]),\
            angles[qi_dian][zhong_dian][0],float(points[zhong_dian][1]),float(points[zhong_dian][2]),\
                float(points[zhong_dian][3]),angles[qi_dian][zhong_dian][1])
        
        # path_dubins_airplane, path_dubins_airplane_r, path_dubins_airplane_l  = ExtractDubinsAirplanePath( LSolution)
        
        # path_dubins_airplane, path_dubins_airplane_r, path_dubins_airplane_l = path_dubins_airplane.T, path_dubins_airplane_r.T, path_dubins_airplane_l.T

        # x1.extend(path_dubins_airplane_r[:, 0])
        # y1.extend(path_dubins_airplane_r[:, 1])
        # z1.extend(path_dubins_airplane_r[:, 2])

        # x4.extend(path_dubins_airplane_l[:, 0])
        # y4.extend(path_dubins_airplane_l[:, 1])
        # z4.extend(path_dubins_airplane_l[:, 2])
        # break
        paths, types = ExtractDubinsAirplanePath(LSolution)

        # paths = paths.T
        for i,p in enumerate(paths):
            p = p.T
            if types[i] == 1:
                ax.plot(p[:,0],p[:,1],p[:,2], marker=',',color='g')
            else:
                ax.plot(p[:,0],p[:,1],p[:,2], marker=',',color='r')


    
    # ax.plot(x1, y1, z1, label='my route', marker=',',color='g')
    # ax.plot(x4, y4, z4, label='my route', marker=',',color='r')



    ax.legend()
    # plt.savefig('./test.jpg')
    # print("here")
    plt.show()



# def gen_x1_to_x2_dubins()


if __name__ == '__main__':
    # graph_list = [ [0, 9, MAX_value, MAX_value, MAX_value, 14,15,MAX_value],
    #                 [9, 0, 24, MAX_value, MAX_value, MAX_value,MAX_value,MAX_value],
    #                 [MAX_value, 24, 0, 6, 2, 18,MAX_value,19],
    #                 [MAX_value, MAX_value, 6, 0, 11,MAX_value,MAX_value, 6],
    #                 [MAX_value,MAX_value, 2, 11, 0, 30,20, 16],
    #                 [14,MAX_value,18,MAX_value,30,0,5,MAX_value],
    #                 [15,MAX_value,MAX_value,MAX_value,20,5,0,44],
    #                 [MAX_value,MAX_value,19,6,16,MAX_value,44,0]]

    
    
    if INCLUDE_DUBIN:
        graph_list, pointsType, points, jiaozhengType, graph_line_points, angles = read_dataset_dubin(DATASET)
    else:
        graph_list, pointsType, points, jiaozhengType = read_dataset(DATASET)

    
    distance, qianqu, tujing_num = dijkstra(graph_list, pointsType, 0)

    if METHOD == 1:
        distance_1, res, tujing_num_1 = dijkstra_with_stepnum(graph_list, pointsType, 0, max(tujing_num),0, jiaozhengType )
    else:   
        res, distance_1 = astar(graph_list, pointsType, 0, jiaozhengType, max(tujing_num))
    # res = astar(graph_list, pointsType, 0, jiaozhengType, max(tujing_num))
    # print(res)

    route = get_route(res, points, pointsType, -1)
    print("\n")
    print("route:",route)
    print("P:",get_safe_P(route,0,0,graph_list,pointsType,jiaozhengType))
    print("distance:",distance_1[-1])

    if INCLUDE_DUBIN:
        plot_route_DUBIN([], res, points, pointsType, [], graph_line_points, angles)
    else:
        plot_route([], res, points, pointsType, [])
    
    # distance, qianqu, tujing_num = dijkstra(graph_list, pointsType, 0)
    
    # distance_1, qianqu_1, tujing_num_1 = dijkstra_with_stepnum(graph_list, pointsType, 0, max(tujing_num),0, jiaozhengType )
    
    # # print(distance)
    # # print(graph_list[0])
    # # print(pointsType)
    # # plot_route(distance, qianqu, points, pointsType, tujing_num)
    # route = get_route(qianqu_1,points,pointsType)

    # # print(route)
    # # print(route)

    # # print("-------++++++++++")
    # print(get_safe_P2(route,0,0,graph_list,pointsType,jiaozhengType))

    # plot_route(distance_1, qianqu_1, points, pointsType, tujing_num_1)
    # # plot_route(distance, qianqu, points, pointsType, tujing_num)
    
    

    
    

    
