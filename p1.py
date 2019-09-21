
import math

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  

MAX_value = 999999
DELTA = 0.001
THETA = 30
ALPHA1 = 25
ALPHA2 = 15
BETA1 = 20
BETA2 = 25

ALPHA = 0.1
BETA = 0.1
GAMA =  0.1

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

def dijkstra_normal(graph, pointType, s):
    # 判断图是否为空，如果为空直接退出
    # pointType: -1:起始点,0:水平点,1:垂直点,2:终点
    if graph is None:
        return None
    dist = [MAX_value]*len(graph)
    dist[s] = 0

    # shuiping = [MAX_value] * len(graph)
    # chuizhi = [MAX_value] * len(graph)

    # shuiping[s] = 0
    # chuizhi[s] = 0

    qianqu = [MAX_value] * len(graph)

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
                    dist[v] = dist[u] + d
                    dist_init[v] = dist[v]
                    qianqu[v] = u
    return dist, qianqu 
     
 

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


def dijkstra_with_stepnum(graph, pointType, s, max_stepnum, func_index):
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
                # if dist[v] > dist[u]+d:
                if ALPHA * (dist[v]/graph[0][-1]) + BETA * funcs[func_index]((tujing_num[v]/max_stepnum)) > ALPHA * ((dist[u] + d)/graph[0][-1]) + BETA * funcs[func_index](((tujing_num[u] + 1)/max_stepnum)):
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
 
 
 

def read_dataset(file):
    f = open("dataset.csv","r")
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

    return graph, pointType, points




def plot_route(distance, qianqu, points, pointsType, tujing_num):
    print(distance)
    print("--------")
    print(qianqu)
    print("++++++++")
    print(tujing_num)
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
    print(qs)
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
    print(z1)
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



if __name__ == '__main__':
    # graph_list = [ [0, 9, MAX_value, MAX_value, MAX_value, 14,15,MAX_value],
    #                 [9, 0, 24, MAX_value, MAX_value, MAX_value,MAX_value,MAX_value],
    #                 [MAX_value, 24, 0, 6, 2, 18,MAX_value,19],
    #                 [MAX_value, MAX_value, 6, 0, 11,MAX_value,MAX_value, 6],
    #                 [MAX_value,MAX_value, 2, 11, 0, 30,20, 16],
    #                 [14,MAX_value,18,MAX_value,30,0,5,MAX_value],
    #                 [15,MAX_value,MAX_value,MAX_value,20,5,0,44],
    #                 [MAX_value,MAX_value,19,6,16,MAX_value,44,0]]
    
    graph_list, pointsType, points = read_dataset("dataset.csv")
    distance, qianqu, tujing_num = dijkstra(graph_list, pointsType, 0)
    
    distance_1, qianqu_1, tujing_num_1 = dijkstra_with_stepnum(graph_list, pointsType, 0, max(tujing_num),0 )
    
    # print(distance)
    # print(graph_list[0])
    # print(pointsType)
    # plot_route(distance, qianqu, points, pointsType, tujing_num)
    plot_route(distance_1, qianqu_1, points, pointsType, tujing_num_1)

    
    

    
