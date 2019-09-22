
import math

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  

MAX_value = 99999999
DELTA = 0.001

ALPHA1 = 25
ALPHA2 = 15
BETA1 = 20
BETA2 = 25
THETA = 30

# ALPHA1 = 20
# ALPHA2 = 10
# BETA1 = 15
# BETA2 = 20
# THETA = 20

ALPHA = 0.1
BETA = 0
GAMA = -10

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
        print(u)
        print(len(Q),'/',len(graph))
        for v, d in enumerate(graph[u]):
            if 0 < d < MAX_value and v in Q:
                # if dist[v] > dist[u]+d:
                # print(u,v)
                route1, route2 = get_two_routes(u,v,qianqu)
                p1,p2 = get_safe_P(route1,0,0,graph,pointType,jiaozhengType),get_safe_P(route2,0,0,graph,pointType,jiaozhengType)
                # print(route1,route2,p1,p2)
                if ALPHA * (dist[v]/graph[0][-1]) + BETA * funcs[func_index]((tujing_num[v]/max_stepnum)) + GAMA * p2 > \
                     ALPHA * ((dist[u] + d)/graph[0][-1]) + BETA * funcs[func_index](((tujing_num[u] + 1)/max_stepnum)) + GAMA * p1:
                    if pointType[v] == 0 and shuiping[u] + d * DELTA < BETA2 \
                        and chuizhi[u] + d * DELTA < BETA1:
                        # if jiaozhengType[v] == 1:
                        dist[v] = dist[u]+d
                        # dist_init[v] = dist[v]
                        dist_init[v] = ALPHA * (dist[v]/graph[0][-1]) + BETA * funcs[func_index]((tujing_num[v]/max_stepnum)) + GAMA * p2
                        shuiping[v] = 0
                        chuizhi[v] = chuizhi[u] + d*DELTA
                        qianqu[v] = u 
                        tujing_num[v] = tujing_num[u] + 1

                    if pointType[v] == 1 and shuiping[u] + d * DELTA < ALPHA2 \
                        and chuizhi[u] + d * DELTA < ALPHA1: 
                        dist[v] = dist[u]+d
                        # dist_init[v] = dist[v]
                        dist_init[v] = ALPHA * (dist[v]/graph[0][-1]) + BETA * funcs[func_index]((tujing_num[v]/max_stepnum)) + GAMA * p2

                        chuizhi[v] = 0
                        shuiping[v] = shuiping[u] + d * DELTA
                        qianqu[v] = u
                        tujing_num[v] = tujing_num[u] + 1

                    if pointType[v] == 2 and shuiping[u] + d* DELTA < THETA \
                        and chuizhi[u] + d * DELTA < THETA:
                        dist[v] = dist[u] + d
                        # dist_init[v] = dist[v]
                        dist_init[v] = ALPHA * (dist[v]/graph[0][-1]) + BETA * funcs[func_index]((tujing_num[v]/max_stepnum)) + GAMA * p2

                        qianqu[v] = u
                        tujing_num[v] = tujing_num[u] + 1
    
    return dist, qianqu, tujing_num
 
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

def get_safe_P2(routes, chushi_shuiping, chushi_chuizhi, graph, pointType, jiaozhengType):
    print(routes)
    if len(routes) <= 1:
        return 1
    
    # if chushi_shuiping + d * DELTA
    d = graph[routes[0]][routes[1]]
    shuiping = chushi_shuiping + d * DELTA
    chuizhi = chushi_chuizhi + d * DELTA
    if pointType[routes[1]] == 0 and shuiping < BETA2 \
        and chuizhi < BETA1:
        if jiaozhengType[routes[1]] == 0: #shuiping
            return 1 * get_safe_P2(routes[1:],0.0, chuizhi, graph, pointType, jiaozhengType)
        else:
            return 0.8* get_safe_P2(routes[1:],0.0, chuizhi, graph, pointType, jiaozhengType) + 0.2 *\
                 get_safe_P2(routes[1:],min(shuiping,5), chuizhi, graph, pointType, jiaozhengType)
    elif pointType[routes[1]] == 1 and shuiping< ALPHA2 and chuizhi < ALPHA1:
        if jiaozhengType[routes[1]] == 0:
            return 1 * get_safe_P2(routes[1:],shuiping, 0.0, graph, pointType, jiaozhengType)
        else:
            return 0.8*get_safe_P2(routes[1:],shuiping, 0.0, graph, pointType, jiaozhengType) +\
                0.2 * get_safe_P2(routes[1:],shuiping, min(chuizhi,5), graph, pointType, jiaozhengType)
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


def get_route(qianqu, points, pointsType):
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


# def gen_route_csv(route, graph_list, pointsType, points, jiaozhengType):

def astar(graph, pointsType, s, jiaozhengType):
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
        min_f = G[min_d] + H[min_d]
        for v,d in enumerate(open_list):
            if G[d] + H[d] < min_f:
                min_f = G[d] + H[d]
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
                if i not in open_list:
                    open_list.append(i)
                    parent[i] = node_current
                    shuiping[i] = 0
                    chuizhi[i] = cz
                    Dist_[i] = Dist_[node_current] + d
                    Cishu_[i] = Cishu_[node_current] + 1
                    G[i] = ALPHA * Dist_[i] + BETA * Cishu_[i]
                    H[i] = graph[i][-1]
                else:
                    G_current = ALPHA * (Dist_[node_current] + d) + BETA * (Cishu_[node_current] + 1)
                    if G_current < G[i]:
                        parent[i] = node_current
                        G[i] = ALPHA * (Dist_[node_current] + d) + BETA * (Cishu_[node_current] + 1)
                        shuiping[i] = 0
                        chuizhi[i] = cz
                        Dist_i = Dist_[node_current] + d
                        Cishu_[i] = Cishu_[node_current] + 1
            elif pointsType[i] == 1 and sp <= ALPHA2 and cz <= ALPHA1:
                #可达
                if i not in open_list:
                    open_list.append(i)
                    parent[i] = node_current
                    shuiping[i] = sp
                    chuizhi[i] = 0
                    Dist_[i] = Dist_[node_current] + d
                    Cishu_[i] = Cishu_[node_current] + 1
                    G[i] = ALPHA * Dist_[i] + BETA * Cishu_[i]
                    H[i] = graph[i][-1]
                else:
                    G_current = ALPHA * (Dist_[node_current] + d) + BETA * (Cishu_[node_current] + 1)
                    if G_current < G[i]:
                        parent[i] = node_current
                        G[i] = ALPHA * (Dist_[node_current] + d) + BETA * (Cishu_[node_current] + 1)
                        shuiping[i] = sp
                        chuizhi[i] = 0
                        Dist_i = Dist_[node_current] + d
                        Cishu_[i] = Cishu_[node_current] + 1
                
            elif pointsType[i] == 2 and sp <= THETA and cz <= THETA:
                #可达
                if i not in open_list:
                    open_list.append(i)
                    parent[i] = node_current
                    shuiping[i] = 0
                    chuizhi[i] = 0
                    Dist_[i] = Dist_[node_current] + d
                    Cishu_[i] = Cishu_[node_current] + 1
                    G[i] = ALPHA * Dist_[i] + BETA * Cishu_[i]
                    H[i] = graph[i][-1]
                else:
                    G_current = ALPHA * (Dist_[node_current] + d) + BETA * (Cishu_[node_current] + 1)
                    if G_current < G[i]:
                        parent[i] = node_current
                        G[i] = ALPHA * (Dist_[node_current] + d) + BETA * (Cishu_[node_current] + 1)
                        shuiping[i] = 0
                        chuizhi[i] = 0
                        Dist_i = Dist_[node_current] + d
                        Cishu_[i] = Cishu_[node_current] + 1                
            else:
                #不可达
                continue
    
    return parent
            

def checklist(graph, open_list, close_list):
    if len(open_list) == 0:
        return False
    if len(graph) - 1 in open_list:
        return False
    return True

if __name__ == '__main__':
    # graph_list = [ [0, 9, MAX_value, MAX_value, MAX_value, 14,15,MAX_value],
    #                 [9, 0, 24, MAX_value, MAX_value, MAX_value,MAX_value,MAX_value],
    #                 [MAX_value, 24, 0, 6, 2, 18,MAX_value,19],
    #                 [MAX_value, MAX_value, 6, 0, 11,MAX_value,MAX_value, 6],
    #                 [MAX_value,MAX_value, 2, 11, 0, 30,20, 16],
    #                 [14,MAX_value,18,MAX_value,30,0,5,MAX_value],
    #                 [15,MAX_value,MAX_value,MAX_value,20,5,0,44],
    #                 [MAX_value,MAX_value,19,6,16,MAX_value,44,0]]
    
    graph_list, pointsType, points, jiaozhengType = read_dataset("dataset.csv")

    res = astar(graph_list, pointsType, 0, jiaozhengType)
    print(res)

    route = get_route(res, points, pointsType)
    print(route)
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
    
    

    
    

    
