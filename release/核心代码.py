def dijkstra_with_stepnum(graph, pointType, s, max_stepnum, func_index, jiaozhengType):
    # 判断图是否为空，如果为空直接退出
    # pointType: -1:起始点,0:水平点,1:垂直点,2:终点
    # 改进的Dijkstra算法
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
    # 改进的A*算法
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
        min_f = G_P *  G[min_d] + H_P * ((H[min_d]+0.0)/graph[0][-1])
        for v,d in enumerate(open_list):
            if G_P * G[d] + H_P * ((H[d] + 0.0)/graph[0][-1]) < min_f:
                min_f = G_P * G[d] + H_P * ((H[d] + 0.0)/graph[0][-1])
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
                    H[i] = ((graph[i][-1] + 0.0)/graph[0][-1])
                else:
                    G_current = ALPHA * ((Dist_[node_current] + d)/graph[0][-1]) + BETA * funcs[FUNC_INDEX](((Cishu_[node_current] + 1.0)/max_time)) + GAMA * sp_
                    if G_current < G[i]:
                        parent[i] = node_current
                        G[i] = G_current
                        shuiping[i] = 0
                        chuizhi[i] = cz
                        Dist_[i] = Dist_[node_current] + d
                        Cishu_[i] = Cishu_[node_current] + 1
                        H[i] = ((graph[i][-1] + 0.0)/graph[0][-1])
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
                    H[i] = ((graph[i][-1] + 0.0)/graph[0][-1])
                else:
                    G_current = ALPHA * ((Dist_[node_current] + d)/graph[0][-1]) + BETA * funcs[FUNC_INDEX](((Cishu_[node_current] + 1.0)/max_time)) + GAMA * sp_
                    if G_current < G[i]:
                        parent[i] = node_current
                        G[i] = G_current
                        shuiping[i] = sp
                        chuizhi[i] = 0
                        Dist_[i] = Dist_[node_current] + d
                        Cishu_[i] = Cishu_[node_current] + 1
                        H[i] = ((graph[i][-1] + 0.0)/graph[0][-1])
                
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
                    H[i] = ((graph[i][-1] + 0.0)/graph[0][-1])
                else:
                    G_current = ALPHA * ((Dist_[node_current] + d)/graph[0][-1]) + BETA * funcs[FUNC_INDEX](((Cishu_[node_current] + 1.0)/max_time)) + GAMA * sp_
                    if G_current < G[i]:
                        parent[i] = node_current
                        G[i] = G_current
                        shuiping[i] = 0
                        chuizhi[i] = 0
                        Dist_[i] = Dist_[node_current] + d+ GAMA * sp_
                        Cishu_[i] = Cishu_[node_current] + 1    
                        H[i] = ((graph[i][-1] + 0.0)/graph[0][-1])
            else:
                #不可达
                continue
    
    return parent, Dist_

def get_safe_P(routes, chushi_shuiping, chushi_chuizhi, graph, pointType, jiaozhengType):
    # 计算一条路径的安全概率
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

