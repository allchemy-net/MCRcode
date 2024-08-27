import itertools
import networkx


def formPaths(prodInfo, allProducts):
    prod, sbses, rxid = prodInfo
    sbsesNotInit = [sbs for sbs in sbses if sbs in allProducts]
    allPathsNoLast = [allProducts[sbs] for sbs in sbsesNotInit]
    allPaths = []
    for pathNoLast in itertools.product(*allPathsNoLast):
        path = []
        for steps in pathNoLast:
            path.extend(steps)
        path.append(prodInfo)
        allPaths.append(path)
    return allPaths


def filterPaths(paths, steps, rxdb, args):
    paths = filterOutDueToVeryFast(paths, steps, rxdb, args)
    paths = filterOutDueToOxRed(paths, rxdb, args)
    paths = filterOutDueCondOrder(paths, rxdb, args)
    return paths


def filterOutDueToVeryFast(paths, steps, rxdb, args):
    veryFast = set()
    for rx in rxdb:
        if rx['RATE'].upper() == 'VF':
            veryFast.add(rx['rxid'])
    if not veryFast:
        return paths
    graph = buildGraph(steps)
    if args.verbose > 10:
        print("detected VERY FAST reactions in database", veryFast)
    rxToRemove = getRxToRemove(graph, veryFast, args)
    if rxToRemove:
        paths = removePathWithGivenSteps(paths, rxToRemove, args)
    return paths


def buildGraph(steps):
    # format step := (products, (direct, substrate, list), rxid)
    graph = networkx.DiGraph()
    for step in steps:
        prod, sbses, rxid = step
        # rxhash = (rxid}::{".".join(sbses)}>>{prod}'
        rxhash = (prod, tuple(sbses), rxid)
        for sbs in sbses:
            graph.add_edge(sbs, rxhash)
        graph.add_edge(rxhash, prod)
    return graph


def getRxToRemove(graph, veryFast, args):
    reactionToRemove = []
    for reaction in graph.nodes:
        if isinstance(reaction, str):
            continue  # compound
        if reaction[-1] not in veryFast:
            continue  # normal reaction
        for sbs in graph.predecessors(reaction):
            rxToRm = [rx for rx in graph.succ[sbs] if rx != reaction]
            if rxToRm:
                reactionToRemove.extend(rxToRm)
            if args.verbose > 10:
                print('reaction to remove due to speed', rxToRm, 'reason:', reaction)
    return reactionToRemove


def removePathWithGivenSteps(paths, rxToRemove, args):
    rxToRm = set(rxToRemove)
    for cmd in list(paths.keys()):
        newPaths = []
        for path in paths[cmd]:
            if set(path).intersection(rxToRm):
                if args.verbose > 5:
                    print(cmd, "REMOVE PATH", path, "HAS", rxToRemove, "OF", len(paths[cmd]))
            else:
                newPaths.append(path)
        if not newPaths:
            if args.verbose:
                print("remove compound", cmd, 'due to competing very fast reaction(s)')
            del paths[cmd]
        else:
            paths[cmd] = newPaths
    return paths


def filterOutDueToOxRed(paths, rxdb, args):
    reductive = set()
    oxidative = set()
    for rx in rxdb:
        if rx['OX/RED'].upper().startswith('O'):
            oxidative.add(rx['rxid'])
        if rx['OX/RED'].upper().startswith('R'):
            reductive.add(rx['rxid'])
    if not (reductive and oxidative):
        return paths
    cmdToRm = []
    for cmd in paths:
        newPaths = []
        for path in paths[cmd]:
            rxOnPath = set([rxinfo[-1] for rxinfo in path])
            if rxOnPath.intersection(reductive) and rxOnPath.intersection(oxidative):
                if args.verbose > 5:
                    print(f'remove path to {cmd} :: {path} due to ox and red on path')
            else:
                newPaths.append(path)
        if not newPaths:
            if args.verbose:
                print("remove compound", cmd, 'due to competing very fast reaction(s)')
            cmdToRm.append(cmd)
        else:
            paths[cmd] = newPaths
    for cmd in cmdToRm:
        if args.verbose > 5:
            print(f"{cmd} removed - no valid paths!")
        del paths[cmd]
    return paths


def filterOutDueCondOrder(paths, rxdb, args):
    rxdbDict = {rx['rxid']: rx for rx in rxdb}
    cmdToRm = []
    for cmd in paths:
        newPaths = []
        for path in paths[cmd]:
            pathGraph = buildPathGraph(path)
            allTopoSort = getRxidFromTopoSort(networkx.all_topological_sorts(pathGraph))
            for topoSort in allTopoSort:
                isWaterOk = checkTopoOrderWater(topoSort, rxdbDict)
                isProticOk = checkTopoOrderProtic(topoSort, rxdbDict)
                isTemperatureOk = checkTopoOrderTemperature(topoSort, rxdbDict)
                isABok = checkTopoOrderAB(topoSort, rxdbDict)
                if isWaterOk and isProticOk and isTemperatureOk and isABok:
                    newPaths.append(path)
                    break  # path has topological sort which can be perform experimentally
        if newPaths:
            if args.verbose > 10:
                print(f"allowed path for {cmd} after conditions check {newPaths}")
            paths[cmd] = newPaths
        else:
            cmdToRm.append(cmd)
    for cmd in cmdToRm:
        if args.verbose > 5:
            print(f"{cmd} removed due to conditions!")
        del paths[cmd]
    return paths


def checkTopoOrderWater(topoSort, rxdb):
    waterCond = [rxdb[rx]['W/WF'] for rx in topoSort]
    if 'W' in waterCond and 'WF' in waterCond:
        W = [poz for poz, wc in enumerate(waterCond) if wc == 'W']
        WF = [poz for poz, wc in enumerate(waterCond) if wc == 'WF']
        if min(W) < max(WF):
            return False
    return True


def checkTopoOrderProtic(topoSort, rxdb):
    proticCond = [rxdb[rx]['protic/aprotic'] for rx in topoSort]
    if 'P' in proticCond and 'AP' in proticCond:
        P = [poz for poz, pc in enumerate(proticCond) if pc == 'P']
        AP = [poz for poz, pc in enumerate(proticCond) if pc == 'AP']
        if min(P) < max(AP):
            return False
    return True


def isExtrema(lst, poz):
    if lst[poz - 1] > lst[poz] and lst[poz + 1] > lst[poz]:
        return True  # local minima
    if lst[poz - 1] < lst[poz] and lst[poz + 1] < lst[poz]:
        return True  # local maxima
    return False


def checkTopoOrderTemperature(topoSort, rxdbDict):
    tempToNum = {'VL': -2, 'L': -1, 'rt': 0, 'H': 1, 'VH': 2}
    tempCond = [rxdbDict[rx]['temperature'].split('.') for rx in topoSort]
    for tempVariant in itertools.product(*tempCond):
        tempNums = [tempToNum[temp] for temp in tempVariant]
        tempNumsNoDupli = []
        for num in tempNums:
            if not tempNumsNoDupli:
                tempNumsNoDupli.append(num)
            elif num != tempNumsNoDupli[-1]:
                tempNumsNoDupli.append(num)
        extrema = [poz for poz in range(1, len(tempNumsNoDupli) - 1) if isExtrema(tempNumsNoDupli, poz)]
        if len(extrema) < 2:
            return True
    return False


def checkTopoOrderAB(topoSort, rxdbDict):
    ABtoNum = {'SB': -3, 'B': -2, 'WB': -1, 'N': 0, 'WA': 1, 'A': 2, 'SA': 3}
    ABcond = [rxdbDict[rx]['conditions'].replace('.LA', '').split('.') for rx in topoSort]
    for ABvariant in itertools.product(*ABcond):
        if 'SB' in ABvariant and ('A' in ABvariant or 'SA' in ABvariant or 'WA' in ABvariant):
            sbPos = ABvariant.index('SB')
            if 'A' in ABvariant and ABvariant.index('A') < sbPos:
                continue  # A then SB is not allowed
            if 'WA' in ABvariant and ABvariant.index('WA') < sbPos:
                continue  # WA then SB not allowed
            if 'SA' in ABvariant and ABvariant.index('SA') < sbPos:
                continue  # SA then SB not allowed
        ABnums = [ABtoNum[ab] for ab in ABvariant]
        ABuniq = []
        for ab in ABnums:
            if not ABuniq:
                ABuniq.append(ab)
            elif ABuniq[-1] != ab:
                ABuniq.append(ab)
        extrema = [poz for poz in range(1, len(ABuniq) - 1) if isExtrema(ABuniq, poz)]
        if not extrema:
            return True
    return False


def buildPathGraph(path):
    graph = networkx.DiGraph()
    for step in path:
        prod, sbses, rx = step
        for sbs in sbses:
            graph.add_edge(sbs, step)
        graph.add_edge(step, prod)
    return graph


def getRxidFromTopoSort(allTopoSorts):
    allTopos = []
    for topoSort in allTopoSorts:
        newTopo = []
        for node in topoSort:
            if isinstance(node, str):
                continue  # compound node
            newTopo.append(node[-1])
        allTopos.append(tuple(newTopo))
    return allTopos
