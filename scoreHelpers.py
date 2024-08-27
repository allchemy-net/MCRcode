from rdkit import Chem
import networkx as nx

problematic_metals = {'[Na', '[K', '[Li', '[Rb', '[Cs', '[Be', '[Mg', '[Ca', '[Sr', '[Ba',
                      '[Sc', '[Y', '[Ti', '[Zr', '[V', '[Cr', '[Mo', '[W', '[Mn', '[Re', '[Fe', '[Ru', '[Os', '[Co',
                      '[Rh', '[Ir', '[Ni', '[Pd', '[Pt', '[Cu', '[Ag', '[Au',
                      '[Zn', '[Cd', '[Hg', '[Al', '[Ga', '[In', '[Tl', '[Ge', '[Sn', '[Pb', '[Bi', '[Sb',
                      '[Ce', '[Pr', '[Nd', '[Pm', '[Sm', '[Eu', '[Gd', '[U',
                      }


def getPathType(path, args):
    allPSB = dict()  # all Products Substrates Byproducts of all steps
    for rx in path:
        prod, sbses, _, byprods, _ = rx
        if prod not in allPSB:
            allPSB[prod] = [0, 0, 0]
        allPSB[prod][0] += 1
        for sbs in sbses:
            if sbs not in allPSB:
                allPSB[sbs] = [0, 0, 0]
            allPSB[sbs][1] += 1
        for byprod in byprods:
            if byprod not in allPSB:
                allPSB[byprod] = [0, 0, 0]
            allPSB[byprod][2] += 1
    # reduce normal intermediate
    for smi in tuple(allPSB.keys()):
        if allPSB[smi] == [1, 1, 0]:
            del allPSB[smi]
    trueSbs = dict()
    trueProd = dict()
    trueByprod = dict()
    others = dict()
    autocat = dict()
    for smi in allPSB:
        if allPSB[smi][0] > 0 and allPSB[smi][1] == 0 and allPSB[smi][2] == 0:
            trueProd[smi] = allPSB[smi][0]
        elif allPSB[smi][0] == 0 and allPSB[smi][1] > 0 and allPSB[smi][2] == 0:
            trueSbs[smi] = allPSB[smi][1]
        elif allPSB[smi][0] == 0 and allPSB[smi][1] == 0 and allPSB[smi][2] > 0:
            trueByprod[smi] = allPSB[smi][2]
        elif allPSB[smi][0] == 0 and allPSB[smi][1] > 0 and allPSB[smi][2] > 0:
            autocat[smi] = allPSB[smi]
        else:
            others[smi] = allPSB[smi]
    if args.verbose > 10:
        print("Substrates:", trueSbs, "Product:", trueProd, "OTHER:", others)
    return {'numberOfComponents': len(trueSbs), 'catalytic':len(autocat)}


def detectProblems(path, rxinfo, args):
    usedMetals = _calcMetalsOnPath(path)
    graph = buildRxGraph(path)
    problems = dict()
    allSbs = set()
    allProds = set()
    for step in path:
        allProds.add(step[0])
        _ = [allSbs.add(smi) for smi in step[1]]
    trueSbs = allSbs - allProds
    for poz, step in enumerate(path):
        problem = _detectProblemInStep(path, poz, rxinfo, graph, usedMetals, trueSbs)
        if problem:
            if args.verbose > 5:
                print(f'problem for MCR detected in {poz} {step}')
            problems[poz] = problem
    metalProblem = _addWarningAboutMetals(usedMetals)
    return problems, metalProblem


def buildRxGraph(path):
    graph = nx.DiGraph()
    for poz, rxdata in enumerate(path):
        prod, sbses = rxdata[0], rxdata[1]
        for sbs in sbses:
            graph.add_edge(sbs, poz)
        graph.add_edge(poz, prod)
    return graph


def _calcMetalsOnPath(histToShow):
    metals = dict()
    for step in histToShow:
        for elem in step[3] + step[4]:
            for met in problematic_metals:
                if met in elem:
                    metal = met[1:]
                    if metal not in metals:
                        metals[metal] = set()
                    metals[metal].add(elem)
    return metals


def _detectProblemInStep(path, poz, rxinfo, graph, usedMetals, trueSbs):
    prevs = [node for node in nx.ancestors(graph, poz) if isinstance(node, int)]
    donerx = [rxinfo[path[poz][2]], ]
    problems = []
    for prev in prevs:
        for smi in path[prev][3]:
            if _isSmilesMakesProblem(smi, donerx):
                problems.append(smi)
    prevRxes = [rxinfo[path[stepoz][2]] for stepoz in prevs]
    for smi in set(path[poz][1]).intersection(trueSbs):
        if _isSmilesMakesProblem(smi, prevRxes):
            problems.append(smi)
    return problems


def _addWarningAboutMetals(usedMetals):
    if len(usedMetals) > 0:
        allMetalsStr = '.'.join(usedMetals.keys())
        allMetalsForm = set()
        for metal in usedMetals:
            allMetalsForm.update(usedMetals[metal])
        allMetalsFormStr = '.'.join(allMetalsForm)
        return f"{len(usedMetals)} metals in mixture: {allMetalsStr} in {len(allMetalsForm)} forms: {allMetalsFormStr} :{usedMetals} "
    return ''


def _isSmilesMakesProblem(sidesmiles, donerxes):
    if not sidesmiles:
        return False
    sidemol = Chem.MolFromSmiles(sidesmiles)
    if not sidemol:
        return False
    for rx in donerxes:
        incopatts = rx['incompatible groups'].split('.')
        for inco in incopatts:
            incomol = Chem.MolFromSmarts(inco)
            if sidemol.HasSubstructMatch(incomol):
                return True
    return False

