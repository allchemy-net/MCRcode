from rdkit import Chem
from rdkit.Chem import AllChem
import scoreHelpers


def scorePaths(paths, steps, rxdb, args):
    rxinfo = formRxdb(rxdb)
    newPathsForAllCmds = dict()
    for cmd in paths:
        #print(cmd,"==>",  paths[cmd])
        newpaths = []
        for path in paths[cmd]:
            newp, info = scoreSinglePath(path, steps, rxinfo, args)
            newpaths.append((newp, info))
        newPathsForAllCmds[cmd] = newpaths
    return newPathsForAllCmds


def formRxdb(rxdb):
    rxinfo = dict()
    for rx in rxdb:
        rxinfo[rx['rxid']] = rx
    return rxinfo


def scoreSinglePath(path, steps, rxinfo, args):
    sbs, prods, manyuse = getNetSbsAndProd(path)
    #print("SBS", sbs, "PR", prods)
    # path = doMappingPath(path, sbs, manyuse, rxinfo)
    path = doAddByprods(path, sbs, rxinfo)
    pathType = calculatePathType(path, rxinfo, args)
    info = {'pathType': pathType}
    return path, info


def getNetSbsAndProd(path):
    rea = dict()
    for rx in path:
        prod, sbses, _ = rx
        if prod in rea:
            rea[prod][1] += 1
        else:
            rea[prod] = [0, 1]
        for sbs in sbses:
            if sbs in rea:
                rea[sbs][0] += 1
            else:
                rea[sbs] = [1, 0]
    sbses, prods, manyuse = dict(), dict(), dict()
    for smi in rea:
        nsbs, nprod = rea[smi]
        if nsbs == nprod:
            continue
        if nsbs and nprod:
            manyuse[smi] = [nsbs, nprod]
        if nsbs:
            sbses[smi] = nsbs
        if nprod:
            prods[smi] = nprod
    return sbses, prods, manyuse


def doAddByprods(path, sbs, rxinfo):
    pathByprod = []
    for step in path:
        prod, sbses, rxid = step
        mols = [Chem.MolFromSmiles(smi) for smi in sbses]
        rxdata = rxinfo[rxid]
        rxn = AllChem.ReactionFromSmarts(rxdata['reaction smarts'])
        prods = rxn.RunReactants(mols)
        prods = prods[0]
        _ = [Chem.SanitizeMol(m) for m in prods]
        prods = [Chem.MolToSmiles(mol) for mol in prods]
        pathByprod.append((prod, sbses, rxid, prods[1:], rxinfo[rxid]['by-products'].split('.')))
    return pathByprod


def calculatePathType(path, rxinfo, args):
    pathType = None
    problems, metals = scoreHelpers.detectProblems(path, rxinfo, args)
    pathType = scoreHelpers.getPathType(path, args)
    #warnings = [step['problems']['other'] for step in pathinfo if step['problems']['other']]
    # values = ( 'onepotweak', 'onepotstrong', 'multicomponentstrong', 'multicomponentweak')
    #ismulti, whynot, trueUsedSbs, allInitSbs, multiUsageDict, reusedSbs = scoreHelpers.pathRestriction.isPathOkForMulticomponent(path, all_rxs, smiles_helper)
    #res = scoreHelpers.getMCRcategory(ismulti, whynot, trueUsedSbs, allInitSbs, byProblem, warnings, multiUsageDict, reusedSbs)
    pathType.update({'crossreactivity': problems, 'metals': metals})
    return pathType
