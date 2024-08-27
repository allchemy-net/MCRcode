import argparse, sys
from rdkit import Chem
import loadDB, performRx, paths, score


def parseArgs():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rxdb", type=str, default='reactions.csv',
                    help='file with reaction database in csv format with tab separated fields (aka tsv)')
    ap.add_argument('--sbs', type=str, help='dot separated list of initial substrates', required=True)
    ap.add_argument("--gens", type=int, default=1, help='number of synthetic generation(s) to perform')
    ap.add_argument("--wrongprods", type=str, default='wrong_smarts_list.txt', help='file with motifs forbidden in reaction product')
    ap.add_argument('--verbose', type=int, default=0, help='verbosity level')
    return ap.parse_args()


def performCalc(args):
    rxdb = loadDB.loadDB(args.rxdb)
    if args.verbose > 1:
        print(f'loaded {len(rxdb)} reaction(s)')
    substrates = [Chem.CanonSmiles(smi) for smi in args.sbs.split('.')]
    initSubstrates = substrates.copy()
    matches = performRx.findMatches(substrates, rxdb, args=args)
    if args.verbose > 1:
        print("MATCHES", matches)
    if not matches:
        return None, None, rxdb, initSubstrates
    newProds = performRx.performRxes(substrates, rxdb, matches, None, args)
    allProducts = dict()
    allSteps = []
    allSteps.extend(newProds)
    for prdRx in newProds:
        prd = prdRx[0]
        if prd not in allProducts:
            allProducts[prd] = [(prdRx, ), ]
        else:  # another path
            allProducts[prd].append((prdRx, ))
    if args.verbose:
        print("=> gen 1 in::", initSubstrates, "MATCHES", matches, "OUT", newProds)
    oldMatches = dict()
    for gen in range(2, args.gens + 1):
        oldMatches = updateMatches(oldMatches, matches)
        # oldMatches.update(matches)
        thisGenBegin = len(substrates)
        newSubstrates = [res[0] for res in newProds]
        substrates.extend(newSubstrates)
        matches = performRx.findMatches(substrates, rxdb, initPos=thisGenBegin, args=args)
        if args.verbose > 5:
            print(f"matches in {gen} generation: {matches} OLD: {oldMatches}")
        newProds = performRx.performRxes(substrates, rxdb, matches, oldMatches, args)
        allSteps.extend(newProds)
        for prdRx in newProds:
            newPaths = paths.formPaths(prdRx, allProducts)
            if prdRx[0] not in allProducts:
                allProducts[prdRx[0]] = []
            allProducts[prdRx[0]].extend(newPaths)
        if args.verbose:
            print("=> gen", gen, "in::", newSubstrates, "MATCHES", matches, "OLD", oldMatches, "OUT", newProds)
        if args.verbose > 10:
            print(f'all products after {gen} generation(s):')
            for prd in allProducts:
                print('    ', prd)
    return allProducts, allSteps, rxdb, initSubstrates


def updateMatches(oldMatches, matches):
    for rx in matches:
        if rx not in oldMatches:
            oldMatches[rx] = matches[rx]
            continue
        for poz in matches[rx]:
            if poz not in oldMatches[rx]:
                oldMatches[rx][poz] = matches[rx][poz]
                continue
            oldMatches[rx][poz].extend(matches[rx][poz])
    return oldMatches


def printPaths(allPaths):
    print('calculation results with paths:')
    print('--------------------------------')
    for cmd in allPaths:
        if not allPaths[cmd]:
            print(cmd, 'is initial substrate')
        else:
            pathLens = [len(path[0]) for path in allPaths[cmd]]
            print(f'\n {cmd} has {len(allPaths[cmd])} path(s) of length {pathLens}')
            for pid, path in enumerate(allPaths[cmd], 1):
                print(f' path {pid}, info: {path[1]["pathType"]};  Steps:')
                for step in path[0]:
                    stepInfo = {'product': step[0], 'substrates': step[1], 'rxid': step[2],
                                 'byproducts': step[3], 'byproductsReactants': step[4]}
                    print(f'   {stepInfo}')


if __name__ == "__main__":
    args = parseArgs()
    allPaths, allSteps, rxdb, initSubstrates = performCalc(args)
    if not allPaths:
        print('no product formed')
        sys.exit()
    allPaths = paths.filterPaths(allPaths, allSteps, rxdb, args)
    if not allPaths:
        print('no product formed (after path selection)')
    allPaths = score.scorePaths(allPaths, allSteps, rxdb, args)
    if args.verbose > 20:
        print('all single steps')
        for prd in allPaths:
            print('  ', prd, allPaths[prd])
    printPaths(allPaths)
