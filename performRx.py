from rdkit import Chem
from rdkit.Chem import AllChem


def isProductValid(smiles, args):
    forbiddenSmarts = []
    for line in open(args.wrongprods):
        if line.strip().startswith('#'):
            continue
        for sma in line.split('.'):
            sma = sma.strip()
            if not sma:
                continue
            forbiddenSmarts.append(sma)
    forbiddenMols = [Chem.MolFromSmarts(sma) for sma in forbiddenSmarts]
    mol = Chem.MolFromSmiles(smiles)
    for patt in forbiddenMols:
        if mol.HasSubstructMatch(patt):
            return False
    return True


def findMatches(substrates, rxinfos, initPos=0, args=None):
    mols = [Chem.MolFromSmiles(sbs) for sbs in substrates]
    matches = dict()
    if args.verbose > 9:
        print(f"find matches of {len(mols)} substrate(s) to {len(rxinfos)} reaction(s)")
    for rxpos, rxinfo in enumerate(rxinfos):
        rxn = AllChem.ReactionFromSmarts(rxinfo['reaction smarts'])
        incoGroups = [Chem.MolFromSmarts(incogr) for incogr in rxinfo['incompatible groups'].split('.')]
        if args and args.verbose > 5:
            print(f"{rxinfo['rxid']} has {len(incoGroups)} not allowed group(s)")
        rxmatches = findMatchingSbsToRx(rxn, incoGroups, mols, initPos, args)
        if not rxmatches:
            if args and args.verbose > 5:
                print(f"no matches for {rxpos}, rxid: {rxinfo['rxid']}")
            continue
        matches[rxpos] = rxmatches
    if args.verbose > 10:
        print("Matches", matches)
    return matches


def findMatchingSbsToRx(rxn, incoGroups, mols, initPos, args=None):
    positions = dict()
    for pattPos, sbsPatt in enumerate(rxn.GetReactants()):
        for molPos, mol in enumerate(mols):
            if molPos < initPos:
                continue
            coreMatches = mol.GetSubstructMatches(sbsPatt)
            if not coreMatches:
                continue
            cores = set.union(*[set(core) for core in coreMatches])
            if args and args.verbose > 9:
                print(f"CORES {cores} detected in {Chem.MolToSmiles(mol)} for pattern {pattPos} i.e. {Chem.MolToSmarts(sbsPatt)}")
            detectedIncompatGroups = 0
            for incoGr in incoGroups:
                incoMatches = mol.GetSubstructMatches(incoGr)
                if not incoMatches:
                    continue
                incos = set.union(*[set(inco) for inco in incoMatches])
                if incos - cores:
                    if args and args.verbose > 5:
                        print(f'detected incompatible group {Chem.MolToSmarts(incoGr)} in substrate{Chem.MolToSmiles(mol)}')
                        if args and args.verbose > 9:
                            print(f' cores: {cores} incos: {incos}')
                    detectedIncompatGroups += 1
            if detectedIncompatGroups > 0:
                continue
            # here we have detected correct matching to reaction
            if pattPos not in positions:
                positions[pattPos] = [molPos, ]
            else:
                positions[pattPos].append(molPos)
    return positions


def performRxes(substrates, rxinfos, newMatching, oldMatching=None, args=None):
    products = []
    if args.verbose > 20:
        print('=== perform reactions ===')
    for rxpos in newMatching:
        if args and args.verbose > 5:
            print('rxPOS', rxpos, "==>", newMatching[rxpos])
        rxid = rxinfos[rxpos]['rxid']
        rxn = AllChem.ReactionFromSmarts(rxinfos[rxpos]['reaction smarts'])
        numTemplates = rxn.GetNumReactantTemplates()
        if numTemplates == 1:
            prods = [perform1cRx(rxn, substrates[sbs], rxid, args) for sbs in newMatching[rxpos][0]]
            products.extend([p for p in prods if p])
            if args and args.verbose > 5:
                print(f"Products for 1c reaction {rxid} {prods}")
        elif numTemplates == 2:
            for pair in get2cPairs(newMatching, oldMatching, rxpos):
                sbses = [substrates[pos] for pos in pair]
                prod = perform2cRx(rxn, sbses, rxid, args)
                if prod:
                    if args and args.verbose > 5:
                        print(f"product of 2c reaction {rxid} {prod}")
                    products.append(prod)
        else:
            print('reactions with 3 or more reactants are impossible at mechanistic level')
            raise NotImplementedError
    return products


def perform1cRx(rxn, sbs, rxid, args):
    mol = Chem.MolFromSmiles(sbs)
    prods = rxn.RunReactants([mol, ])
    if not prods:
        return None
    _ = [Chem.SanitizeMol(prds[0]) for prds in prods]
    prods = [Chem.MolToSmiles(prds[0]) for prds in prods]
    uniqProds = list(set(prods))
    if len(uniqProds) != 1:
        return None
    prod = uniqProds[0]
    if not isProductValid(prod, args):
        return None
    result = (prod, (sbs, ), rxid)
    return result


def perform2cRx(rxn, sbses, rxid, args):
    mols = [Chem.MolFromSmiles(s) for s in sbses]
    prods = rxn.RunReactants(mols)
    if not prods:
        return None
    try:
        _ = [Chem.SanitizeMol(prds[0]) for prds in prods]
    except Chem.rdchem.AtomValenceException:
        print(f"incorrect product from  {Chem.MolToSmiles(prods[0][0])} from {sbses} in {rxid}")
        return None
    prods = [Chem.MolToSmiles(prds[0]) for prds in prods]
    uniqProds = list(set(prods))
    if len(uniqProds) != 1:
        return None
    prod = uniqProds[0]
    if not isProductValid(prod, args):
        return None
    result = (prod, tuple(sbses), rxid)
    return result


def get2cPairs(newMatching, oldMatching, rxpos):
    if rxpos not in newMatching:
        return None
    pairs = []
    # perform new-old
    if oldMatching and oldMatching.get(rxpos):
        for idx1 in newMatching[rxpos].get(0, []):
            for idx2 in oldMatching[rxpos].get(1, []):
                pairs.append((idx1, idx2))
        for idx2 in newMatching[rxpos].get(1, []):
            for idx1 in oldMatching[rxpos].get(0, []):
                pairs.append((idx1, idx2))
    # perform new-new
    for idx1 in newMatching[rxpos].get(0, []):
        for idx2 in newMatching[rxpos].get(1, []):
            pairs.append((idx1, idx2))
    return pairs
