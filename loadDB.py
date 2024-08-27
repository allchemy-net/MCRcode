from errors import WrongDataError

# this module if for loading and chech and convert data from tab-sepated-values file

# module constant variable
requiredFields = ['reaction name', 'reaction smarts', 'conditions', 'conditions described', 'by-products', 'reference', 'OX/RED', 'W/WF', 'irreversible',
                  'RATE', 'incompatible groups', 'protic/aprotic', 'polar/non-polar', 'temperature', 'rxid']

allowedTemperature = ('VL', 'L', 'rt', 'H', 'VH')
allowedAB = ('SB', 'B', 'WB', 'N', 'WA', 'A', 'SA')
allowedWater = ('W', 'WF')
allowedProtic = ('P', 'AP')


# functions:
def loadDB(fn, sep='\t'):
    lines = [line.strip().split(sep) for line in open(fn)]
    rxes = []
    checkFormat(lines[0])
    for rxline in lines[1:]:
        rxinfo = dict()
        for poz, name in enumerate(lines[0]):
            rxinfo[name] = rxline[poz]
        checkField(rxinfo)
        rxes.append(rxinfo)
    checkRxids(rxes)
    return rxes


def checkField(rxinfo):
    hasError = False
    # check acid-base conditions
    for ab in rxinfo['conditions'].split('.'):
        if ab in allowedAB:
            continue
        print(f'not alowed conditions {ab} in {rxinfo["rxid"]}')
        hasError = True
    # check temperature
    for t in rxinfo['temperature'].split('.'):
        if t in allowedTemperature:
            continue
        print(f'not allowed temperature {t} in {rxinfo["rxid"]}')
        hasError = True
    # check water/water-free conditions
    for wf in rxinfo['W/WF'].split('.'):
        if wf in allowedWater:
            continue
        print(f'not allowed water/water-free condition {wf} in {rxinfo["rxid"]}')
        hasError = True
    # check protic conditions
    for p in rxinfo['protic/aprotic'].split('.'):
        if p in allowedProtic:
            continue
        print(f'not allowed prtoci/aprotic conditon {p} in {rxinfo["rxid"]}')
        hasError = True
    if hasError:
        raise WrongDataError
    return True


def checkRxids(rxes):
    rxids = [rx['rxid'] for rx in rxes]
    if len(rxids) != len(set(rxids)):
        for rxid in set(rxids):
            if rxids.count(rxid) > 1:
                print(f'!!ERROR!! reaction id {rxid} present multiple time in reaction database !!')
        raise WrongDataError
    return True


def checkFormat(fields):
    fields = set(fields)
    reqFields = set(requiredFields)
    ignored = fields - reqFields
    missed = reqFields - fields
    if ignored:
        print(f"!WARNING! following field(s) from database will be ignored: {' '.join(ignored)}")
    if missed:
        print(f"!!ERROR!! following field(s) missed: {' '.join(missed)}")
        raise WrongDataError


if __name__ == "__main__":
    rxes = loadDB('reactions.csv')
    print(f'loaded {len(rxes)} reaction(s)')
