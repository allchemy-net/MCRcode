This repo contains code for chemical network propagation to our paper "Systematic, computational discovery of multicomponent reactions and one-pot sequences."
# Requirements:
- [python](https://www.python.org/) 3.x (tested on 3.9.16)
- [rdkit](https://www.rdkit.org)  (tested on 2022.03.2)
- [networkx](https://networkx.org/) 2.x (tested on 2.5.1)

  The easier way to install dependencies is [using conda](https://www.rdkit.org/docs/Install.html) alternatively you can install python & rdkit using packae manager from your operating system and the install networkx using pip. There is no special hardware requirements (any hardware on which you an install rdkit and python is fine). Installation of dependencies may take from several minutes (modern laptop or desktop computer with 'normal' internet connection) up to few hours (very slow machine with very slow internet connection).

# Installation:
There is no special installation procedure it is enough to download/clone code and run commands in directory with source.

# Example usage:
`$ python chemNetwork.py  --sbs "CC(C)CC=O.C=CCN.[C-]#[N+]C1CCCCC1.O=C(O)c1ccccc1"  --gens 5`
Required options:
- `--sbs`  list of substrates in smiles format separated by dot e.g. `CC(C)CC=O.C=CCN.[C-]#[N+]C1CCCCC1.O=C(O)c1ccccc1`
- `--gens`  number of synthetic generations to perform

Other options:
- `--wrongprods` - location of file with defined motifs (in SMARTS format) not allowed in products (default `wrong_smarts_list.txt`)
- `--rxdb` - location of file with reaction databse (default `reactions.csv`)
- `--verbose` - verbosity level, default 0, increase it to obtain information usefull for debuging

# Output explanation:
```
$ python chemNetwork.py  --sbs "CC(C)CC=O.C=CCN.[C-]#[N+]C1CCCCC1.O=C(O)c1ccccc1"  --gens 5  
calculation results with paths:
--------------------------------

 CC(C)CC=[OH+] has 1 path(s) of length [1]
 path 1, info: {'numberOfComponents': 1, 'catalytic': 0, 'crossreactivity': {}, 'metals': ''};  Steps:
   {'product': 'CC(C)CC=[OH+]', 'substrates': ('CC(C)CC=O',), 'rxid': 'rx1', 'byproducts': [], 'byproductsReactants': ['OS(=O)(=O)[O-]']}

 C=CC[NH+]=CCC(C)C has 1 path(s) of length [2]
 path 1, info: {'numberOfComponents': 2, 'catalytic': 0, 'crossreactivity': {}, 'metals': ''};  Steps:
   {'product': 'CC(C)CC=[OH+]', 'substrates': ('CC(C)CC=O',), 'rxid': 'rx1', 'byproducts': [], 'byproductsReactants': ['OS(=O)(=O)[O-]']}
   {'product': 'C=CC[NH+]=CCC(C)C', 'substrates': ('CC(C)CC=[OH+]', 'C=CCN'), 'rxid': 'rx2', 'byproducts': ['O'], 'byproductsReactants': ['']}

 C=CCNC(C#[N+]C1CCCCC1)CC(C)C has 1 path(s) of length [3]
 path 1, info: {'numberOfComponents': 3, 'catalytic': 0, 'crossreactivity': {}, 'metals': ''};  Steps:
   {'product': 'CC(C)CC=[OH+]', 'substrates': ('CC(C)CC=O',), 'rxid': 'rx1', 'byproducts': [], 'byproductsReactants': ['OS(=O)(=O)[O-]']}
   {'product': 'C=CC[NH+]=CCC(C)C', 'substrates': ('CC(C)CC=[OH+]', 'C=CCN'), 'rxid': 'rx2', 'byproducts': ['O'], 'byproductsReactants': ['']}
   {'product': 'C=CCNC(C#[N+]C1CCCCC1)CC(C)C', 'substrates': ('C=CC[NH+]=CCC(C)C', '[C-]#[N+]C1CCCCC1'), 'rxid': 'rx7', 'byproducts': [], 'byproductsReactants': ['']}

 C=CCNC(CC(C)C)C(=NC1CCCCC1)OC(=O)c1ccccc1 has 1 path(s) of length [4]
 path 1, info: {'numberOfComponents': 4, 'catalytic': 0, 'crossreactivity': {}, 'metals': ''};  Steps:
   {'product': 'CC(C)CC=[OH+]', 'substrates': ('CC(C)CC=O',), 'rxid': 'rx1', 'byproducts': [], 'byproductsReactants': ['OS(=O)(=O)[O-]']}
   {'product': 'C=CC[NH+]=CCC(C)C', 'substrates': ('CC(C)CC=[OH+]', 'C=CCN'), 'rxid': 'rx2', 'byproducts': ['O'], 'byproductsReactants': ['']}
   {'product': 'C=CCNC(C#[N+]C1CCCCC1)CC(C)C', 'substrates': ('C=CC[NH+]=CCC(C)C', '[C-]#[N+]C1CCCCC1'), 'rxid': 'rx7', 'byproducts': [], 'byproductsReactants': ['']}
   {'product': 'C=CCNC(CC(C)C)C(=NC1CCCCC1)OC(=O)c1ccccc1', 'substrates': ('C=CCNC(C#[N+]C1CCCCC1)CC(C)C', 'O=C(O)c1ccccc1'), 'rxid': 'rx8', 'byproducts': [], 'byproductsReactants': ['']}

 C=CCN(C(=O)c1ccccc1)C(CC(C)C)C(=O)NC1CCCCC1 has 1 path(s) of length [5]
 path 1, info: {'numberOfComponents': 4, 'catalytic': 0, 'crossreactivity': {}, 'metals': ''};  Steps:
   {'product': 'CC(C)CC=[OH+]', 'substrates': ('CC(C)CC=O',), 'rxid': 'rx1', 'byproducts': [], 'byproductsReactants': ['OS(=O)(=O)[O-]']}
   {'product': 'C=CC[NH+]=CCC(C)C', 'substrates': ('CC(C)CC=[OH+]', 'C=CCN'), 'rxid': 'rx2', 'byproducts': ['O'], 'byproductsReactants': ['']}
   {'product': 'C=CCNC(C#[N+]C1CCCCC1)CC(C)C', 'substrates': ('C=CC[NH+]=CCC(C)C', '[C-]#[N+]C1CCCCC1'), 'rxid': 'rx7', 'byproducts': [], 'byproductsReactants': ['']}
   {'product': 'C=CCNC(CC(C)C)C(=NC1CCCCC1)OC(=O)c1ccccc1', 'substrates': ('C=CCNC(C#[N+]C1CCCCC1)CC(C)C', 'O=C(O)c1ccccc1'), 'rxid': 'rx8', 'byproducts': [], 'byproductsReactants': ['']}
   {'product': 'C=CCN(C(=O)c1ccccc1)C(CC(C)C)C(=O)NC1CCCCC1', 'substrates': ('C=CCNC(CC(C)C)C(=NC1CCCCC1)OC(=O)c1ccccc1',), 'rxid': 'rx9', 'byproducts': [], 'byproductsReactants': ['']}
```
Meaning of fields in dictionary with information about path:
- numberOfComponents - number of unique substrates used to synthetise the product
- catalytic - number of compounds which act as catalyst for the path (i.e. byproduct of one reaction from the path is substrate in other reaction)
- crossreactivity - compounds which may be problematic in other reactions
- metals - information about metals use in path

Meaning of fields in dictionary with information about particular reaction/step:
- product - main product of the reaction
- substrates - substrates of the reaction
- rxid - id of reaction (as defined in reaction database)
- byprods - byprods geneerated from reaction transform, eg., leaving groups (cf. Figure S33)
- byprodsReactants - byprods generated from reactants used in the reaction, eg. protonated bases (cf. Figure S34)
