from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt

# Data set: Drug names and SMILES strings
drug_names = ['Amigdalin', 'Triamcinolone', 'Niclosamide', 'Ethisterone', 'Dapsone', 'Dialifor', 'Talbutal', 'Raffinose',
              'Bensulide', 'Malathion', 'Hydrocortisone', 'Metolazone', 'Mefenacet', 'Linuron', 'Fluvalinate', 'Nimetazepam']

smiles_list = ['OCC3OC(OCC2OC(OC(C#N)c1ccccc1)C(O)C(O)C2O)C(O)C(O)C3O ',
               'CC34CC(O)C1(F)C(CCC2=CC(=O)C=CC12C)C3CC(O)C4(O)C(=O)CO',
               'Oc1ccc(Cl)cc1C(=O)Nc2ccc(cc2Cl)N(=O)=O',
               'CC12CCC(=O)C=C1CCC3C2CCC4(C)C3CCC4(O)C#C',
               'Nc1ccc(cc1)S(=O)(=O)c2ccc(N)cc2',
               'CCOP(=S)(OCC)SC(CCl)N1C(=O)c2ccccc2C1=O',
               'CCC(C)C1(CC=C)C(=O)NC(=O)NC1=O',
               'OCC1OC(CO)(OC2OC(COC3OC(CO)C(O)C(O)C3O)C(O)C(O)C2O)C(O)C1O',
               'CC(C)OP(=S)(OC(C)C)SCCNS(=O)(=O)c1ccccc1',
               'CCOC(=O)CC(SP(=S)(OC)OC)C(=O)OCC',
               'CC12CC(O)C3C(CCC4=CC(=O)CCC34C)C2CCC1(O)C(=O)CO',
               'CC2Nc1cc(Cl)c(cc1C(=O)N2c3ccccc3C)S(N)(=O)=O',
               'CN(C(=O)COc1nc2ccccc2s1)c3ccccc3',
               'CON(C)C(=O)Nc1ccc(Cl)c(Cl)c1',
               'CC(C)C(Nc1ccc(cc1Cl)C(F)(F)F)C(=O)OC(C#N)c2cccc(Oc3ccccc3)c2',
               'CN2C(=O)CN=C(c1ccccc1)c3cc(ccc23)N(=O)=O']


# Create a dictionary matching drug names and SMILES strings
drug_dict = dict(zip(drug_names, smiles_list))

# Visualization of drug molecules
fig, axes = plt.subplots(4, 4, figsize=(12, 12))
for ax, (drug_name, smiles) in zip(axes.ravel(), drug_dict.items()):
    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol, size=(150, 150))
    ax.imshow(img)
    ax.axis('off')
    ax.set_title(drug_name, fontsize=18)

# Setting subheadings
fig.suptitle('Drugs', fontsize=16)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()
