from rdkit import Chem
from rdkit.Chem import Descriptors

# Given SMILES strings
smiles_list = ['OCC3OC(OCC2OC(OC(C#N)c1ccccc1)C(O)C(O)C2O)C(O)C(O)C3O',
               'CC(=O)OC1=CC=CC=C1C(=O)O',
               'Oc1ccc(Cl)cc1C(=O)Nc2ccc(cc2Cl)N(=O)=O',
               'CC12CCC(=O)C=C1CCC3C2CCC4(C)C3CCC4(O)C#C',
               'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
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

# Feature calculations and writing for each "SMILES" string
for smiles in smiles_list:
    mol = Chem.MolFromSmiles(smiles)
    num_atoms = mol.GetNumAtoms()
    molwt = Descriptors.MolWt(mol)
    tpsa = Descriptors.TPSA(mol)
    logp = Descriptors.MolLogP(mol)
    hba = Descriptors.NumHAcceptors(mol)
    hbd = Descriptors.NumHDonors(mol)
    
    print('Molecule: {}'.format(smiles))
    print('Num of atoms: {}'.format(num_atoms))
    print('Molecular weight: {:.2f}'.format(molwt))
    print('TPSA: {:.2f}'.format(tpsa))
    print('LogP: {:.2f}'.format(logp))
    print('HBA: {}'.format(hba))
    print('HBD: {}'.format(hbd))
    print('-' * 30)
