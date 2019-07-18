import untangle
import pandas as pd
from openeye import oechem
from fragmenter import chemi
from multiprocessing import Pool

# discriptors to add to dataframe. Used to filter drugbank

def n_heavy_atoms(smiles):
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    n = sum([not a.IsHydrogen() for a in mol.GetAtoms()])
    return n

def n_rings(smiles):
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    n_rings, parts = oechem.OEDetermineRingSystems(mol)
    return n_rings

def n_aromatic_rings(smiles):
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    nraromsystems, parts = oechem.OEDetermineAromaticRingSystems(mol)
    return nraromsystems

def largest_ring_size(smiles):
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    n_rings, parts = oechem.OEDetermineRingSystems(mol)
    max_i = max(parts)
    l = 0
    for i in range(1, max_i+1):
        n_ring = parts.count(i)
        if n_ring > l:
            l = n_ring
    return l

def n_rotors(smiles):
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    return sum([b.IsRotor() for b in mol.GetBonds()])

def normalize_wbo(smiles, timeout=15):
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    pool = Pool(processes=1)
    result = pool.apply_async(chemi.get_charges, kwds={'molecule': mol, 'strict_stereo': False,
                                                       'strict_types': False})
    try:
        charged = result.get(timeout=timeout)
    except:
        print('process timed out')
        pool.terminate()
        return None
    wbo = 0
    bonds = 0
    for b in charged.GetBonds():
        if 'WibergBondOrder' in b.GetData():
            wbo += b.GetData('WibergBondOrder')
            bonds +=1
        else:
            return None
    pool.terminate()
    return (wbo/bonds)

def determine_connected_components(smiles):
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    count, parts = oechem.OEDetermineComponents(mol)
    return count



#takes 5 minutes
filename="../full_db_9_7_19.xml" # DrugBank Version 5.1.3 (release date: 2019-04-02) Downloaded on 6-7-2019 
obj=untangle.parse(filename)

drugbank=pd.DataFrame(columns=["drugbank_id","name","cas","smiles", "heavy_atoms", "normalized_wbo", "aromatic_rings",
                                 "rotatable_bonds", "n_heterocycles", "fda_approved"])

i=-1
for drug in obj.drugbank.drug:
    drug_type= str(drug["type"])
    
    # select for small molecule drugs
    if drug_type in ["small molecule", "Small Molecule", "Small molecule"]:
        i=i+1    
        
        #Get drugbank_id
        for id in drug.drugbank_id:
            if str(id["primary"])=="true":
                drugbank.loc[i, "drugbank_id"]=id.cdata
        #Drug name
        drugbank.loc[i,"name"]=drug.name.cdata
        
        #Drug CAS
        drugbank.loc[i, "cas"]=drug.cas_number.cdata
        

        if len(drug.calculated_properties.cdata)==0: #If there is no calculated properties
            continue
        
        else:
            for property in drug.calculated_properties.property:
                if property.kind.cdata == "SMILES":
                    drugbank.loc[i, "smiles"]=property.value.cdata
                    
        approval = [g.cdata for g in drug.groups.group]
        if "approved" in approval:
            drugbank.loc[i, "fda_approved"] = True
        else:
            drugbank.loc[i, "fda_approved"] = False

# Only keep entries that have SMILES
drugbank_smiles = drugbank[pd.notnull(drugbank['smiles'])]
# Write out small molecules
drugbank_smiles.to_csv('drugbank_small_mols.csv')
print(drugbank_smiles.shape)

# Apply filters
drugbank_smiles['connected_components'] = drugbank_smiles['smiles'].apply(determine_connected_components)
drugbank_smiles['largest_ring_size'] = drugbank_smiles['smiles'].apply(largest_ring_size)
drugbank_smiles['heavy_atoms'] = drugbank_smiles['smiles'].apply(n_heavy_atoms)
drugbank_smiles['aromatic_rings'] = drugbank_smiles['smiles'].apply(n_aromatic_rings)
drugbank_smiles['n_rings'] = drugbank_smiles['smiles'].apply(n_rings)
drugbank_smiles['rotatable_bonds'] = drugbank_smiles['smiles'].apply(n_rotors)

filtered_drugbank = drugbank_smiles.loc[(drugbank_smiles['largest_ring_size'] <= 14) & 
                                        (drugbank_smiles['largest_ring_size'] >= 3) &
                                        (drugbank_smiles['rotatable_bonds'] >= 4) & 
                                        (drugbank_smiles['aromatic_rings'] >=1) &
                                        (drugbank_smiles['rotatable_bonds'] <= 10) & 
                                        (drugbank_smiles['fda_approved'] == True) &
                                        (drugbank_smiles['connected_components'] == 1)]

# Calculate WBO on filtered set (takes too much memory and time on full drugbank)
filtered_drugbank['normalized_wbo'] = filtered_drugbank['smiles'].apply(normalize_wbo)

# Only keep molecules that WBO calculations does not fail
df = filtered_drugbank.dropna(subset=['normalized_wbo'])
sorted_wbo = df.sort_values(by=['normalized_wbo'], ascending=False)
sorted_wbo.to_csv('drug_bank_filtered.csv')

print('filtered drugbank shape: {}'.format(sorted_wbo.shape))

chemi.to_pdf(drugbank_smiles['smiles'], names=list(drugbank_smiles.name), fname='full_db.pdf')
chemi.to_pdf(sorted_wbo['smiles'], names=list(sorted_wbo.name), fname='filtered_db.pdf')