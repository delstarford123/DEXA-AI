import os
import csv
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors, QED, Lipinski

# --- 1. SETUP PATHS & LOAD DATABASE ---

# This dynamically finds the path to: your_project/data/raw_drug_data.csv
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV_PATH = os.path.join(BASE_DIR, 'data', 'raw_drug_data.csv')

DRUG_DB = []

def load_database():
    """
    Reads the CSV file from data/raw_drug_data.csv and generates RDKit fingerprints.
    """
    global DRUG_DB
    
    if not os.path.exists(CSV_PATH):
        print(f"⚠️ Warning: Database file not found at {CSV_PATH}")
        return

    print(f"⏳ Loading Drug Database from {CSV_PATH}...")
    
    loaded_count = 0
    try:
        with open(CSV_PATH, 'r', encoding='utf-8') as f:
            # Using DictReader to handle columns safely
            reader = csv.DictReader(f)
            
            for row in reader:
                # Flexible column name handling (strips whitespace)
                name = row.get('Drug_Name', '').strip()
                smiles = row.get('SMILES', '').strip()
                cid = row.get('CID', '').strip()

                if name and smiles:
                    try:
                        mol = Chem.MolFromSmiles(smiles)
                        if mol:
                            # Pre-calculate fingerprint for fast searching
                            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                            DRUG_DB.append({
                                "name": name,
                                "smiles": smiles,
                                "cid": cid,
                                "mol": mol,
                                "fp": fp
                            })
                            loaded_count += 1
                    except:
                        # Skip invalid SMILES silently
                        continue
                        
        print(f"✅ Successfully loaded {loaded_count} molecules into memory.")
        
    except Exception as e:
        print(f"❌ Error loading database: {e}")

# Initialize the database immediately when this script is imported
load_database()


# --- 2. DISCOVERY ENGINE LOGIC ---

def analyze_new_compound(smiles):
    """
    Analyzes an unseen compound for drug-likeness (QED, Rule of 5).
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return {"error": "Invalid Chemical Structure"}

        # Calculate quantitative properties
        qed_score = QED.qed(mol) 
        mol_wt = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        
        # Rule of 5 Violations
        ro5_violations = 0
        if mol_wt > 500: ro5_violations += 1
        if logp > 5: ro5_violations += 1
        if Lipinski.NumHDonors(mol) > 5: ro5_violations += 1
        if Lipinski.NumHAcceptors(mol) > 10: ro5_violations += 1

        # Determine if it's a good candidate
        viability = "High" if ro5_violations == 0 and qed_score > 0.6 else "Moderate" if ro5_violations <= 1 else "Low"

        return {
            "type": "novel_analysis",
            "is_viable_candidate": qed_score > 0.5 and ro5_violations <= 1,
            "drug_likeness_score": round(qed_score, 3),
            "rule_of_five_violations": ro5_violations,
            "viability_status": viability,
            "molecular_weight": round(mol_wt, 2),
            "logp": round(logp, 2)
        }
    except Exception as e:
        return {"error": str(e)}

def find_similar_drugs(input_smiles, threshold=0.3):
    """
    Similarity Search: Finds drugs in the CSV DB with similar structures.
    Uses Tanimoto Similarity on Morgan Fingerprints.
    """
    try:
        query_mol = Chem.MolFromSmiles(input_smiles)
        if not query_mol: return {"error": "Invalid SMILES provided"}
        
        # Generate fingerprint for the input drug
        query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=2048)
        
        matches = []
        for drug in DRUG_DB:
            # Compare input vs database drug
            similarity = DataStructs.TanimotoSimilarity(query_fp, drug['fp'])
            
            if similarity > threshold:
                matches.append({
                    "name": drug['name'],
                    "cid": drug['cid'],
                    "smiles": drug['smiles'],
                    "similarity": round(similarity * 100, 1) # Convert to percentage
                })
        
        # Sort by highest similarity
        matches.sort(key=lambda x: x['similarity'], reverse=True)
        
        # Determine if this exact drug exists in our DB
        is_known = any(m['similarity'] == 100.0 for m in matches)
        
        return {
            "type": "similarity_search",
            "matches": matches[:6], # Return top 6 matches
            "count": len(matches),
            "is_known_drug": is_known,
            "best_match": matches[0]['name'] if matches else "Unknown"
        }
    except Exception as e:
        return {"error": str(e)}