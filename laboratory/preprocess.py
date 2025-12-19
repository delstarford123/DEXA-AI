import pandas as pd
import pickle
import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
# New RDKit Generator to silence warnings (Modern Standard)
from rdkit.Chem import rdFingerprintGenerator 
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.preprocessing import LabelEncoder

# --- 🏰 ROYAL PATH CONFIGURATION ---
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, 'data')
LAB_SUB_DIR = os.path.join(DATA_DIR, 'laboratory_data', 'images')
MODEL_DIR = os.path.join(BASE_DIR, 'models_laboratory')

# --- 🧪 CHEMISTRY CONFIGURATION ---
# Initialize the modern Fingerprint Generator once (High Performance)
mfgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

class VirtualLabArchitect:
    def __init__(self):
        self.knowledge_base = {
            'chemicals': {}, 
            'equipment': {}, 
        }
        os.makedirs(MODEL_DIR, exist_ok=True)

    def _find_file(self, filename):
        """
        Smart Search: Looks for the file in multiple likely locations.
        """
        possible_paths = [
            os.path.join(LAB_SUB_DIR, filename), # Check data/laboratory_data/
            os.path.join(DATA_DIR, filename),    # Check data/
            os.path.join(BASE_DIR, filename)     # Check root
        ]
        
        for path in possible_paths:
            if os.path.exists(path):
                return path
        return None

    def _analyze_molecule(self, smiles, name, category="General"):
        """
        Transforms text SMILES into a Digital Molecule Object using Modern RDKit Generators.
        """
        if not isinstance(smiles, str) or len(smiles) < 2:
            return None
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # 1. Generate Mathematical Fingerprint (The Modern Way - No Warnings)
                fp = mfgen.GetFingerprint(mol)
                
                # Convert to Numpy Array for AI Training
                fp_arr = np.zeros((1,))
                Chem.DataStructs.ConvertToNumpyArray(fp, fp_arr)
                
                # 2. Calculate Real-World Properties
                props = {
                    'name': name,
                    'category': category,
                    'smiles': smiles,
                    'molecular_weight': round(Descriptors.MolWt(mol), 3),
                    'logP': round(Descriptors.MolLogP(mol), 3), 
                    'h_donors': Descriptors.NumHDonors(mol),
                    'h_acceptors': Descriptors.NumHAcceptors(mol),
                    'fingerprint_vector': fp_arr
                }
                return props
        except Exception:
            return None
        return None

    def process_chemicals(self):
        print("\n⚗️  Phase 1: Synthesizing Chemical Knowledge Base...")
        
        # 1. Load Drugs
        drug_path = self._find_file('raw_drug_data.csv')
        if drug_path:
            print(f"   📂 Loading Drugs from: {drug_path}")
            # on_bad_lines='skip' prevents crashing if a row is formatted badly
            drugs_df = pd.read_csv(drug_path, on_bad_lines='skip', engine='python')
            for _, row in drugs_df.iterrows():
                name = row.get('Drug_Name', row.get('name', 'Unknown Drug'))
                smiles = row.get('SMILES', row.get('smiles', ''))
                
                data = self._analyze_molecule(smiles, name, "Drug")
                if data:
                    self.knowledge_base['chemicals'][data['name']] = data
        else:
            print(f"   ⚠️ Warning: 'raw_drug_data.csv' could not be found.")

        # 2. Load Reagents
        reagent_path = self._find_file('laboratory_reagents.csv')
        if reagent_path:
            print(f"   📂 Loading Reagents from: {reagent_path}")
            reagents_df = pd.read_csv(reagent_path, on_bad_lines='skip', engine='python')
            for _, row in reagents_df.iterrows():
                name = row.get('Reagent_Name', row.get('name', 'Unknown Reagent'))
                smiles = row.get('SMILES', row.get('structure', ''))
                
                data = self._analyze_molecule(smiles, name, "Reagent")
                if data:
                    self.knowledge_base['chemicals'][data['name']] = data
        else:
            print(f"   ⚠️ Warning: 'laboratory_reagents.csv' could not be found.")

        print(f"   ✅ Successfully indexed {len(self.knowledge_base['chemicals'])} chemical entities.")

    def process_equipment(self):
        print("\n🔬 Phase 2: Calibrating Laboratory Equipment Models...")
        
        equip_path = self._find_file('laboratory_equipments.csv')
        
        if not equip_path:
            print(f"   ❌ CRITICAL ERROR: Could not find 'laboratory_equipments.csv'")
            return

        print(f"   📂 Loading Equipment from: {equip_path}")
        
        # FIX IS HERE: on_bad_lines='skip' ignores the rows with too many commas
        try:
            df = pd.read_csv(equip_path, on_bad_lines='skip', engine='python')
        except Exception as e:
            print(f"   ⚠️ Error reading CSV: {e}")
            return
        
        # Check if empty after skipping
        if len(df) == 0:
            print("   ⚠️ Error: The equipment CSV file seems empty or all lines were bad.")
            return

        # Robust Feature Combination
        df['usage'] = df['usage'].fillna('')
        if 'step_of_experiment_used' in df.columns:
            df['step_of_experiment_used'] = df['step_of_experiment_used'].fillna('')
            df['combined_features'] = df['usage'] + " " + df['step_of_experiment_used']
        else:
            df['combined_features'] = df['usage']
        
        # NLP Vectorization
        tfidf = TfidfVectorizer(stop_words='english', max_features=500)
        tfidf_matrix = tfidf.fit_transform(df['combined_features'])
        
        label_enc = LabelEncoder()
        df['encoded_name'] = label_enc.fit_transform(df['name'])
        
        self.knowledge_base['equipment'] = {
            'dataframe': df,
            'tfidf_model': tfidf,
            'tfidf_matrix': tfidf_matrix,
            'label_encoder': label_enc
        }
        print(f"   ✅ Equipment calibration complete. ({len(df)} items loaded)")

    def save_artifacts(self):
        print("\n💾 Phase 3: Archiving Virtual Laboratory State...")
        
        with open(os.path.join(MODEL_DIR, 'chemical_database.pkl'), 'wb') as f:
            pickle.dump(self.knowledge_base['chemicals'], f)
            
        with open(os.path.join(MODEL_DIR, 'equipment_models.pkl'), 'wb') as f:
            pickle.dump(self.knowledge_base['equipment'], f)
            
        print(f"   🌟 SUCCESS: Virtual Laboratory Environment Ready in '{MODEL_DIR}'")

if __name__ == "__main__":
    print("=======================================================")
    print("   🏥  DEXAGEN B1: VIRTUAL LABORATORY PRE-PROCESSOR    ")
    print("=======================================================")
    
    lab = VirtualLabArchitect()
    lab.process_chemicals()
    lab.process_equipment()
    lab.save_artifacts()