import pandas as pd
import numpy as np
import os
import pickle
import sys

# --- PATH SETUP (Ensures internal imports work) ---
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

# --- IMPORT LOCAL LOGIC ENGINES ---
# We use the pathway_engine's intelligence to classify the drug
# We use the simulation_engine to get physicochemical features
try:
    from simulation_engine import run_molecular_dynamics 
    from pathway_engine import get_target, PATHWAY_DB
except ImportError:
    print("⚠️ CRITICAL: Cannot import necessary engine files from 'scripts/'. Ensure pathway_engine.py and simulation_engine.py exist.")
    sys.exit(1)

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

# --- CONFIGURATION ---
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUT_FILE = os.path.join(BASE_DIR, 'data', 'raw_drug_data.csv')
OUTPUT_DIR = os.path.join(BASE_DIR, 'data', 'processed_phase3')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- DEFINITIVE POSITIVE CLASS LIST (GLUCOCORTICOIDS for Binary Labeling) ---
# This list is used to generate the ground truth 'y' label (1=Steroid, 0=Other)
STEROID_KEYWORDS = [
    "dexamethasone", "prednisone", "hydrocortisone", "betamethasone",
    "triamcinolone", "methylprednisolone", "budesonide", "fluticasone",
    "mometasone", "beclomethasone", "clobetasol", "desonide",
    "difluprednate", "fluocinonide", "halobetasol", "prednicarbate",
    "desoximetasone", "alclometasone", "flunisolide", "ciclesonide",
    "deflazacort", "cortisone", "paramethasone", "loteprednol"
]

def get_structural_fingerprint(smiles):
    """Generates the RDKit Morgan Fingerprint (2048-bit ECFP4)"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return None
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        arr = np.zeros((0,), dtype=np.int8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr
    except Exception:
        return None

def extract_pathway_features(drug_name):
    """
    Uses logic from pathway_engine.py to create numerical pathway features.
    
    Returns: A NumPy array representing binary flags for major drug classes.
    """
    target_key = get_target(drug_name)
    
    # Create Binary Flags for major drug classes (New Features for the AI)
    is_steroid = 1 if target_key == "NR3C1" else 0
    is_nsaid = 1 if target_key == "COX-1/2" else 0
    is_cardio = 1 if target_key == "ACE" or target_key == "ADRB1" else 0
    is_cns = 1 if target_key == "OPRM1" else 0
    is_metabolic = 1 if target_key == "AMPK" else 0
    
    return np.array([is_steroid, is_nsaid, is_cardio, is_cns, is_metabolic], dtype=np.float32)


def preprocess_multiscale():
    print("--- 🧠 DEXA AI: Preprocessing 200+ Compounds for Feature Fusion ---")
    
    if not os.path.exists(INPUT_FILE):
        print(f"❌ Input missing: {INPUT_FILE}. Run fetch_200.py first.")
        return

    df = pd.read_csv(INPUT_FILE)
    print(f"Loaded {len(df)} compounds for processing.")
    
    X_combined = [] # List to hold the final, fused feature vectors
    y = []
    
    positive_count = 0

    for i, row in df.iterrows():
        name = str(row['Drug_Name'])
        smiles = str(row['SMILES'])
        
        # Status Bar
        if i % 20 == 0: print(f"\rProcessing {i}/{len(df)}: {name[:20]}...", end="")

        # 1. Structural Features (The 2048-bit core)
        structural_fp = get_structural_fingerprint(smiles)
        if structural_fp is None: continue

        # 2. Molecular Dynamics/Physicochemical Features (From simulation_engine.py)
        # Assuming run_molecular_dynamics is updated to return float keys (MW, TPSA, LogP)
        md_data = run_molecular_dynamics(smiles) 
        md_features = np.array([
            md_data.get('molecular_weight_float', 0.0), # MW
            md_data.get('tpsa_float', 0.0),              # TPSA
            md_data.get('logp', 0.0)                     # Lipophilicity (LogP)
        ], dtype=np.float32)
        
        # 3. Pathway Features (NEW: Using logic from pathway_engine.py)
        pathway_features = extract_pathway_features(name)
        
        # 4. Feature Fusion: Combine Structural (2048), MD (3), and Pathway (5)
        combined_features = np.concatenate([
            structural_fp.reshape(-1), 
            md_features.reshape(-1),
            pathway_features.reshape(-1)
        ])
        X_combined.append(combined_features)

        # 5. Labeling (Ground Truth: 1=Steroid, 0=Other)
        is_active = 1 if any(s in name.lower() for s in STEROID_KEYWORDS) else 0
        y.append(is_active)
        if is_active: positive_count += 1

    # Final conversion to NumPy arrays
    X_final = np.array(X_combined)
    y_final = np.array(y)

    print(f"\n✅ Processed {len(y_final)} valid samples. Final Feature Dimension: {X_final.shape[1]}")
    print(f"   - Positive (Steroids): {positive_count}")
    print(f"   - Negative (Others):   {len(y_final) - positive_count}")

    # Save Artifacts
    # Saving one combined file is much better for training the final model.
    with open(os.path.join(OUTPUT_DIR, 'X_combined.pkl'), 'wb') as f:
        pickle.dump(X_final, f)
    np.save(os.path.join(OUTPUT_DIR, 'y.npy'), y_final)
    
    print("💾 Saved datasets (X_combined.pkl, y.npy) to data/processed_phase3/")

if __name__ == "__main__":
    preprocess_multiscale() 