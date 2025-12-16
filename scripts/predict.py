import os
import pickle
import numpy as np
import sys
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

# --- CONFIGURATION ---
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL_PATH = os.path.join(BASE_DIR, 'models', 'dexagen_model.pkl')

def get_fingerprint(smiles):
    """
    Generates the exact same fingerprint used during training (ECFP4).
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: return None
        # Must match preprocess.py: Radius 2, 2048 bits
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        arr = np.zeros((0,), dtype=np.int8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr.reshape(1, -1) # Reshape for single sample prediction
    except Exception as e:
        print(f"Featurization Error: {e}")
        return None

def predict_interaction(smiles, drug_name="Unknown Drug"):
    """
    Loads the model and predicts the probability of interaction.
    """
    # 1. Load Model
    if not os.path.exists(MODEL_PATH):
        print(f"❌ Model not found at {MODEL_PATH}")
        print("   Please run 'scripts/train.py' first.")
        return

    with open(MODEL_PATH, 'rb') as f:
        model = pickle.load(f)

    # 2. Featurize Input
    features = get_fingerprint(smiles)
    if features is None:
        print("❌ Invalid SMILES string provided.")
        return

    # 3. Predict
    try:
        # Get probability of Class 1 (Interaction)
        # Handle edge case where model only learned one class
        if hasattr(model, "predict_proba"):
            proba = model.predict_proba(features)
            if proba.shape[1] > 1:
                score = proba[0][1]
            else:
                # If model only knows one class, return 1.0 if that class is 1, else 0.0
                score = 1.0 if model.classes_[0] == 1 else 0.0
        else:
            score = float(model.predict(features)[0])

        # 4. Interpret Result
        threshold = 0.5
        outcome = "POSITIVE (Interaction Detected)" if score > threshold else "NEGATIVE (No Interaction)"
        
        print("-" * 40)
        print(f"🧪 Drug Name:   {drug_name}")
        print(f"📊 Confidence:  {score:.4f} ({score*100:.1f}%)")
        print(f"🧠 AI Verdict:  {outcome}")
        print("-" * 40)
        
        return score

    except Exception as e:
        print(f"❌ Prediction Error: {e}")

if __name__ == "__main__":
    print("--- 🔮 DEXA AI: Inference Engine ---")
    
    # Test Case 1: Dexamethasone (Should be POSITIVE)
    # This is a known Glucocorticoid
    dexa_smiles = "CC1CC2C3CCC4=CC(=O)C=CC4(C3(C(CC2(C1(C(=O)CO)O)C)O)F)C"
    predict_interaction(dexa_smiles, "Dexamethasone (Test)")

    # Test Case 2: Aspirin (Should be NEGATIVE)
    # This is an NSAID, structurally different
    aspirin_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
    predict_interaction(aspirin_smiles, "Aspirin (Test)")