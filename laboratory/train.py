import os
# SILENCE TENSORFLOW WARNINGS (Must be before importing tensorflow)
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

import pickle
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import mean_squared_error, accuracy_score
import joblib

# --- 🏰 PATH CONFIGURATION ---
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL_DIR = os.path.join(BASE_DIR, 'models_laboratory')
CHEM_DB_PATH = os.path.join(MODEL_DIR, 'chemical_database.pkl')
EQUIP_DB_PATH = os.path.join(MODEL_DIR, 'equipment_models.pkl')

# --- 🧠 HYPERPARAMETERS ---
FINGERPRINT_SIZE = 2048  # Must match preprocess.py
BATCH_SIZE = 32
EPOCHS = 50 

def load_knowledge_base():
    """Loads the intelligence artifacts created by preprocess.py"""
    print("📂 Loading Laboratory Knowledge Base...")
    
    if not os.path.exists(CHEM_DB_PATH) or not os.path.exists(EQUIP_DB_PATH):
        raise FileNotFoundError("❌ Preprocessed data not found! Run preprocess.py first.")

    with open(CHEM_DB_PATH, 'rb') as f:
        chem_db = pickle.load(f)
        
    with open(EQUIP_DB_PATH, 'rb') as f:
        equip_db = pickle.load(f)
        
    return chem_db, equip_db

def prepare_chemical_tensors(chem_db):
    """
    Converts the dictionary of chemicals into Mathematical Tensors 
    ready for Deep Learning.
    """
    print("⚗️  Transmuting Chemicals into Tensors...")
    
    X_fingerprints = []
    y_logp = []   # Target 1: Solubility (Hydrophobicity)
    y_molwt = []  # Target 2: Molecular Weight
    
    valid_count = 0
    for name, data in chem_db.items():
        # Ensure we have the vector and the properties
        if 'fingerprint_vector' in data and 'logP' in data:
            # Flatten the vector ensures it's 1D array of 0s and 1s
            fp = data['fingerprint_vector'].flatten() 
            X_fingerprints.append(fp)
            y_logp.append(data['logP'])
            y_molwt.append(data['molecular_weight'])
            valid_count += 1
            
    if valid_count == 0:
        print("❌ Error: No valid chemicals found to train on.")
        return None, None, None

    return (np.array(X_fingerprints), 
            np.array(y_logp), 
            np.array(y_molwt))

def build_neural_chemist_model(input_shape):
    """
    Constructs a Multi-Task Deep Neural Network.
    """
    # Input Layer (The Eye: Sees the 2048-bit Fingerprint)
    inputs = keras.Input(shape=(input_shape,), name="molecular_fingerprint")
    
    # Hidden Layers (The Brain: Processes Chemical Structure)
    x = layers.Dense(512, activation="relu", name="dense_1")(inputs)
    x = layers.Dropout(0.3)(x) 
    x = layers.Dense(256, activation="relu", name="dense_2")(x)
    x = layers.Dropout(0.2)(x)
    x = layers.Dense(128, activation="relu", name="dense_3")(x)

    # Output Branch 1: Predict Solubility (LogP)
    out_logp = layers.Dense(1, name="solubility_prediction")(x)

    # Output Branch 2: Predict Molecular Weight
    out_molwt = layers.Dense(1, name="weight_prediction")(x)

    # Compile the Royal Model
    model = keras.Model(inputs=inputs, outputs=[out_logp, out_molwt], name="DexaGen_Neural_Chemist")
    
    model.compile(
        optimizer=keras.optimizers.Adam(learning_rate=0.001),
        loss={
            "solubility_prediction": "mse",
            "weight_prediction": "mse"
        },
        metrics={
            "solubility_prediction": "mae",
            "weight_prediction": "mae"
        }
    )
    
    return model

def train_equipment_expert(equip_db):
    """
    Trains a Classical Machine Learning model to be the 'Lab Manager'.
    """
    print("\n🔬 Training Equipment Recommendation Engine...")
    
    X = equip_db['tfidf_matrix']
    
    # --- FIX WAS APPLIED HERE ---
    # We must access the dataframe inside the dictionary to get the labels
    y = equip_db['dataframe']['encoded_name']
    
    # Random Forest is excellent for categorical classification
    clf = RandomForestClassifier(n_estimators=100, random_state=42)
    clf.fit(X, y)
    
    acc = clf.score(X, y)
    print(f"   ✅ Equipment Expert Accuracy: {acc*100:.2f}%")
    
    # Save the trained expert
    joblib.dump(clf, os.path.join(MODEL_DIR, 'equipment_classifier.joblib'))
    return clf

def main():
    print("=======================================================")
    print("   🚀  INITIATING DEXAGEN AI TRAINING PROTOCOLS       ")
    print("=======================================================")
    
    try:
        # 1. Load Data
        chem_db, equip_db = load_knowledge_base()
        
        # 2. Train Equipment Model (The "Lab Manager")
        train_equipment_expert(equip_db)
        
        # 3. Train Chemical Model (The "Neural Chemist")
        X, y_logp, y_molwt = prepare_chemical_tensors(chem_db)
        
        if X is None:
            print("❌ Aborting training due to lack of chemical data.")
            return

        print(f"\n🧠 Training Neural Chemist on {len(X)} chemical structures...")
        
        # Split data (Handle small datasets gracefully)
        if len(X) < 5:
             # If very few chemicals, don't split, just train on all (for testing)
             X_train, X_val = X, X
             y_logp_train, y_logp_val = y_logp, y_logp
             y_molwt_train, y_molwt_val = y_molwt, y_molwt
        else:
            X_train, X_val, y_logp_train, y_logp_val, y_molwt_train, y_molwt_val = train_test_split(
                X, y_logp, y_molwt, test_size=0.2, random_state=42
            )
        
        model = build_neural_chemist_model(FINGERPRINT_SIZE)
        
        print("\n⚡ Beginning Deep Learning Epochs...")
        history = model.fit(
            x=X_train,
            y={"solubility_prediction": y_logp_train, "weight_prediction": y_molwt_train},
            validation_data=(
                X_val, 
                {"solubility_prediction": y_logp_val, "weight_prediction": y_molwt_val}
            ),
            epochs=EPOCHS,
            batch_size=BATCH_SIZE,
            verbose=1
        )
        
        # 4. Save the Royal Brain
        model_save_path = os.path.join(MODEL_DIR, 'neural_chemist.keras')
        model.save(model_save_path)
        print(f"\n🏆 TRAINING COMPLETE.")
        print(f"   Saved Neural Chemist to: {model_save_path}")
        print(f"   Saved Equipment Expert to: {os.path.join(MODEL_DIR, 'equipment_classifier.joblib')}")
        print("   The Virtual Laboratory is now ONLINE and INTELLIGENT.")
        
    except Exception as e:
        print(f"\n❌ CRITICAL ERROR DURING TRAINING: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()