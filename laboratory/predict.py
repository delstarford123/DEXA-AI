import os
import pickle
import joblib
import numpy as np
import tensorflow as tf
from rdkit import Chem
from rdkit.Chem import AllChem

# --- 🏰 PATH CONFIGURATION ---
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL_DIR = os.path.join(BASE_DIR, 'models_laboratory')

# Artifact Paths
CHEM_MODEL_PATH = os.path.join(MODEL_DIR, 'neural_chemist.keras')
EQUIP_MODEL_PATH = os.path.join(MODEL_DIR, 'equipment_classifier.joblib')
CHEM_DB_PATH = os.path.join(MODEL_DIR, 'chemical_database.pkl')
EQUIP_DB_PATH = os.path.join(MODEL_DIR, 'equipment_models.pkl')

# Constants
FINGERPRINT_RADIUS = 2
FINGERPRINT_BITS = 2048

class VirtualLab:
    def __init__(self):
        print("⚡ Initializing Virtual Laboratory Interface...")
        self._load_models()
    
    def _load_models(self):
        """Loads the AI brains and databases."""
        try:
            # 1. Load Neural Chemist (Deep Learning)
            self.neural_chemist = tf.keras.models.load_model(CHEM_MODEL_PATH)
            
            # 2. Load Equipment Expert (Random Forest)
            self.equipment_expert = joblib.load(EQUIP_MODEL_PATH)
            
            # 3. Load Databases
            with open(CHEM_DB_PATH, 'rb') as f:
                self.chem_db = pickle.load(f)
            
            with open(EQUIP_DB_PATH, 'rb') as f:
                self.equip_data = pickle.load(f)
                
            self.tfidf = self.equip_data['tfidf_model']
            self.label_encoder = self.equip_data['label_encoder']
            
            print("   ✅ Systems Online: Neural Chemist & Equipment Expert Ready.")
            
        except Exception as e:
            print(f"   ❌ SYSTEM FAILURE: Could not load models. {e}")
            print("   Did you run train.py?")
            exit()

    def _get_fingerprint(self, smiles):
        """Converts a SMILES string into a 2048-bit digital vector."""
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=FINGERPRINT_RADIUS, nBits=FINGERPRINT_BITS)
            arr = np.zeros((1,))
            AllChem.DataStructs.ConvertToNumpyArray(fp, arr)
            return arr.reshape(1, -1) # Reshape for TensorFlow (1 sample, 2048 features)
        return None

    def analyze_chemical(self, smiles, name="Unknown Sample"):
        """
        Performs a deep analysis of a single chemical.
        """
        vector = self._get_fingerprint(smiles)
        if vector is None:
            return "❌ Invalid Chemical Structure"

        # Ask the Neural Chemist
        predictions = self.neural_chemist.predict(vector, verbose=0)
        pred_logp = float(predictions[0][0][0])
        pred_molwt = float(predictions[1][0][0])

        return {
            "name": name,
            "solubility_logp": pred_logp,
            "molecular_weight": pred_molwt,
            "analysis": self._interpret_logp(pred_logp)
        }

    def simulate_reaction(self, smiles_a, smiles_b):
        """
        👑 THE CROWN JEWEL: REACTION SIMULATION
        Mathematically 'mixes' two chemicals and predicts the properties of the result.
        """
        vec_a = self._get_fingerprint(smiles_a)
        vec_b = self._get_fingerprint(smiles_b)

        if vec_a is None or vec_b is None:
            return "❌ One of the reactants is invalid."

        # Simulate Reaction: We create a 'Hybrid Vector'
        # In deep learning terms, this represents the theoretical transition state or product complex
        reaction_vector = (vec_a + vec_b) / 2.0 

        # Predict properties of this new mixture
        predictions = self.neural_chemist.predict(reaction_vector, verbose=0)
        res_logp = float(predictions[0][0][0])
        res_wt = float(predictions[1][0][0])

        return {
            "reaction_type": "Virtual Synthesis",
            "theoretical_weight": res_wt,
            "theoretical_solubility": res_logp,
            "outcome_description": self._interpret_reaction(res_logp)
        }

    def recommend_equipment(self, experiment_description, top_n=5):
        """
        Reads user's intent and suggests the TOP N tools.
        Updated to show MORE results (Top 5) with a lower threshold.
        """
        # 1. Convert text to numbers
        text_vector = self.tfidf.transform([experiment_description])
        
        # 2. Get probabilities for ALL tools (e.g., 80% Microscope, 15% Pipette, 5% Beaker)
        # predict_proba returns a list of scores for every known class
        probabilities = self.equipment_expert.predict_proba(text_vector)[0]
        
        # 3. Sort them to find the highest scores
        # argsort gives us the indices of the sorted scores. [::-1] reverses it to be High->Low
        top_indices = probabilities.argsort()[-top_n:][::-1]
        
        recommendations = []
        for idx in top_indices:
            score = probabilities[idx]
            
            # LOWER THRESHOLD: Show item if AI is even 5% sure (was 10%)
            # This ensures you see more of "what is actually there"
            if score > 0.05:
                name = self.label_encoder.inverse_transform([idx])[0]
                recommendations.append({
                    "name": name,
                    "confidence": float(score)
                })
        
        return recommendations

    def _interpret_logp(self, logp):
        if logp < 0: return "Highly Water Soluble (Hydrophilic)"
        if logp < 2: return "Moderately Soluble"
        if logp < 4: return "Lipophilic (Fat Soluble)"
        return "Highly Lipophilic (Poor Water Solubility)"

    def _interpret_reaction(self, logp):
        if logp > 3:
            return "The resulting product appears highly stable and organic-soluble. Likely precipitates in water."
        else:
            return "The reaction yields a polar, water-soluble compound. Suitable for aqueous solution analysis."

# --- 🖥️ COMMAND LINE INTERFACE (FOR TESTING) ---
if __name__ == "__main__":
    lab = VirtualLab()
    
    print("\n🧪 --- VIRTUAL LAB TEST: SINGLE CHEMICAL ---")
    # Test with Aspirin (SMILES: CC(=O)Oc1ccccc1C(=O)O)
    result = lab.analyze_chemical("CC(=O)Oc1ccccc1C(=O)O", "Aspirin")
    print(f"Analysis of {result['name']}:")
    print(f"  - Predicted Weight: {result['molecular_weight']:.2f} g/mol")
    print(f"  - Predicted LogP:   {result['solubility_logp']:.2f}")
    print(f"  - AI Insight:       {result['analysis']}")

    print("\n⚗️ --- VIRTUAL LAB TEST: REACTION SIMULATION ---")
    # Mixing Aspirin + Water (O) (Simulation)
    # Note: Water SMILES is 'O'
    reaction = lab.simulate_reaction("CC(=O)Oc1ccccc1C(=O)O", "O")
    print("Reaction: Aspirin + Water")
    print(f"  - Result Weight:    {reaction['theoretical_weight']:.2f}")
    print(f"  - Result LogP:      {reaction['theoretical_solubility']:.2f}")
    print(f"  - AI Conclusion:    {reaction['outcome_description']}")

    print("\n🔬 --- VIRTUAL LAB TEST: EQUIPMENT MANAGER ---")
    user_task = "I need to separate the solid precipitate from the liquid mixture rapidly."
    recommendations = lab.recommend_equipment(user_task)
    print(f"User Task: '{user_task}'")
    print("Top Recommendations:")
    for item in recommendations:
        print(f"  - {item['name']} ({int(item['confidence']*100)}% match)")