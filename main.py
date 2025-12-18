import os
import sys
import pickle
import numpy as np
import requests
import pubchempy as pcp
import re
import base64
import io
from flask import Flask, render_template, request, jsonify
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Draw, rdDistGeom
from dotenv import load_dotenv
from google import genai 

# --- 1. PATH & CONFIGURATION ---
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
env_path = os.path.join(BASE_DIR, '.env')
load_dotenv(env_path)

GEMINI_API_KEY = os.environ.get("GEMINI_API_KEY")
MODEL_DIR = os.path.join(BASE_DIR, 'models')
GEMINI_MODEL_NAME = "gemini-2.0-flash" 

# --- 2. GEMINI CLIENT INITIALIZATION ---
client = None
if GEMINI_API_KEY:
    try:
        client = genai.Client(api_key=GEMINI_API_KEY)
        print("✅ Gemini Client Initialized (Global Scope).")
    except Exception as e:
        print(f"❌ Error initializing Gemini client: {e}")

# --- 3. IMPORT ANALYSIS ENGINES ---
try:
    from scripts.further_detail import generate_clinical_report 
    from scripts.simulation_engine import run_molecular_dynamics
    from scripts.pathway_engine import predict_phenotype, analyze_ddi
    from scripts.video_engine import run_video_pipeline
    from scripts.discovery_engine import analyze_new_compound, find_similar_drugs
except ImportError as e:
    print(f"⚠️ Warning: Advanced engines not found. Using fallbacks. Error: {e}")
    # Fallbacks prevent crash if a file is missing
    def run_molecular_dynamics(s): return {"molecular_weight": "0", "logp": "0"}
    def predict_phenotype(n, s): return {"function": "Unknown", "phenotype": "Data unavailable"}
    def analyze_ddi(a, b): return {'risk_level': 'LOW', 'mechanism': 'Unknown', 'advice': 'Monitor.'}
    def generate_clinical_report(n, s, m, p): return f"Report for {n}."
    def run_video_pipeline(c, a, b, p): return {"error": "Video Engine not available."}
    def analyze_new_compound(s): return {"error": "Discovery Engine missing"}
    def find_similar_drugs(s): return {"error": "Discovery Engine missing"}

app = Flask(__name__, template_folder='templates')

# --- 4. LOAD AI MODELS ---
print("--- 🚀 Initializing DexaGen-AI Server ---")
AI_MODELS = {}
try:
    if os.path.exists(os.path.join(MODEL_DIR, 'model_struct.pkl')):
        with open(os.path.join(MODEL_DIR, 'model_struct.pkl'), 'rb') as f:
            AI_MODELS['struct'] = pickle.load(f)
        print("✅ Advanced Ensemble Models Loaded.")
    elif os.path.exists(os.path.join(MODEL_DIR, 'dexagen_model.pkl')):
        with open(os.path.join(MODEL_DIR, 'dexagen_model.pkl'), 'rb') as f:
            AI_MODELS['legacy'] = pickle.load(f)
        print("✅ Standard Model Loaded.")
    else:
        print("⚠️ No models found. Using heuristic mode.")
except Exception as e:
    print(f"❌ Error loading models: {e}")

# --- 5. HELPER FUNCTIONS ---
def get_fingerprint(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: return None
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        arr = np.zeros((0,), dtype=np.int8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr.reshape(1, -1) 
    except: return None

def resolve_smiles(input_str):
    if not input_str: return None, "Empty"
    # Check if input is already SMILES
    if get_fingerprint(input_str) is not None: return input_str, "SMILES"
    try:
        c = pcp.get_compounds(input_str, 'name')
        if c: return c[0].canonical_smiles, "Name"
    except: pass
    return None, "Invalid"

def clean_text_for_speech(text):
    text = re.sub(r'\*+', '', text) 
    text = re.sub(r'#+', '', text)
    return text.strip()

# --- 6. VISUALIZATION ENGINES (2D & 3D) ---
def generate_3d_sdf(smiles):
    """Generates 3D coordinates for the interactive viewer."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return None
        mol = Chem.AddHs(mol)
        rdDistGeom.EmbedMolecule(mol, rdDistGeom.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol)
        return Chem.MolToMolBlock(mol)
    except: return None

def generate_2d_schematic(mol_a, mol_b=None):
    """Generates a 2D structural diagram as Base64 image."""
    try:
        img = None
        if mol_b:
            img = Draw.MolsToGridImage([mol_a, mol_b], molsPerRow=2, subImgSize=(400, 300), 
                                       legends=["Primary Agent", "Secondary Agent"], returnPNG=False)
        else:
            img = Draw.MolToImage(mol_a, size=(800, 600))
        
        img_byte_arr = io.BytesIO()
        img.save(img_byte_arr, format='PNG')
        img_byte_arr.seek(0)
        return f"data:image/png;base64,{base64.b64encode(img_byte_arr.getvalue()).decode('utf-8')}"
    except: return None

# --- 7. ROUTES ---

@app.route('/')
def home():
    return render_template('dashboard.html')

@app.route('/discovery')
def discovery_page():
    """Renders the dedicated Drug Discovery Studio page."""
    return render_template('discovery.html')

@app.route('/predict', methods=['POST'])
def predict_endpoint():
    """
    Single Drug Analysis Endpoint.
    Serves both the Dashboard and Discovery Studio visuals.
    """
    try:
        data = request.get_json()
        input_str = data.get('smiles', '').strip()
        drug_name = data.get('drug_name', input_str)

        # 1. Resolve Structure
        resolved_smiles, _ = resolve_smiles(input_str)
        if not resolved_smiles: 
            return jsonify({'error': 'Invalid Chemical Structure'}), 400

        # 2. Physics Engine Calculation
        md_results = run_molecular_dynamics(resolved_smiles)
        
        # 3. AI Prediction Score
        score = 0.5
        if 'struct' in AI_MODELS:
            fp = get_fingerprint(resolved_smiles)
            if fp is not None: score = float(AI_MODELS['struct'].predict_proba(fp)[0][1])
        elif 'legacy' in AI_MODELS:
            fp = get_fingerprint(resolved_smiles)
            if fp is not None: score = float(AI_MODELS['legacy'].predict_proba(fp)[0][1])

        # 4. Generate Visuals (Critical for Discovery Studio)
        mol_3d = generate_3d_sdf(resolved_smiles)
        mol_obj = Chem.MolFromSmiles(resolved_smiles)
        diagram_url = generate_2d_schematic(mol_obj)

        # 5. Clinical & Phenotype Analysis
        # CRITICAL FIX: Pass the DRUG NAME to get the correct phenotype from your DB
        phenotype = predict_phenotype(drug_name, score)
        
        percentage = score * 100
        risk = "HIGH" if percentage > 75 else "MODERATE" if percentage > 40 else "LOW"
        
        raw_explanation = generate_clinical_report(drug_name, score, md_results, phenotype)
        clean_explanation = clean_text_for_speech(raw_explanation)

        return jsonify({
            'drug_name': drug_name,
            'smiles': resolved_smiles,
            'probability': float(score),
            'percentage': f"{percentage:.1f}%",
            'risk_level': risk,
            'explanation': clean_explanation,
            'md_data': md_results,
            'mol_3d': mol_3d,           # 3D Data for viewer
            'diagram_url': diagram_url, # 2D Image Data
            'phenotype_data': phenotype # Phenotype data for Discovery Studio
        })

    except Exception as e:
        print(f"❌ Analysis Error: {e}")
        return jsonify({'error': str(e)}), 500
# --- ADD THIS IMPORT AT THE TOP OF main.py ---
import itertools

# --- REPLACE THE ddi_endpoint FUNCTION WITH THIS ---
@app.route('/predict_interaction', methods=['POST']) 
def ddi_endpoint():
    """
    Multi-Drug Interaction Analysis (Polypharmacy).
    Handles pairs (A+B) or lists (A+B+C...).
    """
    try:
        data = request.get_json()
        
        # 1. GET DRUGS: Support both new list format and old pair format
        drug_list = data.get('drugs', [])
        
        # Fallback: If 'drugs' list is empty, check for old 'drug_a'/'drug_b' keys
        if not drug_list:
            d_a = data.get('drug_a') or data.get('drugA')
            d_b = data.get('drug_b') or data.get('drugB')
            if d_a and d_b: 
                drug_list = [d_a, d_b]
        
        # 2. VALIDATION
        # Filter out empty strings
        drug_list = [d for d in drug_list if d and d.strip()]
        
        if len(drug_list) < 2:
            return jsonify({'error': 'At least two valid drugs are required for interaction analysis.'}), 400

        # 3. ANALYZE ALL PAIRS
        interactions = []
        highest_risk_score = 0 
        risk_map = {"LOW": 0, "MODERATE": 1, "HIGH": 2, "CRITICAL": 3}
        risk_labels = {0: "LOW", 1: "MODERATE", 2: "HIGH", 3: "CRITICAL"}
        
        # itertools.combinations creates every unique pair: (A,B), (A,C), (B,C)
        pairs = list(itertools.combinations(drug_list, 2))
        
        primary_sdf = None # To store 3D data of the most significant pair
        
        for d1, d2 in pairs:
            # Analyze this specific pair
            res = analyze_ddi(d1, d2)
            
            # Track if this is the riskiest pair found so far
            current_risk_val = risk_map.get(res['risk_level'], 0)
            
            if current_risk_val >= highest_risk_score:
                highest_risk_score = current_risk_val
                # Generate 3D visuals for this "risky" pair
                try:
                    s1, _ = resolve_smiles(d1)
                    s2, _ = resolve_smiles(d2)
                    if s1 and s2:
                        m1 = Chem.MolFromSmiles(s1)
                        m2 = Chem.MolFromSmiles(s2)
                        comb = Chem.CombineMols(m1, m2)
                        comb = Chem.AddHs(comb)
                        rdDistGeom.EmbedMolecule(comb, rdDistGeom.ETKDGv3())
                        primary_sdf = Chem.MolToMolBlock(comb)
                except:
                    pass

            interactions.append({
                "pair": f"{d1} + {d2}",
                "risk": res['risk_level'],
                "mechanism": res['mechanism'],
                "advice": res['advice']
            })

        # 4. FINAL SUMMARY
        final_risk = risk_labels[highest_risk_score]
        
        if highest_risk_score >= 2:
            speech = f"Alert. {final_risk} risk interactions detected in your prescription list. Review the report immediately."
        else:
            speech = f"Analysis complete for {len(drug_list)} drugs. No critical interactions found."

        return jsonify({
            'drugs_analyzed': drug_list,
            'highest_risk': final_risk,
            'interactions': interactions,
            'mol_3d_combined': primary_sdf,
            'summary_speech': speech
        })

    except Exception as e:
        print(f"❌ DDI Error: {e}")
        return jsonify({'error': str(e)}), 500
@app.route('/ask_ai', methods=['POST'])
def ask_ai_endpoint():
    """DEXA AI Proxy"""
    if not client: return jsonify({'error': 'Dexa Client not initialized'}), 500
    data = request.get_json()
    clean_prompt = data['prompt'] + " Summarize concisely in plain text. Do NOT use markdown."
    
    try:
        response = client.models.generate_content(
            model=GEMINI_MODEL_NAME,
            contents=clean_prompt
        )
        ai_text = response.text if response.text else "No response generated."
        clean_text = clean_text_for_speech(ai_text)
        
        return jsonify({
            "candidates": [{"content": {"parts": [{"text": clean_text}]}}]
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/make_video', methods=['POST'])
def make_video_endpoint():
    """Video Pipeline"""
    if not client: return jsonify({"error": "DEXA Client missing."}), 500
    data = request.get_json()
    drug_a = data.get('drugA', 'Drug A')
    drug_b = data.get('drugB', 'Drug B')
    protein = data.get('protein', 'Target') 

    print(f"🎬 Starting video pipeline...")
    result = run_video_pipeline(client, drug_a, drug_b, protein) 
    
    if result.get("success"):
        return jsonify({
            "message": "Success",
            "video_url": result["video_url"],
            "script_preview": result["script_preview"]
        })
    else:
        return jsonify({"error": result["error"]}), 500

@app.route('/discover', methods=['POST'])
def discover_endpoint():
    """
    Handles Discovery Studio requests:
    1. DB Search (Similarity)
    2. Novel Analysis (QED/Ro5)
    """
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        mode = data.get('mode', 'novel') 
        
        if not smiles: return jsonify({'error': 'No SMILES provided'}), 400

        if mode == 'search':
            results = find_similar_drugs(smiles)
        else:
            results = analyze_new_compound(smiles)
        
        if "error" in results:
            return jsonify({'error': results['error']}), 500

        return jsonify({
            "mode": mode,
            "results": results
        })

    except Exception as e:
        print(f"Server Error in /discover: {e}")
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    print(f"✅ Server running on http://127.0.0.1:5000")
    app.run(debug=True, port=5000)