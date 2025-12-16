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
# We try to import your advanced scripts. If missing, we use simple fallbacks.
try:
    from scripts.further_detail import generate_clinical_report 
    from scripts.simulation_engine import run_molecular_dynamics
    from scripts.pathway_engine import predict_phenotype, analyze_ddi
    from scripts.video_engine import run_video_pipeline
except ImportError:
    print("⚠️ Warning: Advanced engines not found. Using fallbacks.")
    def run_molecular_dynamics(s): return {"molecular_weight": "Calc...", "logp": "Calc..."}
    def predict_phenotype(t, s): return "Phenotype prediction pending."
    def analyze_ddi(a, b): return {'risk_level': 'MODERATE', 'mechanism': 'Cytochrome P450 inhibition likely.', 'advice': 'Monitor patient closely.'}
    def generate_clinical_report(n, s, m, p): return f"Clinical report for {n}. Risk score: {s:.2f}."
    def run_video_pipeline(c, a, b, p): return {"error": "Video Engine not available."}

app = Flask(__name__, template_folder='templates')

# --- 4. LOAD AI MODELS ---
print("--- 🚀 Initializing DexaGen-AI Server ---")
AI_MODELS = {}
try:
    # Try loading the ensemble models first
    if os.path.exists(os.path.join(MODEL_DIR, 'model_struct.pkl')):
        with open(os.path.join(MODEL_DIR, 'model_struct.pkl'), 'rb') as f:
            AI_MODELS['struct'] = pickle.load(f)
        if os.path.exists(os.path.join(MODEL_DIR, 'model_sim.pkl')):
             with open(os.path.join(MODEL_DIR, 'model_sim.pkl'), 'rb') as f:
                AI_MODELS['sim'] = pickle.load(f)
        print("✅ Advanced Ensemble Models Loaded.")
    # Fallback to single legacy model
    elif os.path.exists(os.path.join(MODEL_DIR, 'dexagen_model.pkl')):
        with open(os.path.join(MODEL_DIR, 'dexagen_model.pkl'), 'rb') as f:
            AI_MODELS['legacy'] = pickle.load(f)
        print("✅ Standard Model Loaded.")
    else:
        print("⚠️ No models found.")
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
    # If it's already a SMILES string
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
            # Draw interactions side-by-side
            img = Draw.MolsToGridImage([mol_a, mol_b], molsPerRow=2, subImgSize=(400, 300), 
                                     legends=["Primary Agent", "Secondary Agent"], returnPNG=False)
        else:
            # Draw single
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

@app.route('/predict', methods=['POST'])
def predict_endpoint():
    """Single Drug Analysis: Physics + AI + Visuals + Report"""
    data = request.get_json()
    input_str = data.get('smiles', '').strip()
    drug_name = data.get('drug_name', input_str)

    resolved_smiles, _ = resolve_smiles(input_str)
    if not resolved_smiles: return jsonify({'error': 'Invalid Structure'}), 400

    # 1. Physics Engine Calculation
    md_results = run_molecular_dynamics(resolved_smiles)
    
    # 2. AI Prediction Score
    score = 0.5
    if 'struct' in AI_MODELS:
        fp = get_fingerprint(resolved_smiles)
        if fp is not None: score = AI_MODELS['struct'].predict_proba(fp)[0][1]
    elif 'legacy' in AI_MODELS:
        fp = get_fingerprint(resolved_smiles)
        if fp is not None: score = float(AI_MODELS['legacy'].predict_proba(fp)[0][1])

    # 3. Visuals Generation
    mol_3d = generate_3d_sdf(resolved_smiles)
    mol_obj = Chem.MolFromSmiles(resolved_smiles)
    diagram_url = generate_2d_schematic(mol_obj)

    # 4. Detailed Report Generation
    phenotype = predict_phenotype("NR3C1", score)
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
        'mol_3d': mol_3d,          # 3D Data
        'diagram_url': diagram_url # 2D Image Data
    })

@app.route('/predict_ddi', methods=['POST'])
@app.route('/predict_interaction', methods=['POST']) 
def ddi_endpoint():
    """Drug-Drug Interaction Analysis"""
    data = request.get_json()
    drug_a = data.get('drug_a') or data.get('drugA', 'Drug A')
    drug_b = data.get('drug_b') or data.get('drugB', 'Drug B')
    
    # 1. Resolve structures
    smiles_a, _ = resolve_smiles(drug_a)
    smiles_b, _ = resolve_smiles(drug_b)
    
    if not smiles_a or not smiles_b:
        return jsonify({'error': 'Could not resolve structures.'}), 400

    # 2. Run DDI Analysis Engine
    result = analyze_ddi(drug_a, drug_b)
    
    # 3. Generate Visuals
    sdf_a = generate_3d_sdf(smiles_a)
    sdf_b = generate_3d_sdf(smiles_b)
    
    mol_a = Chem.MolFromSmiles(smiles_a)
    mol_b = Chem.MolFromSmiles(smiles_b)
    diagram_url = generate_2d_schematic(mol_a, mol_b)
    
    sdf_combined = None
    if sdf_a and sdf_b:
        try:
            combined = Chem.CombineMols(mol_a, mol_b)
            combined = Chem.AddHs(combined)
            rdDistGeom.EmbedMolecule(combined, rdDistGeom.ETKDGv3())
            sdf_combined = Chem.MolToMolBlock(combined)
        except: pass

    # Return everything needed for frontend
    return jsonify({
        'drug_a': drug_a,
        'drug_b': drug_b,
        'details': result,
        'explanation': result.get('mechanism', 'Interaction analyzed.'),
        'mol_3d_a': sdf_a,
        'mol_3d_b': sdf_b,
        'mol_3d_combined': sdf_combined,
        'diagram_url': diagram_url,
        'flow_active': True
    })

@app.route('/ask_ai', methods=['POST'])
def ask_ai_endpoint():
    """Gemini AI Proxy"""
    if not client: return jsonify({'error': 'Gemini Client not initialized'}), 500
    data = request.get_json()
    clean_prompt = data['prompt'] + " Summarize concisely in plain text. Do NOT use markdown."
    
    try:
        # Generate content via the Client
        response = client.models.generate_content(
            model=GEMINI_MODEL_NAME,
            contents=clean_prompt
        )
        # Extract text safely
        ai_text = response.text if response.text else "No response generated."
        clean_text = clean_text_for_speech(ai_text)
        
        # Return format expected by frontend
        return jsonify({
            "candidates": [{"content": {"parts": [{"text": clean_text}]}}]
        })
    except Exception as e:
        print(f"Gemini Error: {e}")
        return jsonify({'error': str(e)}), 500

@app.route('/make_video', methods=['POST'])
def make_video_endpoint():
    """Video Pipeline"""
    if not client: return jsonify({"error": "Gemini Client missing."}), 500
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

if __name__ == '__main__':
    print(f"✅ Server running on http://127.0.0.1:5000")
    app.run(debug=True, port=5000)