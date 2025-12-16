import os
import sys
import pickle
import numpy as np
import requests
import pubchempy as pcp
import re
from flask import Flask, render_template, request, jsonify
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from dotenv import load_dotenv
from google import genai # Import genai for client initialization

# --- 1. PATH & CONFIGURATION ---
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
env_path = os.path.join(BASE_DIR, '.env')
load_dotenv(env_path)

GEMINI_API_KEY = os.environ.get("GEMINI_API_KEY")
MODEL_DIR = os.path.join(BASE_DIR, 'models')
GEMINI_MODEL_NAME = "gemini-2.0-flash" 

# --- 2. GEMINI CLIENT INITIALIZATION (CRITICAL FIX) ---
# Initialize the client in the global scope so all functions (like make_video_endpoint) can access it.
client = None
if GEMINI_API_KEY:
    try:
        client = genai.Client(api_key=GEMINI_API_KEY)
        print("✅ Gemini Client Initialized (Global Scope).")
    except Exception as e:
        print(f"❌ Error initializing Gemini client: {e}")

# --- 3. IMPORT ENGINES ---
try:
    from scripts.further_detail import generate_clinical_report 
    from scripts.simulation_engine import run_molecular_dynamics
    from scripts.pathway_engine import predict_phenotype, analyze_ddi
    from scripts.video_engine import run_video_pipeline # Import the video pipeline
except ImportError:
    print("⚠️ Warning: One or more engines not found. Check 'scripts/' folder.")
    def run_molecular_dynamics(s): return {}
    def predict_phenotype(t, s): return {}
    def analyze_ddi(a, b): return {'risk_level': 'Unknown', 'mechanism': 'Analysis unavailable'}
    def generate_clinical_report(n, s, m, p): return "Explanation module missing."
    def run_video_pipeline(c, a, b, p): return {"error": "Video Engine not available."} # Dummy for video

app = Flask(__name__, template_folder='templates')

# --- 4. LOAD AI MODELS ---
print("--- 🚀 Initializing DexaGen-AI Server ---")
AI_MODELS = {}
try:
    if os.path.exists(os.path.join(MODEL_DIR, 'model_struct.pkl')):
        with open(os.path.join(MODEL_DIR, 'model_struct.pkl'), 'rb') as f:
            AI_MODELS['struct'] = pickle.load(f)
        if os.path.exists(os.path.join(MODEL_DIR, 'model_sim.pkl')):
             with open(os.path.join(MODEL_DIR, 'model_sim.pkl'), 'rb') as f:
                AI_MODELS['sim'] = pickle.load(f)
        print("✅ Advanced Ensemble Models Loaded.")
    else:
        model_path = os.path.join(MODEL_DIR, 'dexagen_model.pkl')
        with open(model_path, 'rb') as f:
            AI_MODELS['legacy'] = pickle.load(f)
        print("✅ Standard Model Loaded.")
except Exception as e:
    print(f"❌ Error loading models: {e}")
    # Do NOT exit, allow server to run for DDI/Gemini features
    pass

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
    if get_fingerprint(input_str) is not None: return input_str, "SMILES"
    try:
        c = pcp.get_compounds(input_str, 'name')
        if c: return c[0].canonical_smiles, "Name"
    except: pass
    return None, "Invalid"

def clean_text_for_speech(text):
    """Removes Markdown characters (*, #) so the voice is clean"""
    text = re.sub(r'\*+', '', text) # Remove asterisks
    text = re.sub(r'#+', '', text)  # Remove hashes
    return text.strip()

# --- 6. ROUTES ---

@app.route('/')
def home():
    return render_template('dashboard.html')

@app.route('/predict', methods=['POST'])
def predict_endpoint():
    """Single Drug Analysis: Physics + AI + Clinical Report"""
    data = request.get_json()
    input_str = data.get('smiles', '').strip()
    drug_name = data.get('drug_name', 'Unknown')

    resolved_smiles, _ = resolve_smiles(input_str)
    if not resolved_smiles: return jsonify({'error': 'Invalid Structure'}), 400

    # 1. Physics Engine
    md_results = run_molecular_dynamics(resolved_smiles)
    
    # 2. AI Prediction
    score = 0.0
    if 'struct' in AI_MODELS:
        fp = get_fingerprint(resolved_smiles)
        if fp is not None:
            score = AI_MODELS['struct'].predict_proba(fp)[0][1]
    elif 'legacy' in AI_MODELS:
        fp = get_fingerprint(resolved_smiles)
        score = float(AI_MODELS['legacy'].predict_proba(fp)[0][1])

    # 3. Phenotype Engine
    phenotype = predict_phenotype("NR3C1", score)
    percentage = score * 100
    risk = "HIGH" if percentage > 75 else "MODERATE" if percentage > 40 else "LOW"
    
    # 4. Generate Detailed Report
    raw_explanation = generate_clinical_report(drug_name, score, md_results, phenotype)
    clean_explanation = clean_text_for_speech(raw_explanation)

    return jsonify({
        'drug_name': drug_name,
        'smiles': resolved_smiles,
        'probability': float(score),
        'percentage': f"{percentage:.1f}%",
        'risk_level': risk,
        'explanation': clean_explanation,
        'md_data': md_results
    })

@app.route('/predict_ddi', methods=['POST'])
def ddi_endpoint():
    """Drug-Drug Interaction Analysis"""
    data = request.get_json()
    drug_a = data.get('drugA', 'Drug A')
    drug_b = data.get('drugB', 'Drug B')
    
    # Run the advanced DDI engine
    result = analyze_ddi(drug_a, drug_b)
    
    return jsonify({
        'drug_a': drug_a,
        'drug_b': drug_b,
        'details': result 
    })

@app.route('/ask_ai', methods=['POST'])
def ask_ai_endpoint():
    """Gemini AI Proxy"""
    if not client: return jsonify({'error': 'Gemini Client not initialized'}), 500
    data = request.get_json()
    
    clean_prompt = data['prompt'] + " Summarize concisely in plain text. Do NOT use markdown."
    
    url = f"https://generativelanguage.googleapis.com/v1beta/models/{GEMINI_MODEL_NAME}:generateContent?key={GEMINI_API_KEY}"
    try:
        response = requests.post(url, json={ "contents": [{ "parts": [{ "text": clean_prompt }] }] }, timeout=30)
        result = response.json()
        
        try:
            raw = result['candidates'][0]['content']['parts'][0]['text']
            result['candidates'][0]['content']['parts'][0]['text'] = clean_text_for_speech(raw)
        except: pass
            
        return jsonify(result)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/make_video', methods=['POST'])
def make_video_endpoint():
    """Triggers the video script generation and upload simulation."""
    if not client: return jsonify({"error": "Gemini Client not initialized. Cannot generate video script."}), 500

    data = request.get_json()
    drug_a = data.get('drugA', 'Dexamethasone')
    drug_b = data.get('drugB', 'Aspirin')
    protein = data.get('protein', 'Glucocorticoid Receptor (NR3C1)') 

    print(f"🎬 Starting video pipeline for {drug_a} and {drug_b}...")
    
    # CORRECT CALL: Pass the global client object first, then the drug names, and the protein target.
    result = run_video_pipeline(client, drug_a, drug_b, protein) 
    
    if result.get("success"):
        return jsonify({
            "message": "Video script generated and simulated upload successful.",
            "video_url": result["video_url"],
            "script_preview": result["script_preview"]
        })
    else:
        print(f"❌ VIDEO PIPELINE FAILED: {result['error']}")
        return jsonify({"error": result["error"]}), 500


if __name__ == '__main__':
    print(f"✅ Server running on http://127.0.0.1:5000")
    app.run(debug=True, port=5000)