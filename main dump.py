import os
import sys
import pickle
import numpy as np
import requests
import pubchempy as pcp
import re
import base64
import io
import itertools
from flask import Flask, render_template, request, jsonify, send_from_directory
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Draw, rdDistGeom
from dotenv import load_dotenv
from google import genai 

# --- 1. PATH & CONFIGURATION ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(BASE_DIR)

# Load Environment Variables
env_path = os.path.join(BASE_DIR, '.env')
load_dotenv(env_path)

GEMINI_API_KEY = os.environ.get("GEMINI_API_KEY")
GEMINI_MODEL_NAME = "gemini-2.0-flash" 

# Define Directories
MODEL_DIR = os.path.join(BASE_DIR, 'models')
LAB_IMAGES_DIR = os.path.join(BASE_DIR, 'data', 'laboratory_data', 'images')

app = Flask(__name__, template_folder='templates')

# --- 2. INITIALIZE GEMINI CLIENT ---
client = None
if GEMINI_API_KEY:
    try:
        client = genai.Client(api_key=GEMINI_API_KEY)
        print("✅ Gemini Client Initialized.")
    except Exception as e:
        print(f"❌ Error initializing Gemini: {e}")

# --- 3. IMPORT VIRTUAL LABORATORY MODULES ---
lab_system = None
try:
    # Import the Virtual Lab Brain
    from laboratory.predict import VirtualLab
    # Import the Librarian for Manuals & Status
    from laboratory import manual 
    
    print("⚗️  Loading Virtual Laboratory System...")
    lab_system = VirtualLab()
    print("✅ Virtual Laboratory Online.")
except ImportError as e:
    print(f"⚠️ Critical: Could not import Virtual Lab modules. Error: {e}")
except Exception as e:
    print(f"❌ Error initializing Lab System: {e}")

# --- 4. IMPORT ANALYSIS ENGINES (Fallback Safe) ---
try:
    from scripts.further_detail import generate_clinical_report 
    from scripts.simulation_engine import run_molecular_dynamics
    from scripts.pathway_engine import predict_phenotype, analyze_ddi
    from scripts.video_engine import run_video_pipeline
    from scripts.discovery_engine import analyze_new_compound, find_similar_drugs
except ImportError as e:
    print(f"⚠️ Warning: Advanced engines missing ({e}). Using fallbacks.")
    # Fallbacks to prevent crashes
    def run_molecular_dynamics(s): return {"molecular_weight": "0", "logp": "0", "tpsa": "0"}
    def predict_phenotype(n, s): return {"function": "Unknown", "phenotype": "Data unavailable"}
    def analyze_ddi(a, b): return {'risk_level': 'LOW', 'mechanism': 'Unknown', 'advice': 'Monitor.'}
    def generate_clinical_report(n, s, m, p): return f"Report for {n}."
    def run_video_pipeline(c, a, b, p): return {"error": "Video Engine not available."}
    def analyze_new_compound(s): return {"error": "Discovery Engine missing"}
    def find_similar_drugs(s): return {"error": "Discovery Engine missing"}

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
    try:
        mol = Chem.MolFromSmiles(input_str)
        if mol: return input_str, "SMILES"
    except: pass
    try:
        c = pcp.get_compounds(input_str, 'name')
        if c: return c[0].canonical_smiles, "Name"
    except: pass
    return None, "Invalid"

def clean_text_for_speech(text):
    text = re.sub(r'\*+', '', text) 
    text = re.sub(r'#+', '', text)
    return text.strip()

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

def get_equipment_image_url(equip_name):
    """Matches a recommended tool name to its image file in the data folder."""
    if not os.path.exists(LAB_IMAGES_DIR): return None
    
    # Robust matching: "Microscope" matches "Microscope.png" or "Digital Microscope.png"
    clean_target = equip_name.lower().replace(" ", "")
    
    for filename in os.listdir(LAB_IMAGES_DIR):
        clean_file = filename.lower().replace(" ", "")
        if clean_target in clean_file:
            return f"/lab_images/{filename}"
            
    return None 

# --- 6. PAGE ROUTES ---

@app.route('/')
def home():
    return render_template('dashboard.html')

@app.route('/discovery')
def discovery_page():
    """Renders the dedicated Drug Discovery Studio page."""
    return render_template('discovery.html')

@app.route('/laboratory')
def laboratory_page():
    """Renders the Royal Virtual Laboratory Interface."""
    return render_template('laboratory.html')

# --- 7. VIRTUAL LAB API ROUTES ---

@app.route('/lab/analyze', methods=['POST'])
def lab_analyze():
    """Endpoint for the Neural Chemist to analyze a single drug."""
    if not lab_system: return jsonify({'error': 'Lab System Offline'}), 503
    
    data = request.get_json()
    smiles = data.get('smiles')
    name = data.get('name', 'Unknown Sample')
    
    # Use the Lab System
    result = lab_system.analyze_chemical(smiles, name)
    return jsonify(result)

@app.route('/lab/react', methods=['POST'])
def lab_react():
    """Endpoint to Simulate a Chemical Reaction."""
    if not lab_system: return jsonify({'error': 'Lab System Offline'}), 503
    
    data = request.get_json()
    smiles_a, _ = resolve_smiles(data.get('drug_a'))
    smiles_b, _ = resolve_smiles(data.get('drug_b'))
    
    if not smiles_a or not smiles_b:
        return jsonify({'error': 'Invalid chemicals provided for reaction.'}), 400
        
    result = lab_system.simulate_reaction(smiles_a, smiles_b)
    return jsonify(result)

# ... imports ...

@app.route('/lab/equipment', methods=['POST'])
def lab_equipment():
    if not lab_system: return jsonify({'error': 'Lab System Offline'}), 503
    
    data = request.get_json()
    task = data.get('task', '').strip()
    
    # --- FIX 1: EXACT MATCH OVERRIDE ---
    # If user types "Microscope", don't ask AI. Just give Microscope.
    # We check if the input matches any key in our manual.
    direct_match = None
    for equip_name in manual.LIBRARY.keys():
        if equip_name.lower() in task.lower():
            direct_match = equip_name
            break
            
    enriched_results = []
    
    if direct_match:
        # User searched for a specific item -> Return it with 100% confidence
        results = [{'name': direct_match, 'confidence': 1.0}]
    else:
        # User described a task -> Ask AI for Top 5
        results = lab_system.recommend_equipment(task, top_n=5)
    
    for item in results:
        name = item['name']
        img_url = get_equipment_image_url(name)
        details = manual.get_equipment_details(name)
        
        enriched_results.append({
            "name": name,
            "image_url": img_url,
            "status": details['status'],
            "confidence": f"{int(item['confidence']*100)}%"
        })
    
    return jsonify({
        'task': task,
        'recommendations': enriched_results
    })


# --- NEW ROUTES FOR BUTTONS (Availability & Manuals) ---

@app.route('/lab/check_availability', methods=['POST'])
def check_avail():
    """Called when 'Locate' button is clicked."""
    data = request.get_json()
    name = data.get('name')
    # Fetch from manual.py
    info = manual.get_equipment_details(name)
    return jsonify({
        "status": info['status'],
        "location": info['location']
    })

@app.route('/lab/view_manual', methods=['POST'])
def view_manual():
    """Called when 'Manual' button is clicked."""
    data = request.get_json()
    name = data.get('name')
    # Fetch from manual.py
    info = manual.get_equipment_details(name)
    return jsonify({
        "manual_text": info['manual'],
        "safety_text": info['safety']
    })

@app.route('/lab_images/<path:filename>')
def serve_lab_image(filename):
    """Serves the actual image files from data/laboratory_data/images"""
    return send_from_directory(LAB_IMAGES_DIR, filename)

# --- 8. CLINICAL & DISCOVERY AI ROUTES ---

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
        
        # 3. AI Prediction Score (Use Lab System if available)
        score = 0.5
        if lab_system:
            lab_res = lab_system.analyze_chemical(resolved_smiles, drug_name)
            if isinstance(lab_res, dict):
                # Simple heuristic: small molecules < 500 Da are better drugs
                score = 0.9 if lab_res['molecular_weight'] < 500 else 0.5

        # 4. Generate Visuals
        mol_3d = generate_3d_sdf(resolved_smiles)
        mol_obj = Chem.MolFromSmiles(resolved_smiles)
        diagram_url = generate_2d_schematic(mol_obj)

        # 5. Clinical & Phenotype Analysis
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

@app.route('/predict_interaction', methods=['POST']) 
def ddi_endpoint():
    """
    Multi-Drug Interaction Analysis (Polypharmacy).
    """
    try:
        data = request.get_json()
        drug_list = data.get('drugs', [])
        
        if not drug_list:
            d_a = data.get('drug_a') or data.get('drugA')
            d_b = data.get('drug_b') or data.get('drugB')
            if d_a and d_b: drug_list = [d_a, d_b]
        
        drug_list = [d for d in drug_list if d and d.strip()]
        if len(drug_list) < 2:
            return jsonify({'error': 'Need 2+ drugs.'}), 400

        interactions = []
        highest_risk_score = 0 
        risk_map = {"LOW": 0, "MODERATE": 1, "HIGH": 2, "CRITICAL": 3}
        risk_labels = {0: "LOW", 1: "MODERATE", 2: "HIGH", 3: "CRITICAL"}
        
        pairs = list(itertools.combinations(drug_list, 2))
        primary_sdf = None 
        
        for d1, d2 in pairs:
            res = analyze_ddi(d1, d2)
            current_risk_val = risk_map.get(res['risk_level'], 0)
            
            if current_risk_val >= highest_risk_score:
                highest_risk_score = current_risk_val
                # Generate 3D visuals for risky pair
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
                except: pass

            interactions.append({
                "pair": f"{d1} + {d2}",
                "risk": res['risk_level'],
                "mechanism": res['mechanism'],
                "advice": res['advice']
            })

        final_risk = risk_labels[highest_risk_score]
        speech = f"Alert. {final_risk} risk interactions detected." if highest_risk_score >= 2 else "No critical interactions found."

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
    clean_prompt = data['prompt'] + " Summarize concisely in plain text."
    try:
        response = client.models.generate_content(
            model=GEMINI_MODEL_NAME, contents=clean_prompt
        )
        return jsonify({"candidates": [{"content": {"parts": [{"text": clean_text_for_speech(response.text)}]}}]})
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
    print(f"✅ DexaGen Server running on http://127.0.0.1:5000")
    app.run(debug=True, port=5000)