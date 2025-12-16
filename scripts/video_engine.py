import os
import sys
import json
import time
from dotenv import load_dotenv
from google import genai
from google.genai.errors import APIError

# RDKit imports for molecular visualization
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D # Used for generating the molecular image file

# Cloudinary is used for uploading the final video file
import cloudinary
import cloudinary.uploader

# --- CONFIGURATION AND INITIALIZATION ---
# Load environment variables (needed for Cloudinary secrets/URL construction)
load_dotenv(override=True)

# --- CLOUDINARY CREDENTIALS SETUP (SECURELY READS CLOUDINARY_URL) ---
CLOUDINARY_URL = os.environ.get("CLOUDINARY_URL") 
CLOUDINARY_API_SECRET = os.environ.get("CLOUDINARY_API_SECRET") 

# Fallback: If CLOUDINARY_URL is not set, try to construct it from separate variables
if not CLOUDINARY_URL and CLOUDINARY_API_SECRET:
    CLOUD_NAME = "duvn7jgzl"
    API_KEY = "894464824777585"
    CLOUDINARY_URL = f"cloudinary://{API_KEY}:{CLOUDINARY_API_SECRET}@{CLOUD_NAME}"

CLOUDINARY_READY = False
if CLOUDINARY_URL:
    try:
        # Configure Cloudinary using the single environment variable URL
        cloudinary.config(cloudinary_url=CLOUDINARY_URL)
        CLOUDINARY_READY = True
    except Exception:
        pass
# --- Imports needed at the top of scripts/video_engine.py ---
import os
import sys
import json
import time
# ... (other imports)
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import pubchempy as pcp # Added this import to the top of video_engine.py

# ... (Rest of CONFIGURATION and utility functions) ...

def get_smiles_from_name(drug_name):
    """
    Uses PubChemPy to resolve a drug name to its canonical SMILES string.
    This prevents hardcoding and ensures visualization accuracy for any drug in your list.
    """
    try:
        # PubChem is very good, but may take a second
        compounds = pcp.get_compounds(drug_name, 'name')
        if compounds:
            # We must use .canonical_smiles to get the standard, parsable string
            return compounds[0].canonical_smiles
    except Exception:
        # Fallback if PubChem is unavailable or fails
        pass
    return None


def create_molecular_image(drug_a_name, drug_b_name, output_path="temp_complex.png"):
    """
    Generates a static 2D image of the drug complex using PubChem resolution.
    """
    # 1. Resolve SMILES for both drugs dynamically
    smiles_a = get_smiles_from_name(drug_a_name)
    smiles_b = get_smiles_from_name(drug_b_name)
    
    if not smiles_a or not smiles_b:
        print(f"❌ RDKit Error: Could not resolve SMILES for {drug_a_name} or {drug_b_name}.")
        return None

    # Combine with a dot separator to represent a mixture/complex
    smiles_complex = f"{smiles_a}.{smiles_b}"
    
    try:
        mol = Chem.MolFromSmiles(smiles_complex)
        if mol is None:
            print("❌ RDKit Error: Could not parse combined SMILES string.")
            return None 

        # Create drawing canvas (600x400)
        d2d = rdMolDraw2D.MolDraw2DCairo(600, 400)
        
        # Draw options for high visibility
        opts = d2d.drawOptions()
        opts.addAtomIndices = False
        opts.legendFontSize = 20
        
        # Draw the molecules
        d2d.DrawMolecule(mol, legend=f"Molecular Analysis: {drug_a_name} + {drug_b_name}")
        d2d.FinishDrawing()

        # Save to file
        with open(output_path, 'wb') as f:
            f.write(d2d.GetDrawingText())
        
        return output_path
    except Exception as e:
        print(f"❌ Drawing Error: {e}")
        return None
# --- 1. CORE FUNCTION: GENERATE SCRIPT (Requires genai_client) ---
def generate_video_script(genai_client, drug_a, drug_b, protein_target):
    # ... (code for generate_video_script remains the same) ...
    """
    Uses the provided genai_client object to generate a detailed, time-coded script.
    """
    if not genai_client:
        return {"error": "Gemini AI client not initialized in main.py."}
        
    prompt = f"""
    You are an expert scientific animator. Generate a detailed, time-coded script for a 30-second instructional video visualizing the molecular action of {drug_a} and the interaction with {drug_b}.
    
    The video must cover:
    1. Single Drug Structure: {drug_a}
    2. Drug-Drug Interaction: {drug_a} + {drug_b} (Highlight the risk based on the primary mechanism).
    3. Protein Binding: How {drug_a} affects the {protein_target}.
    
    Output the script in a JSON format with a list of "frames".
    
    Example Schema:
    [
        {{
            "time": "00:00-00:05",
            "scene": "...",
            "visuals": "...",
            "narration": "..."
        }}
    ]
    """
    
    try:
        response = genai_client.models.generate_content(
            model='gemini-2.5-flash',
            contents=[prompt],
            config={"response_mime_type": "application/json"}
        )
        
        script_text = response.text.strip()
        
        if script_text.startswith('```json'):
            script_text = script_text[7:].strip()
        if script_text.endswith('```'):
            script_text = script_text[:-3].strip()

        script_json = json.loads(script_text)
        
        title = f"{drug_a}_{drug_b}_VideoScript_{int(time.time())}.json"
        script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), title)
        with open(script_path, 'w') as f:
            json.dump(script_json, f, indent=4)
            
        return {"script_path": title, "narration_preview": script_json[0]['narration']}
        
    except Exception as e:
        return {"error": f"Gemini Script Generation Failed: {e}"}


# --- 2. CORE FUNCTION: UPLOAD SIMULATION (UPDATED) ---
def upload_video_simulation(script_title, drug_a_name, drug_b_name):
    """
    Generates a molecular image and then uploads it to Cloudinary.
    """
    if not CLOUDINARY_READY:
        return {"error": "Cloudinary not configured. Check CLOUDINARY_URL in .env."}
    
    temp_file = "temp_mol_complex.png"
    # Step 1: Generate the image file locally
    image_file_path = create_molecular_image(drug_a_name, drug_b_name, output_path=temp_file)
    
    if image_file_path is None:
        return {"error": "Could not generate molecular visualization image."}
        
    public_id = f"dexagen_visual_{script_title.replace(' ', '_').lower().replace('.json', '')}_{int(time.time())}"
    
    try:
        # Step 2: Upload the generated image file
        upload_result = cloudinary.uploader.upload(
            image_file_path, 
            resource_type="image", 
            public_id=public_id,
            folder="dexagen_ai_assets"
        )
        
        # Step 3: Clean up local file
        if os.path.exists(image_file_path):
             os.remove(image_file_path)

        secure_url = upload_result.get("secure_url")
        return {"url": secure_url, "public_id": public_id}
        
    except Exception as e:
        return {"error": f"Cloudinary Upload Failed: {e}"}

# --- 3. MAIN EXECUTION WRAPPER (Public Function) ---
def run_video_pipeline(genai_client, drug_a, drug_b, protein_target):
    """
    Public function called by main.py. Manages the entire pipeline flow.
    """
    
    # Step 1: Generate the Script
    script_result = generate_video_script(genai_client, drug_a, drug_b, protein_target)
    
    if "error" in script_result:
        return script_result
        
    # Step 2: Upload Visualization (MUST pass drug names)
    title = script_result['script_path']
    # NOTE: Passing drug_a and drug_b so they can be used for visualization!
    upload_result = upload_video_simulation(title, drug_a, drug_b) 
    
    if "error" in upload_result:
        return upload_result
        
    return {
        "success": True,
        "video_url": upload_result['url'], # This is now the URL of the 2D molecule image
        "script_preview": script_result['narration_preview']
    }