import os
import sys
import json
import time
from dotenv import load_dotenv
from google import genai
from google.genai.errors import APIError

# Cloudinary is used for uploading the final video file
import cloudinary
import cloudinary.uploader
import cloudinary.api 

# --- 1. CONFIGURATION AND INITIALIZATION ---

# Load environment variables (GEMINI_API_KEY and CLOUDINARY_URL)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ENV_PATH = os.path.join(BASE_DIR, '..', '.env') 
load_dotenv(ENV_PATH)

# --- GEMINI CLIENT ---
GEMINI_API_KEY = os.environ.get("GEMINI_API_KEY")

if not GEMINI_API_KEY:
    print("❌ ERROR: GEMINI_API_KEY not found in .env file.")
    sys.exit(1)
    
try:
    client = genai.Client(api_key=GEMINI_API_KEY)
    MODEL = 'gemini-2.5-flash'
    print("✅ Gemini Client Initialized.")
except Exception as e:
    print(f"❌ ERROR: Could not initialize Gemini client: {e}")
    sys.exit(1)


# --- CLOUDINARY CREDENTIALS (READ SECURELY FROM .ENV) ---
# FIX: The standard way is to read the full CLOUDINARY_URL variable
CLOUDINARY_URL = os.environ.get("CLOUDINARY_URL") 

if not CLOUDINARY_URL:
    # We will try to construct it from individual parts as a fallback
    API_SECRET = os.environ.get("CLOUDINARY_API_SECRET")
    if API_SECRET:
        CLOUD_NAME = "duvn7jgzl"
        API_KEY = "894464824777585"
        CLOUDINARY_URL = f"cloudinary://{API_KEY}:{API_SECRET}@{CLOUD_NAME}"
    else:
        print("❌ ERROR: CLOUDINARY_URL or CLOUDINARY_API_SECRET not found in .env. Cannot configure Cloudinary.")
        sys.exit(1)

try:
    # Configure Cloudinary using the single environment variable URL
    cloudinary.config(secure=True, cloud_name=None, api_key=None, api_secret=None) # Clear config first
    cloudinary.config(cloudinary_url=CLOUDINARY_URL)
    print("✅ Cloudinary Configured Securely using CLOUDINARY_URL.")
except Exception as e:
    print(f"❌ Cloudinary Configuration Error: {e}")
    sys.exit(1)


# --- 2. GEMINI: VIDEO SCRIPT GENERATION ---

def generate_video_script(drug_a, drug_b, protein_target):
    """
    Uses Gemini AI to generate a detailed, time-coded script for a 3D animation.
    """
    prompt = f"""
    You are an expert scientific animator. Generate a detailed, time-coded script for a 30-second instructional video visualizing drug action.
    
    The video must cover:
    1. Single Drug Structure: {drug_a}
    2. Drug-Drug Interaction: {drug_a} + {drug_b} (Assume HIGH GI Risk for this pair).
    3. Protein Binding: How {drug_a} affects the {protein_target}.
    
    Output the script in a JSON format with a list of "frames".
    """
    print("⏳ Generating video script via Gemini AI...")
    
    try:
        response = client.models.generate_content(
            model=MODEL,
            contents=[prompt],
            config={"response_mime_type": "application/json"}
        )
        
        # Parse the JSON response text
        script_text = response.text.strip()
        
        if script_text.startswith('```json'):
            script_text = script_text[7:].strip()
        if script_text.endswith('```'):
            script_text = script_text[:-3].strip()

        script_json = json.loads(script_text)
        print("✅ Script generated successfully.")
        return script_json
        
    except APIError as e:
        print(f"❌ Gemini API Error during script generation: {e}")
        return None
    except json.JSONDecodeError as e:
        print(f"❌ Error decoding JSON script from AI: {e}")
        return None


# --- 3. CLOUDINARY: ASSET MANAGEMENT (SIMULATION) ---

def upload_final_video_asset(script_title):
    """
    Simulates uploading the rendered video file using the configured Cloudinary account.
    """
    public_id = f"dexagen_video_{script_title.replace(' ', '_').lower()}_{int(time.time())}"
    
    print(f"\n🚀 Simulating video asset upload to Cloudinary (testing API key)...")
    
    try:
        # Using a placeholder image URL for testing Cloudinary config validity
        upload_result = cloudinary.uploader.upload(
            "https://res.cloudinary.com/demo/image/upload/getting-started/shoes.jpg",
            resource_type="auto", 
            public_id=public_id,
            folder="dexagen_ai_assets"
        )
        
        secure_url = upload_result.get("secure_url")
        print(f"✅ Cloudinary Upload Simulation Successful. Public ID: {public_id}")
        return secure_url
        
    except Exception as e:
        print(f"❌ Cloudinary Upload Failed. Check permissions or network connection.")
        print(f"Error details: {e}")
        return None

# --- 4. MAIN EXECUTION FLOW ---

if __name__ == "__main__":
    
    DRUG_A = "Dexamethasone"
    DRUG_B = "Aspirin"
    PROTEIN_T = "Glucocorticoid Receptor (NR3C1)"
    TITLE = f"{DRUG_A}_DDI_Mechanism"
    
    # 1. Generate the Script
    script = generate_video_script(DRUG_A, DRUG_B, PROTEIN_T)
    
    if script:
        # Save the script locally for the animator/renderer
        script_filename = os.path.join(BASE_DIR, f"script_{TITLE}.json")
        with open(script_filename, 'w') as f:
            json.dump(script, f, indent=4)
        print(f"📝 Script saved to {script_filename}")
        
        # --- (Hypothetical step: The actual 3D rendering process happens here) ---

        # 2. Upload the Final Asset (Simulation)
        final_url = upload_final_video_asset(TITLE)
        
        if final_url:
            print("\n-------------------------------------------------")
            print(f"🎉 DEXAGEN VIDEO PIPELINE COMPLETE.")
            print(f"Generated Video Script for: {DRUG_A} + {DRUG_B}")
            print(f"Simulated Asset URL: {final_url}")
            print("-------------------------------------------------")