import os
import requests
from dotenv import load_dotenv

# Load Environment
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
load_dotenv(os.path.join(BASE_DIR, '.env'))
API_KEY = os.environ.get("GEMINI_API_KEY")

if not API_KEY:
    print("❌ API Key not found. Check .env file.")
    exit()

print(f"🔑 Checking available models for key: {API_KEY[:6]}...")

# List Models Endpoint
url = f"https://generativelanguage.googleapis.com/v1beta/models?key={API_KEY}"

try:
    response = requests.get(url)
    if response.status_code == 200:
        models = response.json().get('models', [])
        print("\n✅ AVAILABLE MODELS:")
        for m in models:
            if 'generateContent' in m['supportedGenerationMethods']:
                print(f"  - {m['name'].replace('models/', '')}")
    else:
        print(f"❌ Error listing models: {response.text}")
except Exception as e:
    print(f"❌ Connection Error: {e}")