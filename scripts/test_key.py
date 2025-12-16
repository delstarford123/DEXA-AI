import os
import requests
from dotenv import load_dotenv

# 1. Load Key
load_dotenv()
API_KEY = os.environ.get("GEMINI_API_KEY")

print(f"--- 🔑 DIAGNOSTIC TEST ---")
if not API_KEY:
    print("❌ FAIL: API Key not found in environment.")
    print("   Make sure you have a file named '.env' (no name, just extension) with GEMINI_API_KEY=your_key")
    exit()
else:
    print(f"✅ Key Found: {API_KEY[:6]}.......")

# 2. Test Connection (Using the safest, most standard model)
model = "gemini-1.5-flash"
url = f"https://generativelanguage.googleapis.com/v1beta/models/{model}:generateContent?key={API_KEY}"

payload = {
    "contents": [{
        "parts": [{"text": "Hello, are you online?"}]
    }]
}

print(f"📡 Connecting to Google AI ({model})...")
try:
    response = requests.post(url, json=payload, timeout=10)
    
    if response.status_code == 200:
        print("✅ SUCCESS! Google AI replied:")
        print(response.json()['candidates'][0]['content']['parts'][0]['text'])
    else:
        print(f"❌ GOOGLE ERROR ({response.status_code}):")
        print(response.text)
        print("\nPossible Fixes:")
        print("1. If error is 400 'User location is not supported': You need a VPN (US/UK) or a paid account.")
        print("2. If error is 403 'API_KEY_INVALID': Your key is wrong. Copy it again from Google AI Studio.")
        
except Exception as e:
    print(f"❌ NETWORK ERROR: {e}")