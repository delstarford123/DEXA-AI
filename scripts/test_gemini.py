import os
from dotenv import load_dotenv
from google import genai

# Load environment variables
load_dotenv(override=True) 

API_KEY = os.environ.get("GEMINI_API_KEY")

if not API_KEY:
    print("Test FAILED: API Key not found.")
    exit()

try:
    client = genai.Client(api_key=API_KEY)
    response = client.models.generate_content(
        model='gemini-2.5-flash',
        contents=['Say "Test Successful."'],
    )
    print("---------------------------------")
    print(f"Client Test Result: {response.text.strip()}")
    print("---------------------------------")
except Exception as e:
    print(f"Test FAILED: Exception during connection: {e}")