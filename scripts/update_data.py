import pandas as pd
import os

# Define paths
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_FILE = os.path.join(BASE_DIR, 'data', 'raw_drug_data.csv')

# Load your current data
df = pd.read_csv(DATA_FILE)

# Add Real Scientific Labels (1 = Interacts with GR, 0 = Does not)
# Glucocorticoids (Dex, Pred, Hydro, Beta) -> 1
# NSAIDs (Ibuprofen, Aspirin) -> 0
interaction_map = {
    'Dexamethasone': 1,
    'Prednisone': 1,
    'Hydrocortisone': 1,
    'Betamethasone': 1,
    'Ibuprofen': 0,
    'Aspirin': 0
}

df['Interaction'] = df['Drug_Name'].map(interaction_map)

# Save back
df.to_csv(DATA_FILE, index=False)
print(f"✅ Added 'Interaction' labels to {DATA_FILE}")
print(df[['Drug_Name', 'Interaction']])