import pandas as pd
import pubchempy as pcp
import os
import time

# --- CONFIGURATION ---
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_FILE = os.path.join(BASE_DIR, 'data', 'raw_drug_data.csv')

# --- DRUG LIST (150+ Common Compounds) ---
# Grouped by class to ensure dataset diversity for your AI
target_drugs = [
    # --- Glucocorticoids (The "Positive" Class) ---
    "Triamcinolone", "Methylprednisolone", "Budesonide", "Fluticasone", "Mometasone", 
    "Beclomethasone", "Clobetasol", "Desonide", "Difluprednate", "Fluocinonide", 
    "Halobetasol", "Prednicarbate", "Desoximetasone", "Flurandrenolide", "Alclometasone",
    "Flunisolide", "Ciclesonide", "Deflazacort", "Cortisone", "Paramethasone",

    # --- NSAIDs (Interaction Risk) ---
    "Naproxen", "Diclofenac", "Indomethacin", "Meloxicam", "Piroxicam", 
    "Ketorolac", "Sulindac", "Etodolac", "Nabumetone", "Oxaprozin", 
    "Celecoxib", "Mefenamic acid", "Flurbiprofen", "Diflunisal", "Salsalate",
    "Fenoprofen", "Tolmetin", "Meclofenamate", "Ketoprofen", "Dexketoprofen",

    # --- Analgesics & CNS ---
    "Acetaminophen", "Tramadol", "Morphine", "Codeine", "Oxycodone", 
    "Hydrocodone", "Fentanyl", "Methadone", "Buprenorphine", "Gabapentin",
    "Pregabalin", "Amitriptyline", "Nortriptyline", "Duloxetine", "Venlafaxine",
    "Fluoxetine", "Sertraline", "Paroxetine", "Escitalopram", "Citalopram",
    "Alprazolam", "Clonazepam", "Diazepam", "Lorazepam", "Zolpidem",

    # --- Antibiotics (Interaction Noise) ---
    "Amoxicillin", "Cephalexin", "Ciprofloxacin", "Levofloxacin", "Azithromycin", 
    "Clarithromycin", "Doxycycline", "Minocycline", "Trimethoprim", "Sulfamethoxazole", 
    "Clindamycin", "Metronidazole", "Vancomycin", "Gentamicin", "Tobramycin",
    "Penicillin V", "Ampicillin", "Cefdinir", "Cefuroxime", "Nitrofurantoin",

    # --- Cardiovascular (Metabolic Context) ---
    "Lisinopril", "Enalapril", "Captopril", "Ramipril", "Losartan", 
    "Valsartan", "Irbesartan", "Candesartan", "Amlodipine", "Nifedipine", 
    "Diltiazem", "Verapamil", "Atenolol", "Metoprolol", "Propranolol", 
    "Carvedilol", "Bisoprolol", "Hydrochlorothiazide", "Furosemide", "Spironolactone",
    "Atorvastatin", "Simvastatin", "Rosuvastatin", "Pravastatin", "Lovastatin",
    "Warfarin", "Clopidogrel", "Digoxin", "Isosorbide mononitrate", "Hydralazine",

    # --- Antidiabetics (Crucial for "Metabolic Dysregulation" prediction) ---
    "Metformin", "Glipizide", "Glyburide", "Glimepiride", "Pioglitazone", 
    "Sitagliptin", "Empagliflozin", "Canagliflozin", "Dapagliflozin", "Liraglutide",

    # --- GI & Allergy ---
    "Omeprazole", "Lansoprazole", "Esomeprazole", "Pantoprazole", "Famotidine", 
    "Ranitidine", "Cetirizine", "Loratadine", "Fexofenadine", "Diphenhydramine",
    "Montelukast", "Albuterol", "Salmeterol", "Fluticasone propionate", "Tiotropium"
]

def fetch_and_save():
    print(f"--- 🧪 Expanding Dataset: Fetching {len(target_drugs)} compounds ---")
    
    # 1. Existing Data (Keep your original 6 entries)
    existing_data = [
        {"Drug_Name": "Dexamethasone", "SMILES": "CC1CC2C3CCC4=CC(=O)C=CC4(C3(C(CC2(C1(C(=O)CO)O)C)O)F)C", "CID": 5743},
        {"Drug_Name": "Prednisone", "SMILES": "CC12CC(=O)C3C(C1CCC2(C(=O)CO)O)CCC4=CC(=O)C=CC34C", "CID": 5865},
        {"Drug_Name": "Hydrocortisone", "SMILES": "CC12CCC(=O)C=C1CCC3C2C(CC4(C3CCC4(C(=O)CO)O)C)O", "CID": 5754},
        {"Drug_Name": "Betamethasone", "SMILES": "CC1CC2C3CCC4=CC(=O)C=CC4(C3(C(CC2(C1(C(=O)CO)O)C)O)F)C", "CID": 9782},
        {"Drug_Name": "Ibuprofen", "SMILES": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "CID": 3672},
        {"Drug_Name": "Aspirin", "SMILES": "CC(=O)OC1=CC=CC=C1C(=O)O", "CID": 2244}
    ]
    
    # 2. Fetch New Data
    new_data = []
    for i, drug in enumerate(target_drugs):
        # Skip if already in existing
        if any(d['Drug_Name'].lower() == drug.lower() for d in existing_data):
            continue
            
        try:
            print(f"[{i+1}/{len(target_drugs)}] Fetching {drug}...", end="\r")
            compounds = pcp.get_compounds(drug, 'name')
            if compounds:
                c = compounds[0]
                new_data.append({
                    "Drug_Name": drug,
                    "SMILES": c.canonical_smiles,
                    "CID": c.cid
                })
            time.sleep(0.2) # Polite delay for PubChem API
        except Exception as e:
            print(f"\n❌ Error fetching {drug}: {e}")

    # 3. Combine & Save
    full_data = existing_data + new_data
    df = pd.DataFrame(full_data)
    
    # Reorder columns to match your preferred format: Drug_Name, SMILES, CID
    df = df[['Drug_Name', 'SMILES', 'CID']]
    
    df.to_csv(OUTPUT_FILE, index=False)
    print(f"\n\n✅ Success! Saved {len(df)} compounds to {OUTPUT_FILE}")

if __name__ == "__main__":
    fetch_and_save()