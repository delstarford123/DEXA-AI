import pandas as pd
import pubchempy as pcp
import os
import time

# --- CONFIGURATION ---
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_FILE = os.path.join(BASE_DIR, 'data', 'raw_drug_data.csv')

# --- MASSIVE DRUG LIST (200+ Compounds) ---
target_drugs = [
    # --- 1. Glucocorticoids (The "Positive" Class - Target: NR3C1) ---
    "Dexamethasone", "Prednisone", "Hydrocortisone", "Betamethasone", "Triamcinolone", 
    "Methylprednisolone", "Budesonide", "Fluticasone", "Mometasone", "Beclomethasone", 
    "Clobetasol", "Desonide", "Difluprednate", "Fluocinonide", "Halobetasol", 
    "Prednicarbate", "Desoximetasone", "Flurandrenolide", "Alclometasone", "Flunisolide", 
    "Ciclesonide", "Deflazacort", "Cortisone", "Paramethasone", "Loteprednol",

    # --- 2. NSAIDs (Interaction Risk - Target: COX-1/2) ---
    "Ibuprofen", "Aspirin", "Naproxen", "Diclofenac", "Indomethacin", 
    "Meloxicam", "Piroxicam", "Ketorolac", "Sulindac", "Etodolac", 
    "Nabumetone", "Oxaprozin", "Celecoxib", "Mefenamic acid", "Flurbiprofen", 
    "Diflunisal", "Salsalate", "Fenoprofen", "Tolmetin", "Meclofenamate", 
    "Ketoprofen", "Dexketoprofen", "Etoricoxib", "Parecoxib", "Lumiracoxib",

    # --- 3. Analgesics & Opioids ---
    "Acetaminophen", "Tramadol", "Morphine", "Codeine", "Oxycodone", 
    "Hydrocodone", "Fentanyl", "Methadone", "Buprenorphine", "Tapentadol", 
    "Hydromorphone", "Oxymorphone", "Meperidine", "Butorphanol", "Nalbuphine",

    # --- 4. Antibiotics (Training Noise) ---
    "Amoxicillin", "Cephalexin", "Ciprofloxacin", "Levofloxacin", "Azithromycin", 
    "Clarithromycin", "Doxycycline", "Minocycline", "Trimethoprim", "Sulfamethoxazole", 
    "Clindamycin", "Metronidazole", "Vancomycin", "Gentamicin", "Tobramycin", 
    "Penicillin V", "Ampicillin", "Cefdinir", "Cefuroxime", "Nitrofurantoin",
    "Moxifloxacin", "Ofloxacin", "Tetracycline", "Erythromycin", "Linezolid",

    # --- 5. Cardiovascular (BP, Statins, Heart) ---
    "Lisinopril", "Enalapril", "Captopril", "Ramipril", "Benazepril",
    "Losartan", "Valsartan", "Irbesartan", "Candesartan", "Telmisartan",
    "Amlodipine", "Nifedipine", "Diltiazem", "Verapamil", "Felodipine",
    "Atenolol", "Metoprolol", "Propranolol", "Carvedilol", "Bisoprolol",
    "Hydrochlorothiazide", "Furosemide", "Spironolactone", "Torsemide", "Chlorthalidone",
    "Atorvastatin", "Simvastatin", "Rosuvastatin", "Pravastatin", "Lovastatin",
    "Warfarin", "Clopidogrel", "Digoxin", "Isosorbide mononitrate", "Hydralazine",
    "Rivaroxaban", "Apixaban", "Dabigatran", "Ticagrelor", "Prasugrel",

    # --- 6. Antidiabetics (Crucial for Metabolic Phenotype) ---
    "Metformin", "Glipizide", "Glyburide", "Glimepiride", "Pioglitazone", 
    "Sitagliptin", "Empagliflozin", "Canagliflozin", "Dapagliflozin", "Liraglutide",
    "Saxagliptin", "Linagliptin", "Repaglinide", "Nateglinide", "Acarbose",

    # --- 7. Neuro/Psych (CNS Agents) ---
    "Gabapentin", "Pregabalin", "Amitriptyline", "Nortriptyline", "Duloxetine", 
    "Venlafaxine", "Fluoxetine", "Sertraline", "Paroxetine", "Escitalopram", 
    "Citalopram", "Alprazolam", "Clonazepam", "Diazepam", "Lorazepam", 
    "Zolpidem", "Quetiapine", "Olanzapine", "Risperidone", "Aripiprazole",
    "Lamotrigine", "Levetiracetam", "Topiramate", "Valproic acid", "Carbamazepine",

    # --- 8. GI & Allergy (Common Co-meds) ---
    "Omeprazole", "Lansoprazole", "Esomeprazole", "Pantoprazole", "Famotidine", 
    "Ranitidine", "Cetirizine", "Loratadine", "Fexofenadine", "Diphenhydramine",
    "Montelukast", "Ondansetron", "Metoclopramide", "Dicyclomine", "Loperamide",

    # --- 9. Antivirals & Antifungals ---
    "Acyclovir", "Valacyclovir", "Oseltamivir", "Tenofovir", "Emtricitabine",
    "Fluconazole", "Ketoconazole", "Itraconazole", "Voriconazole", "Terbinafine",

    # --- 10. Oncology & Immunosuppressants ---
    "Methotrexate", "Cyclophosphamide", "Tamoxifen", "Anastrozole", "Letrozole",
    "Tacrolimus", "Cyclosporine", "Mycophenolate", "Azathioprine", "Hydroxychloroquine"
]

def fetch_and_save():
    print(f"--- 🧪 Expanding Dataset: Target {len(target_drugs)} compounds ---")
    
    # Ensure directory exists
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    
    unique_data = []
    seen_cids = set()
    
    print(f"Fetching data from PubChem... (This may take 1-2 minutes)")
    
    for i, drug in enumerate(target_drugs):
        try:
            # Progress bar
            percent = int((i / len(target_drugs)) * 100)
            bar = "█" * (percent // 5) + "-" * (20 - (percent // 5))
            print(f"\r[{bar}] {percent}% | Fetching: {drug:<20}", end="")
            
            compounds = pcp.get_compounds(drug, 'name')
            
            if compounds:
                c = compounds[0]
                # Duplicate check based on CID (Chemical ID)
                if c.cid not in seen_cids:
                    unique_data.append({
                        "Drug_Name": drug,
                        "SMILES": c.canonical_smiles,
                        "CID": c.cid
                    })
                    seen_cids.add(c.cid)
            
            # Polite API delay
            time.sleep(0.15) 
            
        except Exception as e:
            # Silent fail for individual drugs to keep script running
            continue

    print(f"\n\n✅ Fetch Complete!")
    
    # Convert to DataFrame
    df = pd.DataFrame(unique_data)
    
    # Filter for valid SMILES strings
    df = df[df['SMILES'].str.len() > 0]
    
    # Save
    df.to_csv(OUTPUT_FILE, index=False)
    print(f"💾 Saved {len(df)} unique compounds to: {OUTPUT_FILE}")
    print(f"📊 Dataset Size: {len(df)} rows x {len(df.columns)} columns")

if __name__ == "__main__":
    fetch_and_save()