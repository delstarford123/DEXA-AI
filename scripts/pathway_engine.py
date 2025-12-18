import re

# --- 1. BIOLOGICAL TARGET DATABASE ---
# Maps molecular targets to their specific clinical effects.
PATHWAY_DB = {
    # STEROIDS (NR3C1)
    "NR3C1": {
        "function": "Glucocorticoid Receptor Agonist",
        "phenotype": "Potent anti-inflammatory & immunosuppressive. Modulates glucose metabolism. Risk of adrenal suppression."
    },
    # NSAIDS (COX)
    "COX-1/2": {
        "function": "Cyclooxygenase Inhibitor (NSAID)",
        "phenotype": "Reduces prostaglandin synthesis. analgesic, anti-pyretic. Risk of GI ulceration and renal strain."
    },
    # OPIOIDS (OPRM1)
    "OPRM1": {
        "function": "Mu-Opioid Receptor Agonist",
        "phenotype": "Central analgesia and sedation. Inhibits nociceptive signaling. High risk of respiratory depression."
    },
    "Acetaminophen": {
        "function": "Central COX-3 / Cannabinoid Modulator",
        "phenotype": "Central antipyretic and analgesic. Minimal peripheral anti-inflammatory effect. Hepatotoxic in overdose."
    },
    # ANTIBIOTICS
    "Cell_Wall": {
        "function": "Bacterial Cell Wall Synthesis Inhibitor",
        "phenotype": "Bactericidal. Disrupts peptidoglycan lattice. Effective against dividing bacteria."
    },
    "Ribosome": {
        "function": "Bacterial Protein Synthesis Inhibitor",
        "phenotype": "Bacteriostatic/Bactericidal. Inhibits 30S or 50S ribosomal subunits. halts bacterial growth."
    },
    "DNA_Gyrase": {
        "function": "DNA Gyrase / Topoisomerase Inhibitor",
        "phenotype": "Bactericidal. Prevents bacterial DNA replication and transcription."
    },
    # CARDIOVASCULAR
    "ACE": {
        "function": "ACE Inhibitor",
        "phenotype": "Vasodilation via RAAS blockade. Reduces afterload. Renoprotective. Risk of dry cough."
    },
    "AGTR1": {
        "function": "Angiotensin II Receptor Blocker (ARB)",
        "phenotype": "Vasodilation. Lowers BP without bradykinin side effects."
    },
    "CACNA1": {
        "function": "Calcium Channel Blocker (CCB)",
        "phenotype": "Relaxes vascular smooth muscle. Vasodilation. Risk of peripheral edema."
    },
    "ADRB1": {
        "function": "Beta-Adrenergic Blocker",
        "phenotype": "Negative inotrope/chronotrope. Slows heart rate and reduces O2 demand."
    },
    "HMG-CoA": {
        "function": "HMG-CoA Reductase Inhibitor (Statin)",
        "phenotype": "Lowers LDL cholesterol via hepatic enzyme inhibition. Stabilizes plaques."
    },
    # METABOLIC
    "AMPK": {
        "function": "AMPK Activator",
        "phenotype": "Suppresses hepatic glucose production. Increases insulin sensitivity. First-line for T2DM."
    },
    "Sulfonylurea": {
        "function": "K-ATP Channel Blocker",
        "phenotype": "Stimulates pancreatic beta cells to secrete insulin. Risk of hypoglycemia."
    }
}

# --- 2. DRUG MAPPING (The Router) ---
# Explicitly maps your drugs to the targets above.
DRUG_TARGET_MAP = {
    # Steroids
    "Dexamethasone": "NR3C1", "Prednisone": "NR3C1", "Hydrocortisone": "NR3C1", "Betamethasone": "NR3C1",
    "Triamcinolone": "NR3C1", "Methylprednisolone": "NR3C1", "Budesonide": "NR3C1", "Fluticasone": "NR3C1",
    "Mometasone": "NR3C1", "Beclomethasone": "NR3C1", "Clobetasol": "NR3C1", "Desonide": "NR3C1",
    
    # NSAIDs
    "Ibuprofen": "COX-1/2", "Aspirin": "COX-1/2", "Naproxen": "COX-1/2", "Diclofenac": "COX-1/2",
    "Meloxicam": "COX-1/2", "Celecoxib": "COX-1/2", "Indomethacin": "COX-1/2", "Ketorolac": "COX-1/2",
    
    # Opioids
    "Morphine": "OPRM1", "Codeine": "OPRM1", "Oxycodone": "OPRM1", "Hydrocodone": "OPRM1",
    "Fentanyl": "OPRM1", "Tramadol": "OPRM1", "Methadone": "OPRM1", "Buprenorphine": "OPRM1",
    
    # Antibiotics
    "Amoxicillin": "Cell_Wall", "Cephalexin": "Cell_Wall", "Vancomycin": "Cell_Wall", "Penicillin": "Cell_Wall",
    "Azithromycin": "Ribosome", "Doxycycline": "Ribosome", "Clindamycin": "Ribosome", "Gentamicin": "Ribosome",
    "Ciprofloxacin": "DNA_Gyrase", "Levofloxacin": "DNA_Gyrase",
    
    # Cardio/Metabolic
    "Lisinopril": "ACE", "Enalapril": "ACE", "Ramipril": "ACE",
    "Losartan": "AGTR1", "Valsartan": "AGTR1",
    "Amlodipine": "CACNA1", "Nifedipine": "CACNA1",
    "Metoprolol": "ADRB1", "Atenolol": "ADRB1", "Propranolol": "ADRB1",
    "Atorvastatin": "HMG-CoA", "Simvastatin": "HMG-CoA", "Rosuvastatin": "HMG-CoA",
    "Metformin": "AMPK", "Glipizide": "Sulfonylurea"
}

# --- 3. INTELLIGENT RESOLVER ---
def get_target(drug_name):
    """
    Determines biological target using Map > Suffixes > Unknown.
    """
    if not drug_name: return "Unknown"
    
    # 1. Clean and check explicit map
    clean_name = drug_name.strip().capitalize()
    
    # Direct lookup (Fastest)
    if clean_name in DRUG_TARGET_MAP:
        return DRUG_TARGET_MAP[clean_name]
    
    # 2. Check if the input IS ALREADY a target key (e.g., user passed "NR3C1")
    if clean_name in PATHWAY_DB:
        return clean_name
        
    # 3. Smart Suffix Matching (Fallback for typos/unknowns)
    name_lower = clean_name.lower()
    
    # Steroids
    if any(x in name_lower for x in ['sone', 'solone', 'sonide', 'cort']): return "NR3C1"
    # NSAIDs
    if any(x in name_lower for x in ['profen', 'coxib', 'fenac']): return "COX-1/2"
    # Antibiotics
    if any(x in name_lower for x in ['cillin', 'ceph', 'penem']): return "Cell_Wall"
    if any(x in name_lower for x in ['floxacin', 'mycin', 'cycline']): return "Ribosome"
    # BP Meds
    if 'pril' in name_lower: return "ACE"
    if 'sartan' in name_lower: return "AGTR1"
    if 'olol' in name_lower: return "ADRB1"
    if 'dipine' in name_lower: return "CACNA1"
    if 'statin' in name_lower: return "HMG-CoA"
    
    return "Unknown"

# --- 4. PREDICTION ENGINE ---
def predict_phenotype(input_name, score):
    """
    Returns the full clinical prediction dictionary.
    Input can be a Drug Name ('Dexamethasone') OR a Target ('NR3C1').
    """
    # 1. Resolve Target
    target = get_target(input_name)
    
    # 2. Get Data from DB
    info = PATHWAY_DB.get(target)
    
    # 3. Construct Response
    if info:
        return {
            "function": info['function'],
            "phenotype": info['phenotype'],
            "clinical_forecast": f"{info['function']}: {info['phenotype']}"
        }
    
    # 4. Fallback if still unknown
    return {
        "function": "Novel / Unclassified Mechanism",
        "phenotype": "Specific clinical profile requires experimental validation or is not in the current database.",
        "clinical_forecast": "Mechanism analysis pending."
    }

# --- 5. DDI LOGIC ---
def analyze_ddi(drug_a, drug_b):
    target_a = get_target(drug_a)
    target_b = get_target(drug_b)
    
    response = {
        "risk_level": "LOW",
        "mechanism": f"Mechanisms ({target_a} + {target_b}) appear compatible.",
        "advice": "Standard clinical monitoring recommended."
    }

    # CRITICAL DDI RULES
    if (target_a == "COX-1/2" and target_b == "NR3C1") or (target_b == "COX-1/2" and target_a == "NR3C1"):
        response = {"risk_level": "HIGH", "mechanism": "Synergistic GI Toxicity (Ulceration Risk).", "advice": "Avoid combination. Add PPI if necessary."}
    
    elif (target_a == "ACE" and target_b == "COX-1/2") or (target_b == "ACE" and target_a == "COX-1/2"):
        response = {"risk_level": "MODERATE", "mechanism": "Renal Hemodynamic Strain.", "advice": "Monitor Serum Creatinine and Potassium."}
        
    elif (target_a == "OPRM1" and target_b == "OPRM1"):
        response = {"risk_level": "CRITICAL", "mechanism": "Additive Respiratory Depression.", "advice": "Strictly contraindicated (Duplicate Therapy)."}

    return response