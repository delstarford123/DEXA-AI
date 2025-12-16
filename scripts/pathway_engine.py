# --- 1. TARGET DEFINITIONS ---
PATHWAY_DB = {
    "NR3C1": {
        "function": "Glucocorticoid Receptor",
        "phenotype": "Anti-inflammatory effect. Risk of hyperglycemia and adrenal suppression."
    },
    "COX-1/2": {
        "function": "Cyclooxygenase (COX)",
        "phenotype": "Analgesia and antipyresis. Risk of gastric mucosal damage and platelet inhibition."
    },
    "ACE": {
        "function": "Angiotensin-Converting Enzyme",
        "phenotype": "Vasodilation (BP reduction). Risk of dry cough and hyperkalemia."
    },
    "OPRM1": {
        "function": "Mu-Opioid Receptor",
        "phenotype": "CNS depression, pain relief. High risk of sedation and respiratory distress."
    },
    "AMPK": {
        "function": "AMP-Activated Protein Kinase",
        "phenotype": "Improved insulin sensitivity. Risk of lactic acidosis."
    },
    "ADRB1": {
        "function": "Beta-Adrenergic Receptor",
        "phenotype": "Reduced cardiac output. Risk of bradycardia and fatigue."
    }
}

# --- 2. INTELLIGENT TARGET RESOLVER ---
def get_target(drug_name):
    """
    Determines the biological target of a drug using:
    1. Direct Lookup (Fast)
    2. Suffix Analysis (Smart Fallback for new drugs)
    """
    name = drug_name.lower().strip()
    
    # A. Explicit Overrides (Common exceptions)
    if name in ['aspirin', 'ibuprofen', 'naproxen', 'diclofenac', 'celecoxib']: return "COX-1/2"
    if name in ['metformin']: return "AMPK"
    if name in ['morphine', 'codeine', 'fentanyl', 'oxycodone']: return "OPRM1"
    
    # B. Suffix/Substring Matching (Handles the 200+ list automatically)
    # Steroids
    if any(x in name for x in ['sone', 'solone', 'sonide', 'cort']): return "NR3C1"
    # NSAIDs
    if any(x in name for x in ['profen', 'coxib', 'lac']): return "COX-1/2"
    # ACE Inhibitors
    if 'pril' in name: return "ACE"
    # Beta Blockers
    if 'olol' in name: return "ADRB1"
    # Opioids
    if any(x in name for x in ['codone', 'morphone']): return "OPRM1"
    # Statins (Included for better 'Unknown' context)
    if 'statin' in name: return "HMG-CoA"
    
    return "Unknown"

# --- 3. DDI LOGIC ENGINE ---
def analyze_ddi(drug_a, drug_b):
    target_a = get_target(drug_a)
    target_b = get_target(drug_b)
    
    # --- DEFAULT STRUCTURED RESPONSE ---
    # This ensures the dashboard always gets all necessary keys, preventing the 'undefined' error.
    response = {
        "risk_level": "LOW",
        "title": "No Critical Interaction Detected",
        "mechanism": f"Drugs target separate pathways ({target_a} and {target_b}). Standard monitoring protocols apply.",
        "consequence": "Standard monitoring recommended. Potential for additive side effects is low.",
        "proteins": f"{target_a} | {target_b}",
        "advice": "Proceed with standard care, observe for common side effects."
    }

    # --- INTERACTION RULES ---

    # 1. NSAID + Steroid (The "Stomach Bleed" Risk)
    if (target_a == "COX-1/2" and target_b == "NR3C1") or (target_b == "COX-1/2" and target_a == "NR3C1"):
        response = {
            "risk_level": "HIGH",
            "title": "Synergistic Gastrointestinal Toxicity",
            "mechanism": "Double-hit mechanism on Gastric Mucosa. NSAIDs inhibit COX-1 (reducing protective prostaglandins), while Steroids inhibit tissue repair mechanisms.",
            "consequence": "4-15x increased risk of peptic ulceration and GI bleeding.",
            "proteins": "COX-1 (Inhibited) + NR3C1 (Activated)",
            "advice": "AVOID combination. If essential, add PPI prophylaxis (e.g., Omeprazole). Monitor stool for occult blood."
        }

    # 2. ACE Inhibitor + NSAID (The "Kidney" Case)
    elif (target_a == "ACE" and target_b == "COX-1/2") or (target_b == "ACE" and target_a == "COX-1/2"):
        response = {
            "risk_level": "MODERATE",
            "title": "Hemodynamic Renal Impairment",
            "mechanism": "Alteration of Glomerular filtration pressure. NSAIDs constrict Afferent Arteriole, while ACE inhibitors dilate Efferent Arteriole.",
            "consequence": "Drop in Glomerular Filtration Rate (GFR) leading to potential Acute Kidney Injury (AKI).",
            "proteins": "ACE (Inhibited) + COX-2 (Inhibited)",
            "advice": "Monitor Serum Creatinine and Potassium. Ensure adequate hydration. Avoid in elderly or hypovolemic patients."
        }

    # 3. Opioid + Depressant (The "Breathing Stop" Risk)
    elif (target_a == "OPRM1" and target_b == "OPRM1") or \
          ("epam" in drug_a.lower() and target_b == "OPRM1") or \
          ("epam" in drug_b.lower() and target_a == "OPRM1"):
        response = {
            "risk_level": "CRITICAL",
            "title": "Additive CNS & Respiratory Depression",
            "mechanism": "Synergistic inhibition of brainstem respiratory centers, increasing sedation.",
            "consequence": "Profound sedation, respiratory arrest, coma, death.",
            "proteins": "Mu-Opioid Receptor (Agonism) + often GABA-A (Modulation)",
            "advice": "STRICTLY CONTRAINDICATED in unmonitored settings. If necessary, reduce dosages by 50% and have Naloxone available."
        }

    # 4. Beta Blocker + Diabetic Med (The "Silent Low Sugar" Risk)
    elif (target_a == "ADRB1" and target_b == "AMPK") or (target_b == "ADRB1" and target_a == "AMPK"):
        response = {
            "risk_level": "MODERATE",
            "title": "Metabolic Masking Effect",
            "mechanism": "Beta-blockade blunts the adrenergic warning signs (tremors, palpitations) of low blood sugar.",
            "consequence": "Patient may not feel early symptoms of hypoglycemia, delaying intervention.",
            "proteins": "Beta-1 Receptor (Blocked) + AMPK (Activated)",
            "advice": "Educate patient to rely on non-adrenergic signs (e.g., sweating) and monitor glucose frequently."
        }

    return response

def predict_phenotype(target_key, score):
    """Predicts phenotype for single drug analysis."""
    if target_key not in PATHWAY_DB: target_key = get_target(target_key)
    info = PATHWAY_DB.get(target_key)
    if info:
        return {"clinical_forecast": info['phenotype'], "function": info['function']}
    return {"clinical_forecast": "Mechanism analysis unavailable for this compound class."}