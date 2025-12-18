import textwrap
import re

# --- HELPER FUNCTIONS ---

def safe_float(value):
    """Converts varied inputs (strings with units, None) into a clean float."""
    try:
        if value is None: return 0.0
        if isinstance(value, (int, float)): return float(value)
        clean = re.sub(r'[^\d.]', '', str(value))
        return float(clean) if clean else 0.0
    except:
        return 0.0

def get_dosage_guideline(mw):
    """
    (Point 4) Predicts dosage strategy based on molecular weight.
    Heuristic: Smaller molecules often require higher frequency or specific dosing.
    """
    if mw < 300:
        return "Standard oral dosing (e.g., 50-100mg) every 4-6 hours due to likely rapid metabolism."
    elif mw < 500:
        return "Once or twice daily oral dosing (e.g., 10-20mg). Optimal bioavailability range."
    else:
        return "High molecular weight suggests potential injection or controlled-release formulation required."

def get_safe_combinations(target_function):
    """
    (New Feature) Predicts drugs that are likely SAFE to use together.
    """
    if "Glucocorticoid" in target_function:
        return " Proton Pump Inhibitors (e.g., Omeprazole) for stomach protection, and certain Antibiotics (e.g., Amoxicillin)."
    if "Cyclooxygenase" in target_function: # NSAIDs
        return " Acetaminophen (Paracetamol) for synergistic pain relief, and Antihistamines."
    if "Angiotensin" in target_function: # ACE Inhibitors
        return " Calcium Channel Blockers (e.g., Amlodipine) and Thiazide Diuretics."
    if "Opioid" in target_function:
        return " NSAIDs (e.g., Ibuprofen) for multimodal analgesia to reduce opioid requirement."
    if "AMP" in target_function: # Metformin
        return " Sulfonylureas (e.g., Glipizide) or Insulin for enhanced glycemic control."
    
    return " Multivitamins and basic analgesics (consult pharmacist for specific interactions)."

def get_clinical_profile(target_function):
    """
    (Points 6 & 7) Returns 'Who Should Avoid' and 'Important Info'.
    """
    profile = {
        "avoid": "General caution in pregnancy, breastfeeding, and severe liver failure.",
        "storage": "Store at room temperature (20-25°C). Keep away from moisture.",
        "legal": "Prescription Only (Rx)",
        "warnings": "May cause dizziness or drowsiness."
    }

    if "Glucocorticoid" in target_function:
        profile.update({
            "avoid": "Patients with systemic fungal infections, uncontrolled diabetes, or active Tuberculosis.",
            "storage": "Protect from light and heat.",
            "legal": "Prescription / Controlled (depending on potency)",
            "warnings": "Do not stop abruptly (risk of adrenal crisis). May mask signs of infection."
        })
    elif "Cyclooxygenase" in target_function:
        profile.update({
            "avoid": "Patients with active peptic ulcers, bleeding disorders, or severe kidney disease.",
            "legal": "OTC (Low dose) / Rx (High dose)",
            "warnings": "Stop use if signs of GI bleeding (black stools) occur."
        })
    elif "Opioid" in target_function:
        profile.update({
            "avoid": "Patients with respiratory depression, severe asthma, or history of substance abuse.",
            "legal": "Controlled Substance (Schedule II/III)",
            "warnings": "High risk of addiction and respiratory arrest. Avoid alcohol completely."
        })
    elif "Angiotensin" in target_function: # BP Meds
        profile.update({
            "avoid": "Pregnant women (teratogenic risk) and patients with bilateral renal artery stenosis.",
            "warnings": "May cause dry cough or angioedema (swelling). Stop if pregnant."
        })
    
    return profile

# --- MAIN GENERATOR FUNCTION ---

def generate_clinical_report(drug_name, score, md_data, phenotype_data):
    """
    Generates a comprehensive medical report covering Points 1-7.
    """
    
    # 1. Parse Data
    percentage = score * 100
    mw = safe_float(md_data.get('molecular_weight', 0))
    tpsa = safe_float(md_data.get('tpsa', 0))
    target_function = phenotype_data.get('function', 'Unknown Target')
    phenotype_prediction = phenotype_data.get('phenotype', 'No specific data available.')

    # 2. Logic for Pill ID / Discovery (Point 1 & 5)
    # If it fits Lipinski rules, it's a good oral drug candidate
    is_drug_like = (mw < 500 and tpsa < 140)
    discovery_status = "High Drug-Likeness (Oral Candidate)" if is_drug_like else "Low Oral Bioavailability (Likely Injectable/Topical)"
    
    tone = "definitive" if percentage > 75 else "cautious" if percentage > 40 else "negative"
    
    # 3. Get Clinical Details
    safety = get_clinical_profile(target_function)
    safe_combos = get_safe_combinations(target_function)
    dosage = get_dosage_guideline(mw)

    # 4. Construct Report
    report = f"""
CLINICAL INTELLIGENCE REPORT: {drug_name}

1. AI MOLECULAR ANALYSIS
   - Prediction: {percentage:.1f}% Probability of Interaction.
   - Target Identified: {target_function}
   - Predicted Phenotype: {phenotype_prediction}

2. DRUG DISCOVERY & ID (Points 1 & 5)
   - Discovery Classification: {discovery_status}
   - Identification Hint: { "Likely a small tablet/capsule" if mw < 500 else "Likely a large molecule/biologic" }.
   - Molecular Weight: {mw:.2f} g/mol

3. DOSAGE & ADMINISTRATION (Point 4)
   - Guideline: {dosage}

4. SAFETY PROFILE (Point 6 - Who Should Avoid)
   - CONTRAINDICATED IN: {safety['avoid']}

5. IMPORTANT INFORMATION (Point 7)
   - Storage: {safety['storage']}
   - Legal Status: {safety['legal']}
   - Special Warnings: {safety['warnings']}

6. INTERACTION INTELLIGENCE (Point 3)
   - Safe Combinations: Generally safe to co-administer with {safe_combos}
   - Recommendation: Based on the {tone} binding profile, this compound { "is a strong candidate for therapy" if percentage > 40 else "requires structural optimization" }.
    """
    
    return report.strip()