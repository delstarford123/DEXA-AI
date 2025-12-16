import textwrap
import re

def safe_float(value):
    """
    Helper: Converts varied inputs (strings with units, None, etc.) into a clean float.
    Ex: "362.47 g/mol" -> 362.47
    """
    try:
        if value is None: return 0.0
        if isinstance(value, (int, float)): return float(value)
        
        # Remove non-numeric characters (except dot)
        clean = re.sub(r'[^\d.]', '', str(value))
        return float(clean) if clean else 0.0
    except:
        return 0.0

def generate_clinical_report(drug_name, score, md_data, phenotype_data):
    """
    Generates a comprehensive medical report based on Multi-Scale Analysis.
    """
    
    # 1. Interpret AI Confidence
    percentage = score * 100
    if percentage > 75:
        affinity_desc = "High Binding Affinity (Strong Agonist Potential)"
        tone = "definitive"
    elif percentage > 40:
        affinity_desc = "Moderate Binding Affinity (Partial Agonist/Modulator)"
        tone = "cautious"
    else:
        affinity_desc = "Negligible Binding Affinity"
        tone = "negative"

    # 2. Analyze Physicochemical Properties (Lipinski's Rule of 5)
    # FIX: Use safe_float to handle strings like "362.47 g/mol"
    mw = safe_float(md_data.get('molecular_weight', 0))
    tpsa = safe_float(md_data.get('tpsa', 0))
    
    absorption_analysis = []
    
    # MW Rule: < 500 Da is good
    if mw > 500:
        absorption_analysis.append(f"⚠️ High Molecular Weight ({mw:.1f} Da). Oral absorption may be limited.")
    else:
        absorption_analysis.append(f"✅ Molecular Weight ({mw:.1f} Da) is optimal for oral bioavailability.")
        
    # TPSA Rule: < 140 A^2 is good
    if tpsa > 140:
        absorption_analysis.append(f"⚠️ High Polarity (TPSA {tpsa:.1f} Å²). Blood-Brain Barrier (BBB) penetration unlikely.")
    else:
        absorption_analysis.append(f"✅ TPSA ({tpsa:.1f} Å²) suggests good intestinal absorption and membrane permeability.")

    physio_summary = " ".join(absorption_analysis)

    # 3. Generate Clinical Protocol based on Target
    # Extracts the "function" from the dictionary safely
    target_function = phenotype_data.get('function', 'Unknown Target')
    
    protocol = []
    
    if "Glucocorticoid" in target_function:
        protocol = [
            "Baseline Monitoring: Assess Fasting Blood Glucose, HbA1c, and Lipid Profile prior to initiation.",
            "Risk Management: Weigh benefits against risk of metabolic dysregulation (Diabetes) and adrenal suppression.",
            "Patient Education: Instruct patient to report symptoms of hyperglycemia (polydipsia, polyuria) immediately.",
            "Long-term Strategy: If chronic use is required, implement a tapering schedule to prevent adrenal crisis."
        ]
    elif "Cyclooxygenase" in target_function:
        protocol = [
            "Gastrointestinal Safety: Assess history of peptic ulcers. Consider prophylaxis with PPIs for high-risk patients.",
            "Renal Function: Monitor Serum Creatinine and GFR, especially in elderly or hypovolemic patients.",
            "Cardiovascular Risk: Use lowest effective dose for shortest duration to minimize thrombotic risk."
        ]
    elif "Angiotensin" in target_function or "Beta" in target_function:
        protocol = [
            "Hemodynamics: Monitor Blood Pressure and Heart Rate regularly.",
            "Electrolytes: Check Potassium levels (risk of Hyperkalemia).",
            "Adherence: Educate patient on the importance of consistent dosing to prevent rebound hypertension."
        ]
    elif "Opioid" in target_function:
        protocol = [
            "Respiratory Safety: Monitor respiratory rate and sedation level.",
            "Addiction Risk: Screen utilizing tools like SOAPP-R prior to initiation.",
            "Co-prescribing: Avoid concurrent Benzodiazepines or CNS depressants."
        ]
    elif "AMP" in target_function: # Metformin/Diabetes
        protocol = [
            "Renal Function: Contraindicated if eGFR < 30 mL/min/1.73 m² due to risk of lactic acidosis.",
            "Vitamin B12: Periodic monitoring recommended as long-term use can decrease absorption.",
            "Gastrointestinal: Titrate dose slowly to minimize diarrhea and nausea."
        ]
    else:
        protocol = ["General Safety: Monitor for standard adverse drug reactions (ADRs) and hypersensitivity."]

    # 4. Construct the Full Report
    # We format this carefully for the text-to-speech engine
    report = f"""
Clinical Analysis Report for {drug_name}.

1. AI Molecular Intelligence.
Prediction: {percentage:.1f}% Probability of Interaction.
Interpretation: {affinity_desc}.
Target Identified: {target_function}.

2. Physicochemical Suitability.
{physio_summary}

3. Clinical Forecast.
Predicted Phenotype: {phenotype_data.get('phenotype', 'No specific data available.')}

4. Clinical Management Plan.
{'. '.join(protocol)}

5. Recommendation.
Based on the {tone} binding profile and physical properties, this compound { "warrants further investigation" if percentage > 40 else "does not appear to be a viable candidate for this target" }.
    """
    
    return report.strip()