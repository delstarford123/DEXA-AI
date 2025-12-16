DEXA AI: Advanced Pharmacological Intelligence
"From Statistical Guesswork to Mechanistic Certainty."

DexaGen AI is a multi-scale computational platform designed to revolutionize the drug development pipeline. By integrating Structural Bioinformatics, Quantum-Level Molecular Modeling, and Multimodal Deep Learning, DexaGen bridges the "Safety Gap" in pre-clinical screening, predicting adverse Drug-Drug Interactions (DDIs) and binding affinities without the latency of in-vivo experimentation.
Table of Contents
Executive Summary

The Scientific Paradigm

System Architecture

Mathematical & Algorithmic Core

Installation & Usage

Roadmap & Phases

Investigators & Funding

Executive SummaryThe pharmaceutical industry faces a critical bottleneck: the inability to rapidly predict adverse interactions in polypharmacy. Traditional wet-lab assays are capital-intensive and slow.DexaGen AI offers a purely digital, structure-based screening protocol. By analyzing the SMILES (Simplified Molecular Input Line Entry System) string of a molecule, the system provides an instantaneous Risk Probability Score, identifying:Reactive functional groups via Quantum Chemical Modeling.Steric hindrance via Real-time 3D Visualization.Synergistic or Antagonistic outcomes via Random Forest & Graph Neural Networks.⚛️ The Scientific ParadigmDexaGen operates on a Multi-Scale Hybrid Modeling Pipeline, moving beyond simple "Black Box" predictions to explainable, mechanism-forward insights.LayerFunctionTechnology1. Structural AIUnderstands drug chemistry from raw structure.Quantum Chemistry (DFT), SMILES Parsing2. Molecular InteractionPredicts docking, binding affinity, and DDIs.Molecular Dynamics (MD), 3D-CNNs, RDKit3. Pathological OutcomeForecasts disease-modifying effects & pathway switching.Systems Biology, Knowledge Graphs, LLM IntegrationKey CapabilitiesEx-Vivo Chemical Interaction Profiling: Detection of chemical incompatibilities outside biological systems.Systemic Pharmacokinetic Modeling: Simulation of competitive inhibition (e.g., Cortisol displacement at Glucocorticoid Receptors).Metabolic Interference: Prediction of Cytochrome P450 enzyme inhibition.🏗 System ArchitectureThe DexaGen backend is built on a lightweight Flask (Python) micro-framework, ensuring stateless transactions ideal for cloud deployment (AWS/GCP). It integrates Google's Gemini API to translate raw probability scores into natural-language clinical reports.
graph TD
    A[User Input: SMILES String] --> B(Flask API Endpoint);
    B --> C{Feature Extraction};
    C -->|RDKit| D[Morgan Fingerprints ECFP4];
    D --> E[Random Forest / GNN Model];
    E --> F[Interaction Probability Score];
    F --> G[Gemini LLM];
    G --> H[Clinical Narrative Report];
    E --> I[3Dmol.js Renderer];
    H & I --> J[Frontend Dashboard];

    Mathematical & Algorithmic CoreDexaGen transforms abstract chemical notations into high-dimensional vector spaces.1. Feature Extraction (Morgan Fingerprints)We utilize Extended-Connectivity Fingerprints (ECFP4). The algorithm iterates through each atom to generate a unique identifier hash based on its circular neighborhood.The Tanimoto Similarity Coefficient ($T_c$) is used to determine structural similarity between a novel compound ($A$) and known interactors ($B$):$$T_c(A, B) = \frac{|A \cap B|}{|A \cup B|} = \frac{c}{a + b - c}$$Where:$c$ is the number of common bits set in both fingerprints.$a$ is the number of bits set in molecule $A$.$b$ is the number of bits set in molecule $B$.2. Implementation SnippetFrom Technical Appendix A: Feature Extraction Logic

    from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

# Initialize Generator: Radius 2 approximates bond connectivity
MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def smiles_to_fingerprint(smiles):    
    '''Converts SMILES string to 2048-bit Numpy Vector'''    
    mol = Chem.MolFromSmiles(smiles)    
    if mol is not None:
        return MORGAN_GEN.GetFingerprintAsNumPy(mol)    
    return None

    3. Inference Engine
We employ an ensemble learning method (Random Forest) for robust classification of interaction types, prioritizing interpretability over black-box deep learning for the initial screening phase.

💻 Installation & Usage
Prerequisites
Python 3.9+

RDKit

Scikit-Learn

Gemini API Key

Setup

# Clone the repository
git clone https://github.com/MMUST-Medicine/DexaGen-AI.git

# Install dependencies
pip install -r requirements.txt

# Set Environment Variables
export FLASK_APP=app.py
export GEMINI_API_KEY='your_api_key_here'

# Run the Flask Server
flask run
API Usage Example
Send a POST request with a SMILES string (e.g., for Dexamethasone):

curl -X POST http://localhost:5000/predict \
     -H "Content-Type: application/json" \
     -d '{"smiles": "C[C@]12C[C@@H]([C@H]3[C@H]([C@@H]1CC[C@]2(C(=O)CO)O)CCC4=CC(=O)C=C[C@]34C)F"}'

     Roadmap & Phases
The project is structured over a 36-month timeline.

Phase 1 (Months 1-6): Conceptualization, Data Acquisition (ChEMBL, PubChem), and Governance Framework.

Phase 2 (Months 7-15): Development of Molecular Simulation Engine and Quantum Chemical pipelines.

Phase 3 (Months 16-24): Dexa AI Core Integration (API construction and visual dashboards).

Phase 4 (Months 25-30): Large-Scale Validation against FDA drug labels and pharmacovigilance reports.

Phase 5 (Months 31-36): Public Beta Release and Dissemination.

👥 Investigators & Funding
Principal Investigators:

Prof. Barack Otieno Abonyo (MMUST School of Medicine)

Delstarford Isaiah (Project Lead & AI Architect)

Contact:

📍 MMUST School of Medicine, P.O. BOX 190-50100, Kakamega, Kenya

📧 babonyo45@gmail.com | delstarfordisaiah@gmail.com

Funding Agency:

MMUST PPRI Division (Budget: USD 1,987,450)

References & Inspirations:

AtomNet (Wallach et al., 2015)

AlphaFold (Jumper et al., 2021)

DeepDDI (Ryu et al., 2018)

© 2025 Dexa AI. All Rights Reserved. Transforming In-Silico Potential into Real-World Clinical Salvation.