1. Pharmacological Database & Nomenclature
Dual Nomenclature Indexing: The system indexes pharmaceutical agents by both their Proprietary (Market) Names and International Nonproprietary Names (INN/Chemical Names) to ensure accurate identification across different medical standards.

Structural Retrieval: Users can query the database to retrieve precise chemical structures, facilitating the visualization of molecular geometry and functional groups.

2. Drug Interaction & Polypharmacy Analysis
Interaction Prediction Algorithms: When multiple agents are input (e.g., Polypharmacy scenarios), the system analyzes potential Drug-Drug Interactions (DDIs).

Pharmacodynamic Modeling: It predicts the physiological sequelae of these interactions, elucidating mechanisms such as synergistic toxicity, receptor antagonism, or metabolic inhibition (e.g., Cytochrome P450 pathways).

3. De Novo Drug Discovery Module
Structural Optimization: The "Discovery Studio" allows for in silico modification of pharmacophores. Users can append specific functional groups (e.g., Hydroxyl, Carboxyl) to an existing scaffold to synthesize novel analogues.

Combinatorial Synthesis Simulation: Algorithms simulate the conjugation of distinct molecular entities, generating theoretical new drug candidates for viability assessment.

4. Virtual Research Laboratory
Digital Wet Lab: A module dedicated to simulating the physical laboratory environment, cataloging essential Instrumentation (spectrometers, centrifuges) and Reagents required for synthesis.

Research Repository: Integration of a database for retrieving peer-reviewed Drug Research Papers, ensuring users have access to evidence-based clinical data.
# DexaGen-AI: Advanced Pharmacological Intelligence System

![DexaGen Banner](https://via.placeholder.com/1000x300?text=DexaGen-AI+Clinical+Hub)

## 🔬 Overview

**DexaGen-AI** is a state-of-the-art computational platform designed for **in silico drug discovery**, **clinical interaction analysis**, and **pharmacological research**. By integrating molecular dynamics engines with generative AI, DexaGen-AI provides researchers and clinicians with a powerful tool to visualize chemical structures, predict drug-drug interactions (DDIs), and simulate novel compound synthesis.

## 🌟 Key Features

* **Discovery Studio:**
    * **3D Molecular Visualization:** Real-time rendering of small molecules using the `3Dmol.js` engine.
    * **Novelty Analysis:** Evaluate new chemical entities using **QED** (Quantitative Estimate of Drug-likeness) and **Lipinski's Rule of 5** metrics.
    * **"Magic Wand" Editor:** Interactive addition of functional groups for rapid structural prototyping.

* **Clinical Dashboard:**
    * **Polypharmacy Interaction Checker:** Analyze multi-drug regimens for potential contraindications, toxicity risks, and synergistic effects.
    * **AI Clinical Reports:** Generate detailed, medically accurate reports on binding affinity, dosage guidelines, and safety profiles.

* **Virtual Laboratory:**
    * Simulation of physical laboratory environments, including equipment and reagent management.
    * Access to a curated database of drug research papers and clinical trial data.

## 🛠️ Technology Stack

* **Backend:** Python (Flask), RDKit (Cheminformatics), Google Gemini AI (LLM).
* **Frontend:** HTML5, CSS3 (Bootstrap 5, Glassmorphism UI), JavaScript (jQuery, Chart.js).
* **Visualization:** 3Dmol.js (Molecular rendering).
* **Data Source:** PubChemPy, Custom Curated Datasets (`raw_drug_data.csv`).

## ⚙️ Installation & Setup ("Lightworks")

Follow these steps to deploy the DexaGen-AI system locally.

### Prerequisites

* **Python 3.8+**
* **pip** (Python Package Installer)
* **Git**

### 1. Clone the Repository

```bash
git clone [https://github.com/your-username/dexagen-ai.git](https://github.com/your-username/dexagen-ai.git)
cd dexagen-ai

2. Set Up Virtual Environment (Recommended)
Bash

# Windows
python -m venv venv
.\venv\Scripts\activate

# macOS/Linux
python3 -m venv venv
source venv/bin/activate
3. Install Dependencies
Install all required Python libraries.

Bash

pip install -r requirements.txt
4. Configuration
Create a .env file in the root directory to store your API keys safely.

Code snippet

# .env file content
GEMINI_API_KEY=your_google_gemini_api_key_here
FLASK_ENV=development
SECRET_KEY=your_secure_random_string
5. Run the Application
Bash

python main.py
The server will start at http://127.0.0.1:5000/.

Open your web browser and navigate to this URL to access the Clinical Hub.

📦 Dependencies (requirements.txt)
Ensure your requirements.txt file includes the following libraries:

Plaintext

Flask==3.0.0
rdkit-pypi==2023.9.5
pubchempy==1.0.4
google-generativeai==0.3.2
python-dotenv==1.0.0
numpy==1.26.0
pandas==2.1.1
requests==2.31.0
gunicorn==21.2.0
📂 Project Structure
DEXAGEN-AI/
├── models/                # Pre-trained ML models (.pkl)
├── scripts/               # Core logic engines
│   ├── further_details.py # Clinical report generator
│   ├── pathway_engine.py  # Biological target mapping
│   ├── simulation_engine.py # MD simulation logic
│   └── video_engine.py    # AI Video generation
├── static/                # Static assets (CSS, JS, Images, Videos)
│   └── video.mp4          # Background video asset
├── templates/             # HTML Frontend
│   ├── dashboard.html     # Main Clinical Hub
│   └── discovery.html     # Drug Discovery Studio
├── .env                   # Environment variables (GitIgnored)
├── main.py                # Flask Application Entry Point
├── raw_drug_data.csv      # Drug Database
├── requirements.txt       # Python Dependencies
└── README.md              # Project Documentation
🤝 Contribution
Contributions are welcome! Please follow these steps:

Fork the repository.

Create a new feature branch (git checkout -b feature/AmazingFeature).

Commit your changes (git commit -m 'Add some AmazingFeature').

Push to the branch (git push origin feature/AmazingFeature).

Open a Pull Request.

📜 License
Distributed under the MIT License. See LICENSE for more information.

Disclaimer: DexaGen-AI is a research tool. Predictions regarding drug interactions and clinical safety should always be verified by qualified medical professionals