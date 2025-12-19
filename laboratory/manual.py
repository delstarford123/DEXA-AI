import random

# This "Dictionary" acts as your database of manuals.
# Keys match the names in your CSV exactly.
LIBRARY = {
    "Analytical Balance": {
        "manual": "1. Ensure bubble level is centered.\n2. Press TARE to zero.\n3. Add sample gently using a spatula.\n4. Wait for stability indicator (*).",
        "safety": "Max load 200g. Do not overload. Keep clean.",
        "status": "Available",
        "location": "Bench A1"
    },
    "Fume Hood": {
        "manual": "1. Turn on light and fan.\n2. Raise sash to indicated working height.\n3. Keep work 6 inches inside the hood.",
        "safety": "Do not put head inside hood. Keep sash down when not in use.",
        "status": "Running",
        "location": "Wall B"
    },
    "Micro-Pipette": {
        "manual": "1. Set volume dial (do not exceed max).\n2. Press plunger to first stop.\n3. Insert tip into liquid.\n4. Release slowly.",
        "safety": "Never lay pipette flat with liquid in tip. Corrosive hazard.",
        "status": "Available",
        "location": "Drawer C2"
    },
    "Centrifuge": {
        "manual": "1. Balance tubes directly opposite each other.\n2. Secure lid tightly.\n3. Set RPM and Time.\n4. Press Start.",
        "safety": "IMBALANCE HAZARD! Ensure equal weight opposite tubes.",
        "status": "Available",
        "location": "Bench B3"
    },
    "HPLC System": {
        "manual": "1. Degas solvents.\n2. Purge lines.\n3. Check column pressure.\n4. Load sequence in software.",
        "safety": "High pressure! Organic solvents are toxic.",
        "status": "In Use",
        "location": "Analysis Lab 1"
    },
    "Rotary Evaporator": {
        "manual": "1. Attach flask with clip.\n2. Lower into bath.\n3. Turn on vacuum.\n4. Start rotation.",
        "safety": "Implosion risk. Check glassware for cracks.",
        "status": "Available",
        "location": "Synthesis Lab"
    },
    "Magnetic Stirrer": {
        "manual": "1. Place stir bar in beaker.\n2. Center beaker on plate.\n3. Increase speed slowly to create vortex.",
        "safety": "Splash hazard. Wear goggles.",
        "status": "Available",
        "location": "Bench A2"
    },
    "Autoclave": {
        "manual": "1. Check water level.\n2. Load items loosely.\n3. Seal door firmly.\n4. Select 'Sterilize' cycle (121°C).",
        "safety": "HOT STEAM! Wear heat-resistant gloves.",
        "status": "In Use (15m remaining)",
        "location": "Sterilization Room"
    },
    "Spectrophotometer": {
        "manual": "1. Warm up for 15 mins.\n2. Insert blank cuvette and zero.\n3. Insert sample cuvette.\n4. Read Absorbance.",
        "safety": "UV Light hazard. Keep lid closed.",
        "status": "Available",
        "location": "Bench D1"
    },
    "pH Meter": {
        "manual": "1. Remove electrode from storage.\n2. Rinse with distilled water.\n3. Dip in sample.\n4. Wait for smiley face icon.",
        "safety": "Glass electrode is fragile. Handle with care.",
        "status": "Available",
        "location": "Bench A1"
    },
    "Mass Spectrometer": {
        "manual": "1. Check vacuum levels.\n2. Inject sample via inlet.\n3. Monitor ionization peaks.",
        "safety": "High Voltage! Do not open chassis.",
        "status": "Maintenance",
        "location": "Analysis Lab 2"
    },
    "Microscope": {
        "manual": "1. Start with 4x objective.\n2. Place slide on stage.\n3. Adjust coarse focus, then fine focus.",
        "safety": "Cover after use. Do not touch lens.",
        "status": "Available",
        "location": "Microscopy Lab"
    },
    "Incubator": {
        "manual": "1. Set target temp.\n2. Place petri dishes upside down.\n3. Close inner glass door first.",
        "safety": "Biohazard potential. Disinfect spills.",
        "status": "Running",
        "location": "Culture Room"
    },
    "Vortex Mixer": {
        "manual": "1. Set speed dial.\n2. Press tube firmly onto rubber cup.\n3. Pulse for 5-10 seconds.",
        "safety": "Ensure cap is tight on tube.",
        "status": "Available",
        "location": "Bench A2"
    },
    "Freeze Dryer": {
        "manual": "1. Freeze samples solid first.\n2. Turn on condenser.\n3. Turn on vacuum.\n4. Attach flasks.",
        "safety": "Extreme cold. Do not touch condenser coils.",
        "status": "Available",
        "location": "Storage Room"
    },
    "pXRD Instrument": {
        "manual": "1. Prepare powder sample.\n2. Place in holder.\n3. Set 2-theta range.\n4. Start scan.",
        "safety": "X-RAY RADIATION! Door locks automatically.",
        "status": "Available",
        "location": "X-Ray Room"
    },
    "Spectrofluorometer": {
        "manual": "1. Select excitation wavelength.\n2. Insert cuvette.\n3. Close lid completely.\n4. Run scan.",
        "safety": "High intensity light source.",
        "status": "Available",
        "location": "Dark Room"
    },
    "Thermocycler": {
        "manual": "1. Place PCR tubes in wells.\n2. Tighten heated lid.\n3. Select stored protocol.\n4. Run.",
        "safety": "Lid gets very hot (105°C). Do not touch.",
        "status": "In Use",
        "location": "Genomics Lab"
    },
    "Gel Electrophoresis System": {
        "manual": "1. Pour agarose gel.\n2. Load samples into wells.\n3. Connect electrodes (Black to Black, Red to Red).\n4. Run at 100V.",
        "safety": "Electric Shock Hazard! Cover must be on.",
        "status": "Available",
        "location": "Genomics Lab"
    },
    "Water Bath": {
        "manual": "1. Fill with distilled water to fill line.\n2. Set dial to desired temp.\n3. Wait for heating light to toggle.",
        "safety": "Water is HOT. Use tongs for floating racks.",
        "status": "Heating up",
        "location": "Bench B2"
    },
    "Dissolution Tester": {
        "manual": "1. Fill vessels with 900ml media.\n2. Heat to 37°C.\n3. Drop tablets.\n4. Lower paddles/baskets.",
        "safety": "Moving parts. Keep hands clear of shafts.",
        "status": "Available",
        "location": "QC Lab"
    },
    "Disintegration Tester": {
        "manual": "1. Heat water bath to 37°C.\n2. Place tablets in tubes.\n3. Attach disks.\n4. Start basket movement.",
        "safety": "Electrical hazard near water.",
        "status": "Available",
        "location": "QC Lab"
    },
    "Friability Tester": {
        "manual": "1. Weigh tablets (initial).\n2. Load into drum.\n3. Run for 100 revolutions.\n4. Weigh tablets (final).",
        "safety": "Dust hazard. Wear mask.",
        "status": "Available",
        "location": "QC Lab"
    },
    "Tablet Hardness Tester": {
        "manual": "1. Clear debris from jaw.\n2. Place tablet flat.\n3. Press Start.\n4. Record Newtons/Kp.",
        "safety": "Pinch point! Keep fingers away from jaws.",
        "status": "Available",
        "location": "QC Lab"
    },
    "Karl Fischer Titrator": {
        "manual": "1. Check drift value.\n2. Inject sample.\n3. Enter sample weight.\n4. Wait for calculation.",
        "safety": "Reagents contain methanol/iodine. Toxic.",
        "status": "Available",
        "location": "Wet Lab"
    },
    "Melting Point Apparatus": {
        "manual": "1. Fill capillary tube.\n2. Insert into slot.\n3. Ramp temp fast then slow.\n4. Observe melt.",
        "safety": "Surface hot.",
        "status": "Available",
        "location": "Bench A3"
    },
    "Sonicator": {
        "manual": "1. Fill tank with water.\n2. Place beaker in basket.\n3. Set timer.\n4. Lid on.",
        "safety": "High frequency noise. Wear ear protection.",
        "status": "Available",
        "location": "Wash Room"
    },
    "Gas Chromatograph": {
        "manual": "1. Check carrier gas pressure.\n2. Clean syringe.\n3. Inject 1uL sample.\n4. Press Start.",
        "safety": "Hydrogen gas is flammable. Check leaks.",
        "status": "In Use",
        "location": "Analysis Lab 1"
    },
    "Refractometer": {
        "manual": "1. Clean prism with alcohol.\n2. Add drop of liquid.\n3. Look through eyepiece.\n4. Adjust shadow line.",
        "safety": "Do not scratch prism.",
        "status": "Available",
        "location": "Bench A1"
    },
    "Viscometer": {
        "manual": "1. Level the instrument.\n2. Attach correct spindle.\n3. Immerse in fluid.\n4. Set RPM and measure.",
        "safety": "Rotating spindle. Keep hair tied back.",
        "status": "Available",
        "location": "QC Lab"
    },
    "Desiccator": {
        "manual": "1. Slide lid open horizontally.\n2. Place hot crucible inside.\n3. Close lid immediately.",
        "safety": "Vacuum may form. Slide lid, do not pull up.",
        "status": "Available",
        "location": "Bench D2"
    },
    "Hot Plate": {
        "manual": "1. Plug in.\n2. Turn heat dial to setting.\n3. Place beaker on top.",
        "safety": "SURFACE REMAINS HOT AFTER OFF.",
        "status": "Available",
        "location": "Fume Hood"
    },
    "PCR Machine": {
        "manual": "1. Load strip tubes.\n2. Close lid and tighten.\n3. Run 'Standard_35_Cycles'.",
        "safety": "Hot Lid.",
        "status": "Available",
        "location": "Genomics Lab"
    },
    "Biosafety Cabinet": {
        "manual": "1. Turn on blower 15 mins prior.\n2. Wipe surface with 70% ethanol.\n3. Work in center.",
        "safety": "Biohazards inside. Do not block grilles.",
        "status": "Running",
        "location": "Cell Culture"
    },
    "FTIR Spectrometer": {
        "manual": "1. Clean crystal.\n2. Run Background scan.\n3. Place sample on crystal.\n4. Run Sample scan.",
        "safety": "Laser radiation.",
        "status": "Available",
        "location": "Analysis Lab 2"
    },
    "Vacuum Oven": {
        "manual": "1. Place trays inside.\n2. Close door.\n3. Turn on pump.\n4. Open vacuum valve slowly.",
        "safety": "Implosion risk. Do not use cracked glassware.",
        "status": "Available",
        "location": "Synthesis Lab"
    }
}

DEFAULT_MANUAL = {
    "manual": "Standard operating procedures apply. Refer to device sticker.",
    "safety": "General lab safety protocols required.",
    "status": "Unknown",
    "location": "General Storage"
}

def get_equipment_details(name):
    """
    Searches the library for the equipment name.
    If exact match isn't found, it tries to find a partial match.
    """
    # 1. Try exact match
    if name in LIBRARY:
        return LIBRARY[name]
    
    # 2. Try partial match (e.g., "Microscope" matches "Digital Microscope")
    for key in LIBRARY:
        if name.lower() in key.lower() or key.lower() in name.lower():
            return LIBRARY[key]
            
    # 3. Return default if unknown
    return DEFAULT_MANUAL