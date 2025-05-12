import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from rdkit import Chem
from rdkit.Chem import Descriptors
# import pkpdsim # Placeholder, will use specific functions later

# --- Configuration & Constants ---
VERSION = "0.1.0"
APP_TITLE = "Outil de Simulation Pharmacocinétique Clinique"
FRENCH_LABELS = {
    "smiles_input": "SMILES (Structure Moléculaire):",
    "crcl_input": "Clairance de la Créatinine (CrCl, mL/min):",
    "hepatic_function_input": "Fonction Hépatique (Child-Pugh Score):",
    "weight_input": "Poids (kg):",
    "age_input": "Âge (années):",
    "dose_input": "Dose (mg):",
    "bioavailability_input": "Biodisponibilité (F, %):",
    "administration_route_input": "Voie d'administration:",
    "iv_bolus": "IV Bolus",
    "oral": "Orale",
    "run_simulation_button": "Lancer la Simulation",
    "simulation_results_header": "Résultats de la Simulation PK",
    "pk_parameters_header": "Paramètres Pharmacocinétiques Clés",
    "cmax": "Cmax (Concentration Maximale)",
    "tmax": "Tmax (Temps pour Cmax)",
    "t_half": "T1/2 (Demi-vie d'Élimination)",
    "auc": "AUC (Aire Sous la Courbe)",
    "vd": "Vd (Volume de Distribution)",
    "cl": "CL (Clairance)",
    "clinical_summary_header": "Résumé Clinique Standardisé",
    "monitoring_parameters": "Paramètres de Surveillance Clés pour les Prescripteurs:",
    "drug_interactions": "Interactions Médicamenteuses Potentielles:",
    "dose_adjustments": "Ajustements Posologiques Recommandés:",
    "therapeutic_window": "Fenêtre Thérapeutique:",
    "disclaimer_header": "Avertissements et Limites de la Simulation",
    "disclaimer_text": (
        "Cet outil est destiné à des fins éducatives et de recherche uniquement. "
        "Les simulations sont basées sur des modèles simplifiés et des hypothèses générales. "
        "NE PAS utiliser pour prendre des décisions cliniques réelles sans la supervision "
        "d'un professionnel de santé qualifié. Les résultats peuvent ne pas refléter "
        "la variabilité interindividuelle complète ni toutes les interactions médicamenteuses possibles."
    ),
    "plasma_concentration_time_curve": "Courbe Concentration Plasmatique - Temps",
    "time_hours": "Temps (heures)",
    "concentration_units": "Concentration (µg/mL)", # Assuming units, can be made dynamic
    "molecular_weight": "Poids Moléculaire (g/mol)",
    "logp": "LogP (Lipophilie)",
}

# --- Pharmacokinetic Model ---
def one_compartment_model(t, C, Ka, Ke, Vd, Dose, F):
    """Modèle pharmacocinétique à un compartiment."""
    # Pour l'instant, on suppose une administration orale avec absorption de premier ordre
    # ou un bolus IV. La différenciation se fera dans la fonction de simulation.
    if Ka is None: # IV Bolus
        dCdt = -Ke * C
    else: # Oral
        # This is a simplified representation for dC/dt after oral admin.
        # A more accurate model would involve a gut compartment.
        # For now, let's assume input into central compartment is F*Dose*Ka*exp(-Ka*t)
        # and elimination is Ke*C.
        # However, solve_ivp expects dC/dt.
        # Let A_gut be amount in gut, A_central be amount in central.
        # dA_gut/dt = -Ka * A_gut
        # dA_central/dt = Ka * A_gut - Ke * A_central
        # C = A_central / Vd
        # This requires a system of ODEs.
        # For a simpler approach for now, let's assume the input is directly affecting concentration
        # This is not entirely correct but will be refined.
        # A common simplification for plotting is to calculate absorption phase separately
        # or use a combined function.
        # For now, let's use a placeholder for oral, focusing on IV first.
        # This part needs significant refinement for oral model.
        absorption_rate = F * Dose * Ka * np.exp(-Ka * t) / Vd # Incorrect, placeholder
        elimination_rate = Ke * C
        dCdt = absorption_rate - elimination_rate # Placeholder for oral
    return dCdt

def simulate_pk_profile(smiles, crcl, hepatic_func, weight, age, dose, bioavailability, admin_route):
    """Simule le profil PK et calcule les paramètres."""
    # Placeholder for actual simulation logic
    # 1. Estimate Vd, CL from SMILES and patient data (complex part, needs literature/models)
    # For now, use dummy values or simple correlations
    
    mol = Chem.MolFromSmiles(smiles)
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)

    # Simplified estimations (THESE ARE VERY CRUDE - FOR DEMO ONLY)
    Vd_base = 0.7 * weight # L/kg * kg = L
    CL_base = 0.1 * weight # L/hr/kg * kg = L/hr (very rough estimate)

    # Adjust CL based on CrCl (simplified Cockcroft-Gault like adjustment)
    # Normal CrCl ~100-120 mL/min.
    # CL_renal_fraction = 0.7 # Assume 70% renal clearance for this hypothetical drug
    # CL_adjusted = CL_base * ( (crcl / 100) * CL_renal_fraction + (1 - CL_renal_fraction) )
    # CL_adjusted = max(CL_adjusted, 0.01 * weight) # Ensure CL is not zero

    # More direct: Ke is influenced by CL and Vd (Ke = CL/Vd)
    # Let's try to estimate Ke directly or CL first.
    # For a 1-compartment model: CL = Ke * Vd
    
    # Example: If we assume a T1/2 of 4 hours for a standard patient
    # Ke_standard = np.log(2) / 4  # per hour
    # CL_standard = Ke_standard * Vd_base
    
    # Adjust Ke based on CrCl (very simplified)
    # This is highly drug-specific. For now, a linear scaling.
    Ke_estimated = (np.log(2) / 8) * (crcl / 100 if crcl > 0 else 0.1) # Highly simplified, assumes 8hr T1/2 at normal CrCl
    Ke_estimated = max(Ke_estimated, 0.01) # Avoid zero Ke

    Vd_estimated = Vd_base # L
    CL_estimated = Ke_estimated * Vd_estimated # L/hr

    F = bioavailability / 100.0

    t_eval = np.linspace(0, 48, 500) # Simuler sur 48 heures
    C0 = 0 # Initial concentration in plasma

    params_iv = (None, Ke_estimated, Vd_estimated, dose, 1.0) # Ka=None for IV, F=1 for IV
    params_oral = (1.0, Ke_estimated, Vd_estimated, dose, F) # Assuming Ka=1.0 hr^-1 for oral (NEEDS REFINEMENT)

    if admin_route == FRENCH_LABELS["iv_bolus"]:
        # For IV bolus, C0 is Dose / Vd, and then it declines.
        # The ODE dC/dt = -Ke*C starts from C(0) = Dose/Vd.
        C0_iv = dose / Vd_estimated
        sol = solve_ivp(one_compartment_model, [t_eval[0], t_eval[-1]], [C0_iv],
                        args=params_iv, dense_output=True, t_eval=t_eval, method='LSODA')
        C_t = sol.y[0]
        # Adjust Tmax for IV Bolus
        tmax_val = 0
        cmax_val = C0_iv
    else: # Oral
        # This requires a proper 2-compartment model (gut and central) or a 1-compartment with absorption phase
        # Using a simplified approach for now, will need pkpdsim or more complex ODE for oral
        st.warning("Le modèle oral est hautement simplifié et en cours de développement.")
        # For oral, C(0) = 0. The model needs to handle absorption.
        # Let's use a very basic analytical solution for 1-comp oral if possible, or stick to numerical
        # C(t) = (F * Dose * Ka) / (Vd * (Ka - Ke)) * (np.exp(-Ke * t) - np.exp(-Ka * t))
        Ka_oral = params_oral[0] # 1.0 hr^-1
        if Ka_oral == Ke_estimated: # Avoid division by zero, use a slightly different Ka
            Ka_oral += 0.01

        C_t = (F * dose * Ka_oral) / (Vd_estimated * (Ka_oral - Ke_estimated)) * \
              (np.exp(-Ke_estimated * t_eval) - np.exp(-Ka_oral * t_eval))
        C_t = np.maximum(C_t, 0) # Ensure no negative concentrations due to numerical issues

        if len(C_t) > 0:
            tmax_idx = np.argmax(C_t)
            tmax_val = t_eval[tmax_idx]
            cmax_val = C_t[tmax_idx]
        else:
            tmax_val = 0
            cmax_val = 0


    # Calculate PK parameters
    if len(C_t) > 0 and cmax_val > 0:
        # T1/2
        t_half_val = np.log(2) / Ke_estimated if Ke_estimated > 0 else np.inf
        # AUC (trapezoidal rule)
        auc_val = np.trapz(C_t, t_eval)
    else:
        t_half_val = np.inf
        auc_val = 0
        cmax_val = 0
        tmax_val = 0


    pk_params = {
        FRENCH_LABELS["cmax"]: f"{cmax_val:.2f} µg/mL",
        FRENCH_LABELS["tmax"]: f"{tmax_val:.2f} h",
        FRENCH_LABELS["t_half"]: f"{t_half_val:.2f} h",
        FRENCH_LABELS["auc"]: f"{auc_val:.2f} µg.h/mL",
        FRENCH_LABELS["vd"]: f"{Vd_estimated:.2f} L",
        FRENCH_LABELS["cl"]: f"{CL_estimated:.2f} L/h ({CL_estimated*1000/60:.2f} mL/min)",
        FRENCH_LABELS["molecular_weight"]: f"{mw:.2f} g/mol",
        FRENCH_LABELS["logp"]: f"{logp:.2f}",
    }

    return t_eval, C_t, pk_params

# --- Streamlit UI ---
st.set_page_config(page_title=APP_TITLE, layout="wide")
st.title(f"{APP_TITLE} (v{VERSION})")

st.sidebar.header("Paramètres du Patient et du Médicament")
smiles_input = st.sidebar.text_input(FRENCH_LABELS["smiles_input"], "CCO") # Ethanol example
dose_mg = st.sidebar.number_input(FRENCH_LABELS["dose_input"], min_value=0.1, value=100.0, step=10.0)
admin_route_input = st.sidebar.selectbox(FRENCH_LABELS["administration_route_input"],
                                         [FRENCH_LABELS["iv_bolus"], FRENCH_LABELS["oral"]])
bioavailability_pct = st.sidebar.slider(FRENCH_LABELS["bioavailability_input"], 0, 100, 70,
                                        disabled=(admin_route_input == FRENCH_LABELS["iv_bolus"]))


st.sidebar.subheader("Données Cliniques du Patient")
weight_kg = st.sidebar.number_input(FRENCH_LABELS["weight_input"], min_value=1.0, value=70.0, step=1.0)
age_years = st.sidebar.number_input(FRENCH_LABELS["age_input"], min_value=1, value=50, step=1)
crcl_ml_min = st.sidebar.number_input(FRENCH_LABELS["crcl_input"], min_value=0.0, value=100.0, step=5.0)
# hepatic_func_score = st.sidebar.selectbox(FRENCH_LABELS["hepatic_function_input"], ["A (5-6 points)", "B (7-9 points)", "C (10-15 points)"])
# For now, hepatic function is not directly used in the simplified model, but good to have for future.
st.sidebar.info("La fonction hépatique n'est pas encore intégrée dans ce modèle simplifié.")


if st.sidebar.button(FRENCH_LABELS["run_simulation_button"]):
    error_message = ""
    if not smiles_input:
        error_message = "Veuillez entrer une structure SMILES valide."
    else:
        try:
            mol = Chem.MolFromSmiles(smiles_input)
            if mol is None:
                error_message = "SMILES invalide. Impossible de générer la molécule."
        except Exception as e:
            error_message = f"Erreur de validation SMILES: {e}"

    if error_message:
        st.error(error_message)
    else:
        st.header(FRENCH_LABELS["simulation_results_header"])
        
        # Effective bioavailability
        effective_bioavailability = 100.0 if admin_route_input == FRENCH_LABELS["iv_bolus"] else bioavailability_pct

        t_sim, c_sim, pk_params_results = simulate_pk_profile(
            smiles=smiles_input,
            crcl=crcl_ml_min,
            hepatic_func=None, # hepatic_func_score,
            weight=weight_kg,
            age=age_years,
            dose=dose_mg,
            bioavailability=effective_bioavailability,
            admin_route=admin_route_input
        )

        # Plotting
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(t_sim, c_sim, label="Concentration Plasmatique", color="blue", linewidth=2)
        ax.set_title(FRENCH_LABELS["plasma_concentration_time_curve"], fontsize=14)
        ax.set_xlabel(FRENCH_LABELS["time_hours"], fontsize=12)
        ax.set_ylabel(FRENCH_LABELS["concentration_units"], fontsize=12)
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.legend()

        # Annotations (basic for now)
        cmax_val_num = float(pk_params_results[FRENCH_LABELS["cmax"]].split(" ")[0])
        tmax_val_num = float(pk_params_results[FRENCH_LABELS["tmax"]].split(" ")[0])

        if cmax_val_num > 0:
            ax.plot(tmax_val_num, cmax_val_num, 'ro') # Mark Cmax
            ax.annotate(f"Cmax: {cmax_val_num:.2f}", (tmax_val_num, cmax_val_num),
                        textcoords="offset points", xytext=(0,10), ha='center', color="red")

        st.pyplot(fig)

        st.subheader(FRENCH_LABELS["pk_parameters_header"])
        cols = st.columns(2)
        i = 0
        for key, value in pk_params_results.items():
            cols[i % 2].metric(label=key, value=str(value))
            i += 1
        
        st.subheader(FRENCH_LABELS["clinical_summary_header"])
        st.markdown(f"**{FRENCH_LABELS['monitoring_parameters']}**")
        st.info("Surveiller les signes d'efficacité et de toxicité. Ajuster la dose en fonction de la réponse clinique et des paramètres PK.")

        st.markdown(f"**{FRENCH_LABELS['drug_interactions']}**")
        st.warning("Non implémenté. Consulter une base de données d'interactions médicamenteuses (ex: Lexicomp, Micromedex).")

        st.markdown(f"**{FRENCH_LABELS['dose_adjustments']}**")
        st.info(f"Une clairance rénale de {crcl_ml_min} mL/min a été utilisée. "
                f"Des ajustements peuvent être nécessaires pour les patients avec une fonction rénale altérée. "
                f"Le Vd estimé est de {pk_params_results[FRENCH_LABELS['vd']]}, CL de {pk_params_results[FRENCH_LABELS['cl']]}.")

        st.markdown(f"**{FRENCH_LABELS['therapeutic_window']}**")
        st.warning("Non implémenté. Dépend de chaque médicament spécifique.")

st.markdown("---")
st.subheader(FRENCH_LABELS["disclaimer_header"])
st.warning(FRENCH_LABELS["disclaimer_text"])
st.markdown(f"<p style='text-align: center;'>{APP_TITLE} v{VERSION}</p>", unsafe_allow_html=True)

# To run: streamlit run app.py
