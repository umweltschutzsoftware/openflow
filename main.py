import streamlit as st
from math import pi, log10, sqrt

# -----------------------------
# Physikalische Konstanten / Annahmen
# -----------------------------
GAMMA = 1.30          # Isentropenexponent Erdgas ~ Methan
R = 520.0             # spezifische Gaskonstante J/(kg·K)
T0 = 288.0            # Umgebungstemperatur (K) ~ 15 °C
CD = 0.8              # Ausflussbeiwert
RHO_N = 0.8           # Dichte Erdgas bei Normbedingungen (kg/m³)
ETA_ACOUSTIC = 1e-3   # Anteil Strömungsleistung -> Schallleistung
P_REF = 1e-12         # Referenz-Schallleistung (W)


# -----------------------------
# Hilfsfunktionen
# -----------------------------
def diameter_to_area(diameter_value: float, unit: str) -> float:
    """
    Wandelt einen Durchmesser in m² Querschnitt um.
    unit: "inch" oder "mm"
    """
    if unit == "inch":
        d_m = diameter_value * 25.4e-3  # inch -> m
    elif unit == "mm":
        d_m = diameter_value / 1000.0   # mm -> m
    else:
        raise ValueError("Unbekannte Einheit für Durchmesser")

    area = pi * (d_m ** 2) / 4.0
    return area


def pressure_to_absolute(p_value: float, kind: str) -> float:
    """
    Wandelt Druck in bar(abs) um.
    kind: "gauge" = Überdruck, "absolute" = Absolutdruck
    """
    if kind == "gauge":
        return p_value + 1.0  # 1 bar Atmosphärendruck addieren
    elif kind == "absolute":
        return p_value
    else:
        raise ValueError("Unbekannte Druckart")


def compute_flows_and_acoustics(p_bar_abs: float, area_m2: float, distances_m):
    """
    Berechnet Massenstrom, Volumenströme und akustische Größen
    für kritische Ausströmung eines idealisierten Erdgasstroms.
    p_bar_abs: Absolutdruck in bar
    area_m2: Querschnitt in m²
    distances_m: Liste von Abständen in m, für die L_p berechnet werden soll
    """

    # Druck in Pa
    p0 = p_bar_abs * 1e5

    # Kritische Ausströmung: Massenstrom (ideales Gas, choked flow)
    # Formel: m_dot = C_d * A * p0 * sqrt( gamma / (R*T0) ) * (2/(gamma+1))^((gamma+1)/(2*(gamma-1)))
    critical_factor = (2.0 / (GAMMA + 1.0)) ** ((GAMMA + 1.0) / (2.0 * (GAMMA - 1.0)))
    m_dot = CD * area_m2 * p0 * sqrt(GAMMA / (R * T0)) * critical_factor  # kg/s

    # Kritischer Zustand (T*, p*, rho*)
    T_star = T0 * 2.0 / (GAMMA + 1.0)
    p_star = p0 * (2.0 / (GAMMA + 1.0)) ** (GAMMA / (GAMMA - 1.0))
    rho_star = p_star / (R * T_star)

    # Volumenströme
    V_dot_crit_m3_s = m_dot / rho_star
    V_dot_crit_m3_h = V_dot_crit_m3_s * 3600.0

    V_dot_N_m3_s = m_dot / RHO_N
    V_dot_N_m3_h = V_dot_N_m3_s * 3600.0

    # Jet-Leistung
    a_star = sqrt(GAMMA * R * T_star)  # Schallgeschwindigkeit im kritischen Zustand
    P_jet = 0.5 * m_dot * (a_star ** 2)  # W

    # Akustische Leistung
    P_ac = ETA_ACOUSTIC * P_jet
    L_W = 10.0 * log10(P_ac / P_REF)  # dB

    # Schalldruckpegel an gewünschten Entfernungen
    L_p = {}
    for r in distances_m:
        if r <= 0:
            continue
        # Lp(r) = Lw - 10 log10(2*pi*r^2) Halbkugel-Ausbreitung
        L_p[r] = L_W - 10.0 * log10(2.0 * pi * (r ** 2))

    return {
        "m_dot": m_dot,
        "V_dot_crit_m3_h": V_dot_crit_m3_h,
        "V_dot_N_m3_h": V_dot_N_m3_h,
        "P_jet": P_jet,
        "P_ac": P_ac,
        "L_W": L_W,
        "L_p": L_p,
        "T_star": T_star,
        "p_star_bar": p_star / 1e5,
    }


def fmt(x, digits=0):
    """Einfache Formatierung mit fester Nachkommastellenzahl."""
    return f"{x:.{digits}f}".replace(".", ",")


def generate_report(
    p_input: float,
    p_kind: str,
    diameter: float,
    d_unit: str,
    distances_m,
    results: dict,
) -> str:
    """
    Erzeugt einen gutachtentauglichen Berichtstext auf Basis
    der Eingabedaten und berechneten Ergebnisse.
    """

    p_abs = pressure_to_absolute(p_input, p_kind)

    # Kleine sprechende Strings
    druck_bez = "Überdruck" if p_kind == "gauge" else "Absolutdruck"
    d_unit_str = 'Zoll' if d_unit == 'inch' else 'mm'

    m_dot = results["m_dot"]
    V_crit = results["V_dot_crit_m3_h"]
    V_N = results["V_dot_N_m3_h"]
    L_W = results["L_W"]
    L_p = results["L_p"]

    # Text-Baustein für die Abstände
    if L_p:
        pegelliste = "\n".join(
            f"- im Abstand {fmt(r, 1)} m: ca. {fmt(L_p[r], 0)} dB"
            for r in sorted(L_p.keys())
        )
    else:
        pegelliste = "- (keine gültigen Abstände angegeben)"

    text = f"""
Überschlägige Abschätzung von Gasstrom und Geräuschemissionen bei offenem Gasaustritt

1. Beschreibung des Szenarios
Im vorliegenden Fall wird ein plötzlicher Austritt von Erdgas aus einer offenen Rohrleitung betrachtet. 
Ziel ist eine überschlägige Ermittlung der ausströmenden Gasmenge sowie der zu erwartenden Geräuschemissionen. 
Die Berechnung dient einer ersten technischen Bewertung und erfolgt konservativ.

2. Eingangsdaten

- Medium: Erdgas (vereinfachend als Methan behandelt)
- Rohrdurchmesser: {fmt(diameter, 1)} {d_unit_str}
- Betriebsdruck in der Leitung: {fmt(p_input, 1)} bar ({druck_bez})
- Umgebung: Atmosphärischer Druck (ca. 1 bar), Umgebungstemperatur ca. 15 °C
- Zur Berechnung wird der Betriebsdruck in einen Absolutdruck von {fmt(p_abs, 1)} bar(abs) überführt.

3. Grundannahmen der Berechnung

Für eine praxisnahe und gleichzeitig konservative Abschätzung werden folgende Annahmen zugrunde gelegt:

- Aufgrund des vorliegenden Druckverhältnisses wird eine kritische Ausströmung angenommen. 
  Das Erdgas verlässt die Öffnung dabei mit sehr hoher Strömungsgeschwindigkeit.
- Verluste durch Rauigkeit oder Detaileffekte an der Rohröffnung werden pauschal über einen Ausflussbeiwert berücksichtigt.
- Für die Geräuschentstehung wird angenommen, dass ein kleiner Teil der Strömungsenergie in Schall umgewandelt wird. 
  Dieser Ansatz bildet reale Gasstrahlen erfahrungsgemäß in guter Näherung ab.
- Die Schallausbreitung wird als Freifeldsituation betrachtet. Abschirmungen oder Reflexionen werden nicht berücksichtigt, 
  sodass die berechneten Werte in der Regel als konservativ angesehen werden können.

4. Ergebnisse der überschlägigen Berechnung

Auf Basis der kinetischen Strömungsleistung des Gasstrahls und eines konservativen Ansatzes zur Umsetzung 
von Strömungsenergie in Schall ergibt sich überschlägig:

- Schallleistungspegel der Gasstrahlquelle: ca. {fmt(L_W, 0)} dB

Unter Annahme einer halbkugelförmigen Ausbreitung im Freifeld ergeben sich folgende Schalldruckpegel:

{pegelliste}

Die A-bewerteten Pegel liegen in der Regel um wenige Dezibel darunter, 
da Gasstrahlen einen nennenswerten Anteil ihrer Schallenergie im Mittel- und Hochfrequenzbereich aufweisen.

5. Einordnung

Die berechneten Werte sind als Größenordnungsabschätzung zu verstehen. 
Aufgrund der konservativen Annahmen ist davon auszugehen, dass die tatsächlichen Werte eher etwas niedriger liegen. 
Für detaillierte Immissionsprognosen, die z. B. in der Genehmigungsplanung anzuwenden sind, 
sind ergänzend weitere Dämpfungsmechanismen (z. B. Abschirmung, Gebäude, Bodeneffekte) zu berücksichtigen 
oder alternativ Messungen bzw. weiterführende Berechnungen (z. B. frequenzaufgelöste Betrachtungen) durchzuführen.
"""
    return text.lstrip()


# -----------------------------
# Streamlit UI
# -----------------------------
def main():
    st.set_page_config(page_title="Gasaustritt – Überschlägige Berechnung", layout="wide")

    st.title("Überschlägige Berechnung der Geräuschemissionen und Geräuschimmissionen in einem Open Flow Szenario bei Gasaustritt")
    st.markdown(
        "Diese Anwendung erstellt aus Druck- und Durchmesserangaben "
        "einen standardisierten Textbaustein zur überschlägigen Bewertung der" \
        "Geräuschimmissionen in angegebenen Abständen."
    )

    with st.form("input_form"):
        col1, col2 = st.columns(2)

        with col1:
            p_value = st.number_input("Druck in der Leitung (bar)", value=17.0, min_value=0.1, step=0.5)
            #p_kind = st.radio(
            #    "Art des Drucks",
            #    options=["gauge", "absolute"],
            #    format_func=lambda x: "Überdruck [bar(g)]" if x == "gauge" else "Absolutdruck [bar(abs)]",
            #)
            p_kind="gauge"

        with col2:
            diameter = st.number_input("Rohrdurchmesser (Zoll)", value=2.0, min_value=0.1, step=0.1)
            #d_unit = st.radio(
            #    "Einheit des Durchmessers",
            #    options=["inch", "mm"],
            #    format_func=lambda x: "Zoll" if x == "inch" else "mm",
            #)
            d_unit="inch"

        distances_str = st.text_input(
            "Abstände für Geräuschimmission [m] (durch Komma trennen)",
            value="80, 85, 100",
            help="Beispieleingabe: 80, 85, 100"
        )

        submitted = st.form_submit_button("Berechnen und Bericht erstellen")

    if submitted:
        try:
            # Abstände parsen
            raw_parts = [p.strip() for p in distances_str.replace(";", ",").split(",") if p.strip() != ""]
            distances = [float(p.replace(",", ".")) for p in raw_parts]
            distances = [d for d in distances if d > 0]

            if not distances:
                st.error("Bitte geben Sie mindestens einen positiven Abstand in Metern an.")
                return

            # Geometrie & Druck aufbereiten
            area = diameter_to_area(diameter, d_unit)
            p_abs = pressure_to_absolute(p_value, p_kind)

            # Physikalische Berechnung
            results = compute_flows_and_acoustics(p_abs, area, distances)

            # Kurze tabellarische Übersicht
            st.subheader("Berechnung – Übersicht")
            col_a, col_b, col_c = st.columns(3)

            with col_a:
                st.markdown("**Gasstrom**")
                st.write(f"Massenstrom: {results['m_dot']:.1f} kg/s")
                st.write(f"Volumenstrom (Betrieb): {results['V_dot_crit_m3_h']:.0f} m³/h")
                st.write(f"Normvolumenstrom: {results['V_dot_N_m3_h']:.0f} Nm³/h")

            with col_b:
                st.markdown("**Akustik**")
                st.write(f"Schallleistung: {results['P_ac']:.1f} W")
                st.write(f"Schallleistungspegel: {results['L_W']:.0f} dB")

            with col_c:
                st.markdown("**Schalldruckpegel (Freifeld)**")
                for r in sorted(results["L_p"].keys()):
                    st.write(f"{r:.1f} m: {results['L_p'][r]:.0f} dB")

            # Berichtstext generieren
            report_text = generate_report(
                p_input=p_value,
                p_kind=p_kind,
                diameter=diameter,
                d_unit=d_unit,
                distances_m=distances,
                results=results,
            )

            st.subheader("Berichtstext")
            st.text_area(
                "Der folgende Text kann direkt in ein Gutachten übernommen werden:",
                value=report_text,
                height=500,
            )

            st.download_button(
                label="Bericht als TXT herunterladen",
                data=report_text.encode("utf-8"),
                file_name="bericht_gasaustritt.txt",
                mime="text/plain",
            )

        except ValueError:
            st.error("Die Abstände konnten nicht interpretiert werden. Bitte geben Sie Zahlen, getrennt durch Kommas, ein.")
        except Exception as e:
            st.error(f"Bei der Berechnung ist ein Fehler aufgetreten: {e}")


if __name__ == "__main__":
    main()
