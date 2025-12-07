import streamlit as st
import warnings
import pytz
from datetime import datetime, time as dt_time
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astroplan import Observer
from astropy.utils.exceptions import AstropyWarning

from ui import setup_sidebar, display_dashboard, display_targets, display_search_results, display_external_maps

# ==========================================
# 0. KONFIGURACJA
# ==========================================
st.set_page_config(page_title="RedCat Planner Ultimate", layout="wide", page_icon="🔭")
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.filterwarnings("ignore")

TZ_PL = pytz.timezone('Europe/Warsaw')

st.markdown("""
    <style>
    .stApp { background-color: #050505; color: #b0b0b0; }
    
    /* KOLOR I ROZMIAR NAGŁÓWKÓW (SQM, Bortle...) */
    .stMetricLabel { 
        color: #e0e0e0 !important; 
        font-size: 14px !important; /* <-- TU ZMIEŃ wielkość napisu "SQM" */
    }
    
    /* KOLOR I ROZMIAR WARTOŚCI (21.50, 4...) */
    .stMetricValue { 
        color: #ffffff !important; 
        font-size: 1.8rem !important; /* <-- TU ZMIEŃ wielkość liczby (1.1rem to mało, 1.8rem to dużo) */
    }
    
    /* WYGLĄD SAMEGO BOKSU (Ramka) */
    div[data-testid="stMetric"] { 
        background-color: #111; 
        border: 1px solid #333; 
        padding: 15px; /* <-- TU ZMIEŃ margines wewnętrzny (zwiększy/zmniejszy boks) */
        border-radius: 10px; /* <-- Zaokrąglenie rogów */
    }
    
    div[data-testid="stExpander"] { background-color: #0e0e0e; border: 1px solid #333; }
    iframe { border: 1px solid #333; }
    </style>
""", unsafe_allow_html=True)

# ==========================================
# 1. GŁÓWNA APLIKACJA
# ==========================================

def main():
    lat, lon, fov_w, fov_h, date_obs, min_alt, search_query = setup_sidebar()

    # --- SETUP ASTRONOMICZNY ---
    loc = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=100*u.m)
    obs = Observer(location=loc)
    
    t_noon = Time(datetime.combine(date_obs, dt_time(12, 0)))
    t_midnight = t_noon + 12*u.hour

    if search_query:
        display_search_results(search_query)
    else:
        display_dashboard(obs, t_noon, t_midnight, lat, lon)
        display_targets(obs, t_midnight, fov_w, fov_h, date_obs, min_alt)
        
    display_external_maps(lat, lon)

if __name__ == "__main__":
    main()
