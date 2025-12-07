import streamlit as st
import pandas as pd
import requests
import re
from bs4 import BeautifulSoup
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u

@st.cache_data
def load_targets():
    """Wczytuje bazę obiektów z pliku CSV."""
    try:
        df = pd.read_csv("targets.csv")

        # Naprawa jasności (NaN -> 99.9)
        df['mag'] = pd.to_numeric(df['mag'], errors='coerce').fillna(99.9)

        # Naprawa rozmiaru (NaN -> 0.1). To zapobiegnie błędowi Matplotlib!
        df['size'] = pd.to_numeric(df['size'], errors='coerce').fillna(0.1)

        # --- NOWOŚĆ: Obsługa nazwy zwyczajowej ---
        # Jeśli kolumny nie ma (stary plik), tworzymy ją pustą.
        # Jeśli jest, zamieniamy puste wartości (NaN) na pusty tekst "".
        if 'common_name' not in df.columns:
            df['common_name'] = ""
        else:
            df['common_name'] = df['common_name'].fillna("")
        # -----------------------------------------

        return df
    except FileNotFoundError:
        st.error("Błąd: Nie znaleziono pliku 'targets.csv'. Uruchom najpierw generator!")
        return pd.DataFrame()

def get_coordinates_from_city(city_name):
    try:
        url = f"https://geocoding-api.open-meteo.com/v1/search?name={city_name}&count=1&language=pl&format=json"
        r = requests.get(url, timeout=5)
        data = r.json()
        if "results" in data and len(data["results"]) > 0:
            res = data["results"][0]
            return res["latitude"], res["longitude"], res["name"], res.get("country", "")
    except: pass
    return None, None, None, None

@st.cache_data(ttl=3600)
def scrape_clear_outside_meta(lat, lon):
    url = f"https://clearoutside.com/forecast/{lat}/{lon}"
    try:
        r = requests.get(url, headers={'User-Agent': 'Mozilla/5.0'}, timeout=10)
        if r.status_code != 200: return {}
        soup = BeautifulSoup(r.content, 'html.parser')
        text = soup.get_text(separator=" ", strip=True)
        meta = {}
        meta["sqm"] = (re.search(r"([\d\.]+)\s*Magnitude", text) or [None, "-"])[1]
        meta["bortle"] = (re.search(r"Class\s*(\d+)\s*Bortle", text) or [None, "-"])[1]
        meta["brightness"] = (re.search(r"([\d\.]+)\s*mcd", text) or [None, "-"])[1]
        return meta
    except: return {}
@st.cache_data(ttl=1800)
def fetch_weather_api(lat, lon):
    try:
        # ZMIANA: Upewniamy się, że w URL są parametry 'apparent_temperature' (odczuwalna) i 'visibility' (mgła)
        url = (
            f"https://api.open-meteo.com/v1/forecast?latitude={lat}&longitude={lon}"
            f"&hourly=temperature_2m,apparent_temperature,dewpoint_2m,precipitation_probability,cloud_cover,cloud_cover_low,cloud_cover_mid,cloud_cover_high,visibility,wind_speed_10m,wind_direction_10m"
            f"&forecast_days=3&timezone=auto"
        )
        return requests.get(url, timeout=5).json()
    except: return None

@st.cache_data(ttl=3600*24)
def fetch_stars_vizier(ra_deg, dec_deg, fov_deg):
    try:
        # Oryginalne zapytanie, które działało
        v = Vizier(columns=['_RA', '_DE', 'VTmag'], column_filters={'VTmag': '<11'}, row_limit=2500)
        coord = SkyCoord(ra=ra_deg, dec=dec_deg, unit=(u.deg, u.deg), frame='icrs')
        radius = (fov_deg / 2) * 1.4
        result = v.query_region(coord, radius=radius*u.deg, catalog='I/259/tyc2')
        if len(result) > 0: return result[0].to_pandas()
    except: pass
    return None
