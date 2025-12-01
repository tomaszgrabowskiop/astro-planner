import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.dates as mdates
import requests
import warnings
import re
from datetime import datetime, timedelta, time as dt_time
from bs4 import BeautifulSoup
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, get_sun, get_body
from astroplan import Observer
from astroquery.vizier import Vizier
from astropy.utils.exceptions import AstropyWarning
import pytz
import streamlit.components.v1 as components

# ==========================================
# 0. KONFIGURACJA
# ==========================================
st.set_page_config(page_title="RedCat Planner Ultimate", layout="wide", page_icon="🔭")
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.filterwarnings("ignore")
plt.style.use('dark_background')

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
# 1. FUNKCJE DANYCH
# ==========================================

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

def process_weather_data(data, location):
    if not data or "hourly" not in data: return None
    h = data['hourly']
    
    # 1. PRZYGOTOWANIE DANYCH
    times = pd.to_datetime(h['time'])
    
    # Filtrowanie przeszłości
    now = pd.Timestamp.now() - pd.Timedelta(hours=1)
    mask = times >= now
    
    # Przycinanie tablic
    times = times[mask]
    temps = np.array(h['temperature_2m'])[mask]
    feels = np.array(h['apparent_temperature'])[mask] # <--- Pobranie odczuwalnej
    rains = np.array(h['precipitation_probability'])[mask]
    dews = np.array(h['dewpoint_2m'])[mask]
    clouds = np.round(np.array(h['cloud_cover'])[mask]).astype(int)
    lows = np.round(np.array(h['cloud_cover_low'])[mask]).astype(int)
    mids = np.round(np.array(h['cloud_cover_mid'])[mask]).astype(int)
    highs = np.round(np.array(h['cloud_cover_high'])[mask]).astype(int)
    # Konwersja widoczności z metrów na kilometry
    vis = [v/1000 for v in np.array(h['visibility'])[mask]] 
    winds = np.array(h['wind_speed_10m'])[mask]
    wind_dirs = np.array(h['wind_direction_10m'])[mask]
    
    sun_display, labels, wind_fmt = [], [], []
    day_map = {0:'Pn', 1:'Wt', 2:'Sr', 3:'Cz', 4:'Pt', 5:'Sb', 6:'Nd'}
    prev_sun_alt = None
    
    # Pętla obliczająca ciemność i etykiety (bez zmian)
    for i, t in enumerate(times):
        t_aware = t.tz_localize(TZ_PL) if t.tzinfo is None else t
        t_utc = Time(t_aware.tz_convert('UTC').to_pydatetime())
        sun_alt = get_sun(t_utc).transform_to(AltAz(obstime=t_utc, location=location)).alt.degree
        
        if sun_alt > -0.833: s_txt = "Dzien"
        elif sun_alt > -6: s_txt = "Cywilny"
        elif sun_alt > -12: s_txt = "Zeglarski"
        elif sun_alt > -18: s_txt = "Astro"
        else: s_txt = "Noc"
        
        if prev_sun_alt is not None:
            h_geom = -0.833
            if prev_sun_alt > h_geom and sun_alt < h_geom:
                mins = int(((prev_sun_alt - h_geom) / (prev_sun_alt - sun_alt)) * 60)
                s_txt = f"↓ {t.hour:02d}:{mins:02d}"
            elif prev_sun_alt < h_geom and sun_alt > h_geom:
                mins = int(((h_geom - prev_sun_alt) / (sun_alt - prev_sun_alt)) * 60)
                s_txt = f"↑ {t.hour:02d}:{mins:02d}"
        
        prev_sun_alt = sun_alt
        sun_display.append(s_txt)
        d_name = day_map[t.dayofweek]
        labels.append(f"{d_name} {t.strftime('%H')}")
        
        dirs = ["N", "NE", "E", "SE", "S", "SW", "W", "NW", "N"]
        w_dir = dirs[int(round(wind_dirs[i] / 45)) % 8]
        wind_fmt.append(f"{w_dir} {int(winds[i])}")

    # 2. BUDOWANIE TABELI KOŃCOWEJ (ZMIANA KOLEJNOŚCI)
    # Tutaj ustalamy co jest pod czym
    final_df = pd.DataFrame({
        "Ciemność": sun_display,
        "Chmury (%)": df_to_str(clouds),
        "Niskie (%)": df_to_str(lows),
        "Średnie (%)": df_to_str(mids),
        "Wysokie (%)": df_to_str(highs),
        
        # SEKCJ 1: WODA
        "Deszcz (%)": df_to_str(rains),
        "Mgła / Wid. (km)": [f"{x:.1f}" for x in vis], # <--- MGŁA POD DESZCZEM
        
        "Wiatr": wind_fmt,
        
        # SEKCJA 2: TEMPERATURA
        "Temp (°C)": [f"{x:.1f}" for x in temps],
        "Odczuwalna (°C)": [f"{x:.1f}" for x in feels], # <--- ODCZUWALNA POD TEMP
        "Rosa (°C)": [f"{x:.1f}" for x in dews]
    })
    
    final_df.index = labels
    return final_df.transpose()

# Pomocnicza funkcja, żeby kod był czystszy
def df_to_str(arr):
    return [str(x) for x in arr]

def style_weather(df):
    # Ponieważ dane są teraz stringami, musimy rzutować na int/float wewnątrz map
    def bg_clouds(v):
        try:
            val = int(float(v))
            if val <= 10: return 'background-color: #004d00; color: #cfc'
            if val <= 40: return 'background-color: #8B8000; color: white'
            return 'background-color: #8b0000; color: white; font-weight: bold'
        except: return ''
    
    def bg_darkness(v):
        s = str(v)
        if "↓" in s or "↑" in s: return 'background-color: #d35400; color: white; font-weight: bold'
        colors = {"Cywilny": "#a0522d", "Zeglarski": "#483d8b", "Astro": "#191970", "Noc": "#000000", "Dzien": "#4682b4"}
        return f'background-color: {colors.get(s, "black")}; color: transparent'
    
    def bg_rain(v):
        try:
            val = int(float(v))
            if val > 50: return 'background-color: #00008b; color: white; font-weight: bold'
            if val > 20: return 'background-color: #4169e1; color: white'
            return ''
        except: return ''

    styler = df.style
    idx = pd.IndexSlice
    
    styler.map(bg_darkness, subset=idx[["Ciemność"], :])
    styler.map(bg_clouds, subset=idx[["Chmury (%)", "Niskie (%)", "Średnie (%)", "Wysokie (%)"], :])
    styler.map(bg_rain, subset=idx[["Deszcz (%)"], :])
    styler.map(lambda v: 'background-color: #1a2533; color: #b0c4de', subset=idx[["Wiatr", "Temp (°C)", "Odczuwalna (°C)", "Rosa (°C)"], :])
    return styler

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

# ==========================================
# 2. WYKRESY
# ==========================================

def generate_star_chart(ra_center, dec_center, fov_w, fov_h, obj_size_arcmin):
    # 1. Zabezpieczenie rozmiaru (żeby nie było błędu NaN)
    if obj_size_arcmin is None or np.isnan(obj_size_arcmin) or obj_size_arcmin <= 0:
        obj_size_arcmin = 1.0

    # 2. Pobieramy gwiazdy (bierzemy nieco większy margines, żeby złapać rogi przy rotacji sfery)
    max_fov = max(fov_w, fov_h)
    stars = fetch_stars_vizier(ra_center, dec_center, max_fov * 1.5)
    
    fig, ax = plt.subplots(figsize=(4, 4))
    fig.patch.set_facecolor('#0e0e0e')
    ax.set_facecolor('#000000') 
    
    # 3. PROJEKCJA SFERYCZNA NA PŁASZCZYZNĘ (To naprawia "rozjeżdżanie")
    if stars is not None and not stars.empty:
        try:
            col_ra = next((c for c in stars.columns if '_RA' in c or 'RA' in c), None)
            col_de = next((c for c in stars.columns if '_DE' in c or 'DE' in c), None)
            col_mag = next((c for c in stars.columns if 'VTmag' in c or 'mag' in c), None)
            
            if col_ra and col_de:
                st_ra = stars[col_ra].values
                st_dec = stars[col_de].values
                
                # A. Naprawa przejścia przez 0h (RA Wrapping)
                # Obliczamy różnicę RA i normalizujemy ją do zakresu -180...180
                delta_ra = st_ra - ra_center
                delta_ra = (delta_ra + 180) % 360 - 180 

                # B. Projekcja (Przybliżenie płaszczyzny stycznej)
                # Skracamy oś X w zależności od deklinacji (im wyżej, tym bardziej skracamy)
                x_proj = delta_ra * np.cos(np.radians(dec_center))
                y_proj = st_dec - dec_center

                # C. Rysowanie gwiazd
                mags = stars[col_mag].fillna(10).values if col_mag else np.full(len(stars), 10)
                sizes = (11 - mags) * 4
                sizes = np.clip(sizes, 0.5, 40)
                
                # Oś X odwrócona (-x_proj), bo na mapach nieba Wschód jest z lewej
                ax.scatter(-x_proj, y_proj, s=sizes, c='white', alpha=0.9, edgecolors='none')
        except Exception: 
            pass

    # 4. Rysowanie elementów (Teraz wszystko jest wycentrowane w 0,0)
    
    # Ramka kadru (FOV)
    rect = patches.Rectangle((-fov_w/2, -fov_h/2), fov_w, fov_h, 
                             linewidth=1.5, edgecolor='#ff3333', facecolor='none')
    ax.add_patch(rect)
    
    # Kółko oznaczające obiekt
    r_deg = (obj_size_arcmin / 60.0) / 2
    dso_circle = patches.Ellipse((0, 0), width=r_deg*2, height=r_deg*2, 
                                 linewidth=1, edgecolor='#00ccff', facecolor='#00ccff', alpha=0.2)
    ax.add_patch(dso_circle)
    
    # Krzyżyk na środku
    ax.plot(0, 0, '+', color='#00ccff', markersize=8)

    # 5. Ustawienie limitów i proporcji
    margin = 1.2
    limit_w = (fov_w / 2) * margin
    limit_h = (fov_h / 2) * margin
    
    # Limity są teraz symetryczne względem 0
    ax.set_xlim(limit_w, -limit_w) # Odwracamy X (konwencja astronomiczna)
    ax.set_ylim(-limit_h, limit_h)
    
    ax.set_aspect('equal') # KLUCZOWE: Teraz piksele są kwadratowe, bo zrobiliśmy projekcję ręcznie
    ax.axis('off')
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    return fig

def plot_yearly_chart(obs, target_coord, tz_pl):
    now = datetime.now()
    current_year = now.year
    
    # 1. Obliczamy dokładną pozycję "Teraz" na osi X
    # Wykres ma punkty co miesiąc (1, 2, 3...). Zakładamy, że punkt to środek miesiąca (15-ty).
    # Przesuwamy linię w zależności od dnia: (dzień - 15) / 30
    # Np. 1 grudnia = 12 + (-0.5) = 11.5 (połowa między Listopadem a Grudniem)
    current_x = now.month + (now.day - 15) / 30.0
    
    months = []
    alts = []
    for m in range(1, 13):
        d = datetime(current_year, m, 15, 0, 0)
        dt_local = tz_pl.localize(d)
        t_utc = Time(dt_local.astimezone(pytz.utc))
        alt = obs.altaz(t_utc, target_coord).alt.degree
        alts.append(max(0, alt))
        months.append(m)

    fig, ax = plt.subplots(figsize=(4, 1.5))
    fig.patch.set_facecolor('#0e0e0e')
    ax.set_facecolor('#0e0e0e')
    
    ax.plot(months, alts, color='#00ccff', linewidth=1.5)
    ax.fill_between(months, 0, alts, color='#00ccff', alpha=0.2)
    
    # 2. Rysujemy czerwoną linię
    ax.axvline(current_x, color='#ff3333', linestyle='--', linewidth=1.5)
    
    # 3. NAPRAWA: Zwiększamy zakres X (od 0.5 do 12.5)
    # Dzięki temu linie dla Stycznia (1) i Grudnia (12) nie będą ucinane przez krawędź okna
    ax.set_xlim(0.5, 12.5)
    
    ax.set_ylim(0, 90)
    ax.axis('off')
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    return fig

def plot_night_chart(obs, target_coord, date_obs, tz_pl):
    t_noon = Time(datetime.combine(date_obs, dt_time(12, 0)))
    try:
        t_start = obs.twilight_evening_astronomical(t_noon, which='next')
        t_end = obs.twilight_morning_astronomical(t_noon, which='next')
    except:
        try:
            t_start = obs.twilight_evening_nautical(t_noon, which='next')
            t_end = obs.twilight_morning_nautical(t_noon, which='next')
        except: return None
        
    if t_end < t_start: t_end = t_end + 1*u.day
    
    start_dt, end_dt = t_start.to_datetime(), t_end.to_datetime()
    minutes = int((end_dt - start_dt).total_seconds() / 60)
    if minutes <= 0: return None
    
    time_points = [start_dt + timedelta(minutes=i) for i in range(0, minutes, 15)]
    astro_times = Time(time_points)
    alt_target = obs.altaz(astro_times, target_coord).alt.degree
    
    moon_coord = get_body("moon", astro_times, obs.location)
    alt_moon = obs.altaz(astro_times, moon_coord).alt.degree
    
    times_pl = [t.replace(tzinfo=pytz.utc).astimezone(tz_pl) for t in time_points]

    fig, ax = plt.subplots(figsize=(5, 2))
    fig.patch.set_facecolor('#0e0e0e')
    ax.set_facecolor('#0e0e0e')
    
    ax.plot(times_pl, alt_target, color='#00ccff', linewidth=2)
    ax.fill_between(times_pl, 0, alt_target, color='#00ccff', alpha=0.1)
    
    moon_vis = np.maximum(0, alt_moon)
    if np.max(moon_vis) > 0:
        ax.fill_between(times_pl, 0, moon_vis, color='#aaaaaa', alpha=0.2)
        ax.plot(times_pl, moon_vis, color='#aaaaaa', linestyle='--', linewidth=1)

    ax.axhline(30, color='#444', linestyle=':')
    ax.set_ylim(0, 90)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M', tz=tz_pl))
    ax.tick_params(colors='#888', labelsize=7)
    ax.grid(color='#222', linestyle=':', linewidth=0.5)
    
    return fig

def render_skymap_decimal(ra_deg, dec_deg):
    ra_hours = ra_deg / 15.0
    url = f"https://server1.sky-map.org/skywindow?ra={ra_hours}&de={dec_deg}&zoom=4&img_source=DSS2"
    components.iframe(url, height=350, scrolling=False)
    return url

# ==========================================
# 3. GŁÓWNA APLIKACJA
# ==========================================

def main():
    with st.sidebar:
        st.header("1. Lokalizacja")
        if 'lat' not in st.session_state: 
            st.session_state['lat'] = 52.43
            st.session_state['lon'] = 17.11
            st.session_state['city'] = "Gortatowo"

        city_search = st.text_input("Szukaj miasta:", placeholder="Wpisz nazwę...")
        if city_search:
            slat, slon, sname, _ = get_coordinates_from_city(city_search)
            if slat:
                st.session_state['lat'] = slat
                st.session_state['lon'] = slon
                st.session_state['city'] = sname
                st.success(f"Ustawiono: {sname}")
        
        st.caption(f"📍 {st.session_state['city']} ({st.session_state['lat']:.2f}, {st.session_state['lon']:.2f})")
        lat = float(st.session_state['lat'])
        lon = float(st.session_state['lon'])

        st.header("2. Teleskop / Kadr")
        fl = st.number_input("Ogniskowa (mm)", value=300)
        sw = st.number_input("Sensor W (mm)", value=23.5)
        sh = st.number_input("Sensor H (mm)", value=15.7)
        fov_w = (2 * np.arctan((sw / (2 * fl)))) * (180 / np.pi)
        fov_h = (2 * np.arctan((sh / (2 * fl)))) * (180 / np.pi)
        st.info(f"FOV: {fov_w:.1f}° x {fov_h:.1f}°")
        
        st.header("3. Filtracja")
        date_obs = st.date_input("Data obserwacji", pd.Timestamp.now())
        min_alt = st.slider("Min. Wysokość o północy", 10, 80, 25)

    # --- SETUP ASTRONOMICZNY ---
    loc = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=100*u.m)
    obs = Observer(location=loc)
    
    t_noon = Time(datetime.combine(date_obs, dt_time(12, 0)))
    t_midnight = t_noon + 12*u.hour

    st.title(f"🔭 RedCat Planner: {st.session_state['city']}")

    # --- DASHBOARD POGODOWY (SEKCJA ROZBUDOWANA) ---
    def fmt_time(t_obj):
        if hasattr(t_obj, 'to_datetime'):
            dt = t_obj.to_datetime()
            if dt.tzinfo is None: dt = dt.replace(tzinfo=pytz.utc)
            return dt.astimezone(TZ_PL).strftime('%H:%M')
        return "-"

    s_rise = fmt_time(obs.sun_rise_time(t_noon, which='nearest'))
    s_set = fmt_time(obs.sun_set_time(t_noon, which='nearest'))
    dark_start = fmt_time(obs.twilight_evening_astronomical(t_noon, which='nearest'))
    dark_end = fmt_time(obs.twilight_morning_astronomical(t_noon, which='next'))

    try:
        gh_start = fmt_time(obs.sun_set_time(t_noon, which='nearest', horizon=6*u.deg))
        gh_end = fmt_time(obs.sun_set_time(t_noon, which='nearest', horizon=-4*u.deg))
        golden_str = f"{gh_start} - {gh_end}"
    except: golden_str = "-"

    moon_illum = obs.moon_illumination(t_midnight) * 100
    moon_phase = obs.moon_phase(t_midnight)
    deg = moon_phase.to(u.deg).value
    m_rise = fmt_time(obs.moon_rise_time(t_midnight, which='next'))
    m_set = fmt_time(obs.moon_set_time(t_midnight, which='next'))

    if 350 < deg or deg <= 10: p_desc = "Nów"
    elif 10 < deg <= 80: p_desc = "Sierp (Rosnący)"
    elif 80 < deg <= 100: p_desc = "I Kwadra"
    elif 100 < deg <= 170: p_desc = "Garb (Rosnący)"
    elif 170 < deg <= 190: p_desc = "Pełnia"
    elif 190 < deg <= 260: p_desc = "Garb (Malejący)"
    elif 260 < deg <= 280: p_desc = "III Kwadra"
    else: p_desc = "Sierp (Malejący)"

    with st.spinner("Pobieranie danych o niebie..."):
        meta = scrape_clear_outside_meta(lat, lon)

  # --- WIZUALIZACJA DASHBOARDU ---
    st.markdown("### 🌌 Warunki i Pogoda")
    
    # 1. Metryki Jakości Nieba
    c_sqm1, c_sqm2, c_sqm3, c_sqm4,  = st.columns(4)
    c_sqm1.metric("SQM", meta.get("sqm", "-"))
    c_sqm2.metric("Bortle", meta.get("bortle", "-"))
    c_sqm3.metric("Jasność", meta.get("brightness", "-"))
    # ZMIANA: Usunięto p_desc (opis słowny) stąd, zostaje sam procent
    c_sqm4.metric("Księżyc", f"{moon_illum:.0f}%")

    # 2. Słońce i Księżyc
    c_sun, c_moon = st.columns(2)
    with c_sun:
        st.markdown("#### ☀️ Słońce")
        cc1, cc2 = st.columns(2)
        cc1.metric("Wschód", s_rise)
        cc2.metric("Zachód", s_set)
        st.info(f"📸 Złota Godzina: **{golden_str}**")
        st.markdown("##### 🌌 Noc Astronomiczna")
        d1, d2 = st.columns(2)
        d1.metric("Start", dark_start)
        d2.metric("Koniec", dark_end)
        

    with c_moon:
        st.markdown("#### 🌕 Księżyc")
        mm1, mm2 = st.columns(2)
        mm1.metric("Wschód", m_rise)
        mm2.metric("Zachód", m_set)
        st.success(f"🌑 Faza: **{p_desc}** | Geometria: **{deg:.1f}°**")

    st.markdown("#### 🌤️ Prognoza (3 Dni)")
    weather_json = fetch_weather_api(lat, lon)
    df_weather = process_weather_data(weather_json, loc)
    
    if df_weather is not None:
        st.dataframe(
            style_weather(df_weather), 
            height=425,                 # Wysokość tabeli (żeby nie było paska przewijania pionowego)
            width='content',   # Rozciąga tabelę na całą szerokość strony (żeby nie była wąska). For `use_container_width=True`, use `width='stretch'`. For `use_container_width=False`, use `width='content'`.
            column_config={
                "_index": st.column_config.Column(width="20"), # "small", "medium", "large" lub liczba np. 150
            }
        )
    else:
        st.warning("Nie udało się pobrać prognozy pogody.")
        
    st.markdown("---")

    # --- LOGIKA BIZNESOWA (TARGETS) ---
    df_targets = load_targets()
    
    if df_targets.empty:
        st.stop()

    valid_targets = []
    
    with st.spinner(f"Analiza {len(df_targets)} obiektów dla Twojej lokalizacji..."):
        dec_limit = lat - 90 + min_alt - 5 
        df_filtered = df_targets[df_targets['dec'] > dec_limit].copy()

        for _, row in df_filtered.iterrows():
            try:
                coord = SkyCoord(ra=row['ra'], dec=row['dec'], unit=(u.deg, u.deg))
                alt_mid = obs.altaz(t_midnight, coord).alt.degree
                
                if alt_mid < min_alt:
                    continue
                
                obj_area = (row['size'] / 60.0) ** 2
                fov_area = fov_w * fov_h
                fill_pct = min(100, (obj_area / fov_area) * 100)
                
                score = fill_pct * 1.5 
                mag = row['mag']
                if mag < 90:
                    score += (15 - mag) * 5
                else:
                    score += 20 

                valid_targets.append({
                    "data": row,
                    "coord": coord,
                    "alt": alt_mid,
                    "fill": fill_pct,
                    "score": score
                })
            except Exception:
                continue

    valid_targets.sort(key=lambda x: x['score'], reverse=True)
    top_targets = valid_targets[:15]

    st.subheader(f"Znaleziono {len(valid_targets)} obiektów. Oto Top 15 na dziś:")

    moon_coord = get_body("moon", t_midnight, loc)

    type_map_pl = {
        'Gx': 'Galaktyka', 'Nb': 'Mgławica', 'Pl': 'Mgławica Planetarna',
        'OC': 'Gromada Otwarta', 'Gb': 'Gromada Kulista', 'C+N': 'Gromada + Mgławica', '?': 'Nieznany'
    }

    for i, t in enumerate(top_targets):
        row = t['data']
        coord = t['coord']
        
        # Logika Księżyca
        sep = coord.separation(moon_coord).degree
        warn = "⚠️ BLISKO KSIĘŻYCA!" if sep < 30 else ""
        moon_icon = "🌑" if sep > 60 else "🌔"
        
        # Formatowanie jasności
        if row['mag'] < 90:
            mag_str = f"{row['mag']:.1f}"
        else:
            mag_str = "?"
        
        # --- NOWA LOGIKA ETYKIETY (ID + NAZWA) ---
        name_display = f"**{row['id']}**"
        
        # Jeśli mamy nazwę zwyczajową (np. Heart Nebula), dopisujemy ją
        if row['common_name']:
            name_display += f" – {row['common_name']}"
            
        label = f"{name_display} | Mag: {mag_str} | Alt: {t['alt']:.0f}° | Kadr: {t['fill']:.0f}%"
        
        with st.expander(label, expanded=(i == 0)):
            c1, c2, c3 = st.columns([1, 1, 1])
            
            with c1:
                st.markdown("**Symulacja Kadru**")
                st.pyplot(generate_star_chart(row['ra'], row['dec'], fov_w, fov_h, row['size']))
                
                st.markdown("**Sezonowość**")
                st.pyplot(plot_yearly_chart(obs, coord, TZ_PL))

            with c2:
                st.markdown("**Podgląd DSS2**")
                url_dss = render_skymap_decimal(row['ra'], row['dec'])
                # --- PRZYCISK 1: Pełna Mapa (NIEBIESKI) ---
                st.markdown(
                    f"""
                    <a href="{url_dss}" target="_blank" style="text-decoration: none;">
                        <div style="
                            background-color: #1c83e1; 
                            color: white; 
                            padding: 8px; 
                            border-radius: 5px; 
                            text-align: center; 
                            font-weight: bold;
                            margin-top: 10px;
                            margin-bottom: 5px;
                            border: 1px solid #1668b2a;
                        ">
                            🗺️ Pełna Mapa (DSS2)
                        </div>
                    </a>
                    """,
                    unsafe_allow_html=True
                )
                # --- PRZYCISK 2: Telescopius (ZIELONY) ---
                # 1. Przygotowanie ID dla URL (np. "NGC 7000" -> "ngc-7000", "M42" -> "m-42")
                tid = row['id'].lower().strip()
                # FIX: Usunięcie końcówki ".0" (np. "sh2-260.0" -> "sh2-260")
                tid = tid.replace(".0", "")
                # Jeśli to Messier bez myślnika (np. "m42"), dodaj myślnik ("m-42")
                if re.match(r"^m\d+$", tid):
                    tid = tid.replace("m", "m-")
                # Zamień spacje na myślniki (dla NGC/IC/Sh2)
                tid = tid.replace(" ", "-")
                
                tele_url = f"https://telescopius.com/deep-sky-objects/{tid}"
                
                # 2. Wyświetlenie jako zielony przycisk (HTML)
                st.markdown(
                    f"""
                    <a href="{tele_url}" target="_blank" style="text-decoration: none;">
                        <div style="
                            background-color: #2eb050; 
                            color: white; 
                            padding: 8px; 
                            border-radius: 5px; 
                            text-align: center; 
                            font-weight: bold;
                            margin-top: 5px;
                            border: 1px solid #248f41;
                        ">
                            🔭 Zobacz na Telescopius
                        </div>
                    </a>
                    """,
                    unsafe_allow_html=True
                )
            with c3:
                st.markdown("**Wykres Nocy (dziś)**")
                fig_night = plot_night_chart(obs, coord, date_obs, TZ_PL)
                if fig_night:
                    st.pyplot(fig_night)
                else:
                    st.warning("Obiekt widoczny zbyt krótko.")
                
                st.markdown("---")
                st.markdown("**Info:**")
                
                obj_type_pl = type_map_pl.get(row['type'], row['type'])
                st.write(f"Typ: **{obj_type_pl}**")
                
                st.write(f"RA: `{coord.ra.to_string(u.hour, sep=':', precision=0)}`")
                st.write(f"DEC: `{coord.dec.to_string(u.deg, sep=':', precision=0)}`")
                
                st.write(f"Mag: **{mag_str}**")
                st.write(f"Kadr: **{t['fill']:.1f}%**")
                st.write(f"Odl. od Księżyca: **{sep:.0f}°** {moon_icon}")
                
                if warn: st.error(warn)

if __name__ == "__main__":
    main()