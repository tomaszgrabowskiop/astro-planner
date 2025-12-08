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
        df['mag'] = pd.to_numeric(df['mag'], errors='coerce').fillna(99.9)
        df['size'] = pd.to_numeric(df['size'], errors='coerce').fillna(0.1)
        if 'common_name' not in df.columns:
            df['common_name'] = ""
        else:
            df['common_name'] = df['common_name'].fillna("")
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
def scrape_lightpollutionmap_data(lat, lon):
    url = f"https://lightpollutionmap.app/pl/#lat={lat}&lon={lon}&zoom=10"
    try:
        # Używamy aiohttp, aby obsłużyć dynamicznie ładowaną treść
        import asyncio
        from pyppeteer import launch

        async def fetch():
            browser = await launch(headless=True, args=['--no-sandbox'])
            page = await browser.newPage()
            await page.goto(url, {'waitUntil' : 'networkidle0'})
            content = await page.content()
            await browser.close()
            return content

        html_content = asyncio.get_event_loop().run_until_complete(fetch())

        soup = BeautifulSoup(html_content, 'html.parser')

        sqm_value = soup.select_one('span.sqm-value').text.strip()
        bortle_value = soup.select_one('span.bortle-value').text.strip()

        return {"sqm": sqm_value, "bortle": bortle_value}
    except Exception as e:
        return {"sqm": "-", "bortle": "-"}


@st.cache_data(ttl=1800)
def fetch_weather_api(lat, lon):
    try:
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
        v = Vizier(columns=['_RA', '_DE', 'VTmag'], column_filters={'VTmag': '<11'}, row_limit=2500)
        coord = SkyCoord(ra=ra_deg, dec=dec_deg, unit=(u.deg, u.deg), frame='icrs')
        radius = (fov_deg / 2) * 1.4
        result = v.query_region(coord, radius=radius*u.deg, catalog='I/259/tyc2')
        if len(result) > 0: return result[0].to_pandas()
    except: pass
    return None
