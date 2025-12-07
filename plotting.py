import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.dates as mdates
from datetime import datetime, timedelta, time as dt_time
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import get_sun, get_body
import pytz
from data_fetching import fetch_stars_vizier
import streamlit.components.v1 as components

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
