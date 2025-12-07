import streamlit as st
import pandas as pd
import numpy as np
from datetime import datetime, time as dt_time, timedelta
import pytz
import re
import astropy.units as u
from astropy.coordinates import EarthLocation, SkyCoord, get_body, AltAz
from astroplan import Observer
from astropy.time import Time
from data_fetching import load_targets, scrape_clear_outside_meta, fetch_weather_api
from plotting import generate_star_chart, plot_yearly_chart, plot_night_chart, render_skymap_decimal
from preferences import load_preferences, add_favorite, remove_favorite, hide_object, get_visible_objects
import streamlit.components.v1 as components

def setup_sidebar():
    with st.sidebar:
        st.header("1. Lokalizacja")
        if 'lat' not in st.session_state:
            st.session_state['lat'] = 52.43
            st.session_state['lon'] = 17.11
            st.session_state['city'] = "Gortatowo"

        city_search_input = st.text_input("Szukaj miasta:", placeholder="Wpisz nazwę...")
        if city_search_input:
            from data_fetching import get_coordinates_from_city
            slat, slon, sname, _ = get_coordinates_from_city(city_search_input)
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

    return lat, lon, fov_w, fov_h, date_obs, min_alt

def display_dashboard(obs, date_obs, lat, lon):
    tz = pytz.timezone('Europe/Warsaw')
    start_of_day = tz.localize(datetime.combine(date_obs, dt_time(0, 0)))
    time_utc = Time(start_of_day)

    s_rise = obs.sun_rise_time(time_utc, which='next').to_datetime(timezone=tz)
    s_set = obs.sun_set_time(time_utc, which='next').to_datetime(timezone=tz)

    astro_dawn = obs.twilight_morning_astronomical(time_utc, which='next').to_datetime(timezone=tz)
    astro_dusk = obs.twilight_evening_astronomical(time_utc, which='next').to_datetime(timezone=tz)

    m_rise = obs.moon_rise_time(time_utc, which='next').to_datetime(timezone=tz)
    m_set = obs.moon_set_time(time_utc, which='next').to_datetime(timezone=tz)

    t_midnight = Time(tz.localize(datetime.combine(date_obs, dt_time(23, 59))))
    moon_illum = obs.moon_illumination(t_midnight) * 100

    yesterday = t_midnight - timedelta(days=1)
    illum_yesterday = obs.moon_illumination(yesterday) * 100
    phase_trend = "Ubywa" if moon_illum < illum_yesterday else "Przybywa"

    moon_phase_angle = obs.moon_phase(t_midnight).to(u.deg).value

    if 0 <= moon_phase_angle < 10 or 350 < moon_phase_angle <= 360:
        p_desc = "Nów"
    elif 170 < moon_phase_angle <= 190:
        p_desc = "Pełnia"
    elif 80 < moon_phase_angle <= 100:
        p_desc = "I Kwadra"
    elif 260 < moon_phase_angle <= 280:
        p_desc = "III Kwadra"
    else:
        p_desc = f"Księżyc ({phase_trend})"

    lp_url = f"https://www.lightpollutionmap.info/#zoom=8&lat={lat}&lon={lon}"
    st.title(f"🔭 RedCat Planner: {st.session_state['city']} [🗺️]({lp_url})")

    with st.expander("Podsumowanie Nocy", expanded=True):
        c1, c2 = st.columns(2)
        with c1:
            st.markdown("#### ☀️ Słońce i Zmierzchy")
            st.write(f"**Wschód Słońca:** {s_rise.strftime('%H:%M')}")
            st.write(f"**Zachód Słońca:** {s_set.strftime('%H:%M')}")
            st.write(f"**Noc astronomiczna:** {astro_dusk.strftime('%H:%M')} - {astro_dawn.strftime('%H:%M')}")

        with c2:
            st.markdown(f"#### 🌕 Księżyc ({p_desc})")
            st.write(f"**Wschód Księżyca:** {m_rise.strftime('%H:%M')}")
            st.write(f"**Zachód Księżyca:** {m_set.strftime('%H:%M')}")
            st.write(f"**Oświetlenie:** {moon_illum:.2f}%")

        st.markdown("[✨ Prognoza Zorzy Polarnej](https://www.spaceweatherlive.com/en/auroral-activity/aurora-forecast.html)", unsafe_allow_html=True)

    st.markdown("#### 🌤️ Prognoza (3 Dni)")
    weather_json = fetch_weather_api(lat, lon)
    df_weather = process_weather_data(weather_json, obs.location)

    if df_weather is not None:
        st.dataframe(
            style_weather(df_weather),
            height=425,
            width='content',
            column_config={
                "_index": st.column_config.Column(width="20"),
            }
        )
    else:
        st.warning("Nie udało się pobrać prognozy pogody.")

def display_object_details(t, obs, fov_w, fov_h, date_obs, moon_coord, i):
    row = t['data']
    coord = t['coord']

    sep = coord.separation(moon_coord).degree
    warn = "⚠️ BLISKO KSIĘŻYCA!" if sep < 30 else ""
    moon_icon = "🌑" if sep > 60 else "🌔"

    if row['mag'] < 90:
        mag_str = f"{row['mag']:.1f}"
    else:
        mag_str = "?"

    name_display = f"**{row['id']}**"
    if row['common_name']:
        name_display += f" – {row['common_name']}"

    label = f"{'⭐' if t['is_favorite'] else ''} {name_display} | Mag: {mag_str} | Alt: {t['alt']:.0f}° | Kadr: {t['fill']:.0f}%"

    with st.expander(label, expanded=(i == 0)):
        c1, c2, c3 = st.columns([1, 1, 1])

        with c1:
            st.markdown("**Symulacja Kadru**")
            st.pyplot(generate_star_chart(row['ra'], row['dec'], fov_w, fov_h, row['size']))

            st.markdown("**Sezonowość**")
            st.pyplot(plot_yearly_chart(obs, coord, pytz.timezone('Europe/Warsaw')))

        with c2:
            st.markdown("**Akcje**")

            if t['is_favorite']:
                if st.button("💔 Usuń z Ulubionych", key=f"fav_{row['id']}"):
                    remove_favorite(row['id'])
                    st.experimental_rerun()
            else:
                if st.button("⭐ Dodaj do Ulubionych", key=f"fav_{row['id']}"):
                    add_favorite(row['id'])
                    st.experimental_rerun()

            if st.button("🙈 Ukryj na 90 dni", key=f"hide_{row['id']}"):
                hide_object(row['id'])
                st.experimental_rerun()

        with c3:
            st.markdown("**Wykres Nocy (dziś)**")
            fig_night = plot_night_chart(obs, coord, date_obs, pytz.timezone('Europe/Warsaw'))
            if fig_night:
                st.pyplot(fig_night)
            else:
                st.warning("Obiekt widoczny zbyt krótko.")

            st.markdown("**Info:**")

            type_map_pl = {'Gx': 'Galaktyka', 'Nb': 'Mgławica', 'Pl': 'Mgławica Planetarna', 'OC': 'Gromada Otwarta', 'Gb': 'Gromada Kulista', 'C+N': 'Gromada + Mgławica', '?': 'Nieznany'}
            obj_type_pl = type_map_pl.get(row['type'], row['type'])
            st.write(f"Typ: **{obj_type_pl}**")

            st.write(f"RA: `{coord.ra.to_string(u.hour, sep=':', precision=0)}`")
            st.write(f"DEC: `{coord.dec.to_string(u.deg, sep=':', precision=0)}`")

            st.write(f"Mag: **{mag_str}**")
            st.write(f"Kadr: **{t['fill']:.1f}%**")
            st.write(f"Odl. od Księżyca: **{sep:.0f}°** {moon_icon}")

            if warn: st.error(warn)

def display_targets(obs, t_midnight, fov_w, fov_h, date_obs, min_alt):
    df_targets = load_targets()

    if df_targets.empty:
        st.stop()

    df_targets = get_visible_objects(df_targets)

    search_query = st.text_input("Filtruj listę obiektów:", placeholder="np. M31, Heart Nebula...")

    if search_query:
        query = search_query.lower()
        mask = (df_targets['id'].str.lower().str.contains(query)) | \
               (df_targets['common_name'].str.lower().str.contains(query))
        df_targets = df_targets[mask]

    valid_targets = []

    with st.spinner(f"Analiza {len(df_targets)} obiektów dla Twojej lokalizacji..."):
        dec_limit = obs.location.lat.deg - 90 + min_alt - 5
        df_filtered = df_targets[df_targets['dec'] > dec_limit].copy()

        for _, row in df_filtered.iterrows():
            try:
                coord = SkyCoord(ra=row['ra'], dec=row['dec'], unit=(u.deg, u.deg))
                alt_mid = obs.altaz(t_midnight, coord).alt.degree

                if not search_query and alt_mid < min_alt:
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
                    "data": row, "coord": coord, "alt": alt_mid,
                    "fill": fill_pct, "score": score,
                    "is_favorite": row['id'] in load_preferences()['favorites']
                })
            except Exception:
                continue

    valid_targets.sort(key=lambda x: (x['is_favorite'], x['score']), reverse=True)

    if not search_query:
        top_targets = valid_targets[:15]
        st.subheader(f"Top 15 obiektów na dziś:")
    else:
        top_targets = valid_targets
        st.subheader(f"Znaleziono {len(top_targets)} obiektów:")

    moon_coord = get_body("moon", t_midnight, obs.location)

    for i, t in enumerate(top_targets):
        display_object_details(t, obs, fov_w, fov_h, date_obs, moon_coord, i)

def process_weather_data(data, location):
    if not data or "hourly" not in data: return None
    h = data['hourly']

    times = pd.to_datetime(h['time'])

    now = pd.Timestamp.now() - pd.Timedelta(hours=1)
    mask = times >= now

    times = times[mask]
    temps = np.array(h['temperature_2m'])[mask]
    feels = np.array(h['apparent_temperature'])[mask]
    rains = np.array(h['precipitation_probability'])[mask]
    dews = np.array(h['dewpoint_2m'])[mask]
    clouds = np.round(np.array(h['cloud_cover'])[mask]).astype(int)
    lows = np.round(np.array(h['cloud_cover_low'])[mask]).astype(int)
    mids = np.round(np.array(h['cloud_cover_mid'])[mask]).astype(int)
    highs = np.round(np.array(h['cloud_cover_high'])[mask]).astype(int)
    vis = [v/1000 for v in np.array(h['visibility'])[mask]]
    winds = np.array(h['wind_speed_10m'])[mask]
    wind_dirs = np.array(h['wind_direction_10m'])[mask]

    sun_display, labels, wind_fmt = [], [], []
    day_map = {0:'Pn', 1:'Wt', 2:'Sr', 3:'Cz', 4:'Pt', 5:'Sb', 6:'Nd'}
    prev_sun_alt = None

    for i, t in enumerate(times):
        t_aware = t.tz_localize(pytz.timezone('Europe/Warsaw')) if t.tzinfo is None else t
        t_utc = Time(t_aware.tz_convert('UTC').to_pydatetime())
        sun_alt = get_body('sun', t_utc, location).transform_to(AltAz(obstime=t_utc, location=location)).alt.degree

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

    final_df = pd.DataFrame({
        "Ciemność": sun_display,
        "Chmury (%)": df_to_str(clouds),
        "Niskie (%)": df_to_str(lows),
        "Średnie (%)": df_to_str(mids),
        "Wysokie (%)": df_to_str(highs),
        "Deszcz (%)": df_to_str(rains),
        "Mgła / Wid. (km)": [f"{x:.1f}" for x in vis],
        "Wiatr": wind_fmt,
        "Temp (°C)": [f"{x:.1f}" for x in temps],
        "Odczuwalna (°C)": [f"{x:.1f}" for x in feels],
        "Rosa (°C)": [f"{x:.1f}" for x in dews]
    })

    final_df.index = labels
    return final_df.transpose()

def df_to_str(arr):
    return [str(x) for x in arr]

def style_weather(df):
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
