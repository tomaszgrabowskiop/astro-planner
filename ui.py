import streamlit as st
import pandas as pd
import numpy as np
from datetime import datetime, time as dt_time
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

        st.header("3. Filtracja i Wyszukiwanie")
        date_obs = st.date_input("Data obserwacji", pd.Timestamp.now())
        min_alt = st.slider("Min. Wysokość o północy", 10, 80, 25)

        search_query = st.text_input("Szukaj obiektu w bazie:", placeholder="np. M31, Heart Nebula...")

    return lat, lon, fov_w, fov_h, date_obs, min_alt, search_query

def display_dashboard(obs, t_noon, t_midnight, lat, lon):
    st.title(f"🔭 RedCat Planner: {st.session_state['city']}")

    def fmt_time(t_obj):
        if hasattr(t_obj, 'to_datetime'):
            dt = t_obj.to_datetime()
            if dt.tzinfo is None: dt = dt.replace(tzinfo=pytz.utc)
            return dt.astimezone(pytz.timezone('Europe/Warsaw')).strftime('%H:%M')
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

    st.markdown("### 🌌 Warunki i Pogoda")

    c_sqm1, c_sqm2, c_sqm3, c_sqm4,  = st.columns(4)
    c_sqm1.metric("SQM", meta.get("sqm", "-"))
    c_sqm2.metric("Bortle", meta.get("bortle", "-"))
    c_sqm3.metric("Jasność", meta.get("brightness", "-"))
    c_sqm4.metric("Księżyc", f"{moon_illum:.0f}%")

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

    st.markdown("---")

def display_targets(obs, t_midnight, fov_w, fov_h, date_obs, min_alt):
    df_targets = load_targets()

    if df_targets.empty:
        st.stop()

    df_targets = get_visible_objects(df_targets)

    valid_targets = []

    with st.spinner(f"Analiza {len(df_targets)} obiektów dla Twojej lokalizacji..."):
        dec_limit = obs.location.lat.deg - 90 + min_alt - 5
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
                    "score": score,
                    "is_favorite": row['id'] in load_preferences()['favorites']
                })
            except Exception:
                continue

    valid_targets.sort(key=lambda x: (x['is_favorite'], x['score']), reverse=True)
    top_targets = valid_targets[:15]

    st.subheader(f"Znaleziono {len(valid_targets)} obiektów. Oto Top 15 na dziś:")

    moon_coord = get_body("moon", t_midnight, obs.location)

    type_map_pl = {
        'Gx': 'Galaktyka', 'Nb': 'Mgławica', 'Pl': 'Mgławica Planetarna',
        'OC': 'Gromada Otwarta', 'Gb': 'Gromada Kulista', 'C+N': 'Gromada + Mgławica', '?': 'Nieznany'
    }

    for i, t in enumerate(top_targets):
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
                st.markdown("**Podgląd DSS2**")
                url_dss = render_skymap_decimal(row['ra'], row['dec'])
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
                tid = row['id'].lower().strip()
                tid = tid.replace(".0", "")
                if re.match(r"^m\d+$", tid):
                    tid = tid.replace("m", "m-")
                tid = tid.replace(" ", "-")

                tele_url = f"https://telescopius.com/deep-sky-objects/{tid}"

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
                fig_night = plot_night_chart(obs, coord, date_obs, pytz.timezone('Europe/Warsaw'))
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

                st.markdown("---")

                col_fav, col_hide = st.columns(2)

                with col_fav:
                    if t['is_favorite']:
                        if st.button("💔 Usuń z Ulubionych", key=f"fav_{row['id']}"):
                            remove_favorite(row['id'])
                            st.experimental_rerun()
                    else:
                        if st.button("⭐ Dodaj do Ulubionych", key=f"fav_{row['id']}"):
                            add_favorite(row['id'])
                            st.experimental_rerun()

                with col_hide:
                    if st.button("🙈 Ukryj na 90 dni", key=f"hide_{row['id']}"):
                        hide_object(row['id'])
                        st.experimental_rerun()

def display_search_results(query):
    df_targets = load_targets()
    if df_targets.empty:
        return

    query = query.lower()
    mask = (df_targets['id'].str.lower().str.contains(query)) | \
           (df_targets['common_name'].str.lower().str.contains(query))

    results = df_targets[mask]

    st.subheader(f"Znaleziono {len(results)} obiektów pasujących do '{query}':")

    if results.empty:
        st.warning("Brak wyników.")
        return

    favorites = load_preferences()['favorites']

    for _, row in results.iterrows():
        is_favorite = row['id'] in favorites
        label = f"**{row['id']}**"
        if row['common_name']:
            label += f" – {row['common_name']}"

        st.markdown(f"--- \n {label}")

        col1, col2 = st.columns([3, 1])

        with col1:
            st.write(f"Typ: {row['type']} | RA: {row['ra']:.2f} | DEC: {row['dec']:.2f} | Mag: {row['mag']:.1f} | Rozmiar: {row['size']:.1f}'")

        with col2:
            if not is_favorite:
                if st.button("⭐ Dodaj do Ulubionych", key=f"search_add_{row['id']}"):
                    add_favorite(row['id'])
                    st.experimental_rerun()
            else:
                st.success("⭐ Ulubiony")
def display_external_maps(lat, lon):
    st.markdown("### 🗺️ Mapy Zewnętrzne")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("#### 💡 Zanieczyszczenie Światłem")
        lp_url = f"https://lightpollutionmap.app/?lat={lat}&lon={lon}&zoom=7"
        components.iframe(lp_url, height=400)

    with col2:
        st.markdown("#### ✨ Prognoza Zorzy Polarnej")
        aurora_url = f"https://lightpollutionmap.app/?lat={lat}&lon={lon}&zoom=4&layer=aurora"
        components.iframe(aurora_url, height=400)


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
