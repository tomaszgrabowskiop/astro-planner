import pandas as pd
import numpy as np
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord, match_coordinates_sky
import astropy.units as u

# ==========================================
# KONFIGURACJA (STRICTE WIDEFIELD)
# ==========================================
OUTPUT_FILENAME = "targets.csv"

# Agresywne filtry dla NGC/IC (nie dotyczą Messiera)
MIN_SIZE_ARCMIN = 12.0
MAX_MAGNITUDE = 10.5
MIN_SIZE_SHARPLESS = 20.0

# 1. MAPOWANIE NGC -> MESSIER
NGC_TO_MESSIER = {
    "NGC 1952": "M1 Krab", "NGC 7078": "M15", "NGC 224": "M31 Andromeda",
    "NGC 598": "M33 Trójkąt", "NGC 1039": "M34", "NGC 1976": "M42 Orion",
    "NGC 1982": "M43", "NGC 1432": "M45 Plejady", "NGC 2632": "M44 Żłóbek",
    "NGC 6705": "M11 Dzika Kaczka", "NGC 6611": "M16 Orzeł", "NGC 6618": "M17 Omega",
    "NGC 6523": "M8 Laguna", "NGC 6514": "M20 Trójlistna", "NGC 3031": "M81 Bode",
    "NGC 3034": "M82 Cygaro", "NGC 5194": "M51 Wir", "NGC 5457": "M101 Wiatraczek",
    "NGC 4594": "M104 Sombrero", "NGC 6853": "M27 Hantle", "NGC 6720": "M57 Pierścień",
    "NGC 2099": "M37", "NGC 1912": "M38", "NGC 2287": "M41", "NGC 2168": "M35",
    "NGC 6205": "M13 Herkules", "NGC 6229": "M92", "NGC 6266": "M62", 
    "NGC 6273": "M19", "NGC 6333": "M9", "NGC 6341": "M92", "NGC 6402": "M14"
}

# --- NOWOŚĆ: SŁOWNIK NAZW (ANGIELSKI) ---
COMMON_NAMES = {
    # Messier Objects
    "M1": "Crab Nebula", "M2": "Resurrection Nebula", "M6": "Butterfly Cluster", 
    "M7": "Ptolemy's Cluster",  "M8": "Lagoon Nebula", "M11": "Wild Duck Cluster",
    "M13": "Great Hercules Cluster", "M16": "Eagle Nebula", "M17": "Omega Nebula",
    "M20": "Trifid Nebula", "M27": "Dumbbell Nebula", "M31": "Andromeda Galaxy",
    "M33": "Triangulum Galaxy", "M42": "Orion Nebula", "M44": "Beehive Cluster",
    "M45": "Pleiades", "M51": "Whirlpool Galaxy", "M57": "Ring Nebula", 
    "M63": "Sunflower Galaxy", "M64": "Black Eye Galaxy", 
    "M81": "Bode's Galaxy", "M82": "Cigar Galaxy", "M97": "Owl Nebula",  
    "M101": "Pinwheel Galaxy", "M106": "Majestic Galaxy",
    "M104": "Sombrero Galaxy", 
       
    # Famous Nebulae
     "IC 1396": "Elephant's Trunk Nebula", "IC 1805": "Heart Nebula", "IC 1848": "Soul Nebula", 
     "IC 353": "Molecular Cloud near M45", "IC 360": "Molecular Cloud in Taurus", 
     "IC 2118": "Witch Head Nebula", "IC 405": "Flaming Star Nebula", "IC 410": "Tadpoles Nebula", 
     "IC 434": "Horsehead Nebula", "IC 5070": "Pelican Nebula", "IC 5146": "Cocoon Nebula", 
     "NGC 1499": "California Nebula", "NGC 1977": "Running Man Nebula", 
     "NGC 2024": "Flame Nebula", "NGC 2237": "Rosette Nebula", 
     "NGC 2244": "Rosette Cluster", "NGC 2264": "Cone Nebula / Christmas Tree", 
     "NGC 253": "Sculptor Galaxy", "NGC 281": "Pacman Nebula", "NGC 2903": "Lion's Tongue", 
     "NGC 4565": "Needle Galaxy", "NGC 5128": "Centaurus A", "NGC 6888": "Crescent Nebula", 
     "NGC 6960": "Western Veil (Witch's Broom)", "NGC 6992": "Eastern Veil", 
     "NGC 7000": "North America Nebula", "NGC 7023": "Iris Nebula", "NGC 7293": "Helix Nebula", 
     "NGC 7380": "Wizard Nebula", "NGC 7635": "Bubble Nebula", 
     "NGC 869": "Double Cluster (h Per)", "NGC 884": "Double Cluster (chi Per)", 
     "Sh2-101": "Tulip Nebula", "Sh2-132": "Lion Nebula", "Sh2-155": "Cave Nebula", 
     "Sh2-171": "The Teddy Bear Nebula", 
     "Sh2-190": "Heart Nebula", "Sh2-230": "Obszar emisyjny w Woźnicy", 
     "Sh2-240": "Spaghetti Nebula", "Sh2-245": "Eridanus Loop/Fishhook Nebula", 
     "Sh2-264": "Angelfish Nebula",  "Sh2-276": "Barnard’s Loop"
}

def get_common_name(obj_id):
    # Proste czyszczenie ID, żeby znaleźć pasującą nazwę w słowniku
    # Np. "M31 Andromeda" -> szukamy "M31" w kluczach
    for key, val in COMMON_NAMES.items():
        if key in obj_id:
            return val
    return ""

def fetch_ngc_data():
    print("1. Pobieranie katalogu NGC/IC (VII/118)...")
    v = Vizier(row_limit=-1, columns=['Name', 'Type', 'mag', 'size', 'RAB2000', 'DEB2000'])
    catalogs = v.get_catalogs('VII/118')
    return catalogs[0].to_pandas() if catalogs else pd.DataFrame()

def fetch_sharpless_data():
    print("2. Pobieranie katalogu Sharpless (VII/20)...")
    v = Vizier(row_limit=-1, columns=['Sh2', 'Diam', '_RAJ2000', '_DEJ2000'])
    catalogs = v.get_catalogs('VII/20')
    return catalogs[0].to_pandas() if catalogs else pd.DataFrame()

def process_catalogs(df_ngc, df_sh2):
    # --- ETAP 1: Przetwarzanie NGC/IC ---
    print("3. Filtrowanie NGC (Restrykcyjne)...")
    ngc_temp = []
    
    for _, row in df_ngc.iterrows():
        try:
            raw_name = str(row['Name'])
            base_name = f"IC {raw_name[1:].strip()}" if raw_name.startswith('I') else f"NGC {raw_name.strip()}"
            name = NGC_TO_MESSIER.get(base_name, base_name)
            
            otype = str(row['Type']).strip()
            if otype in ['', 'Non', '*', '**', 'Ast', '?']: continue
            
            # Konwersja danych
            try: mag = float(row['mag'])
            except: mag = 99.0
            try: size = float(str(row['size']).replace('<', '').replace('>', ''))
            except: size = 0.0

            # --- LOGIKA FILTRACJI ---
            is_messier = base_name in NGC_TO_MESSIER
            
            keep = False
            if is_messier:
                keep = True
            else:
                if size >= MIN_SIZE_ARCMIN and mag <= MAX_MAGNITUDE:
                    keep = True
                elif size >= 50.0:
                    keep = True

            if not keep: continue

            # Pobranie współrzędnych (Twoja działająca metoda)
            ra, dec = row['RAB2000'], row['DEB2000']
            if pd.isna(ra) or pd.isna(dec): continue
            
            c = SkyCoord(ra=ra, dec=dec, unit=(u.hour, u.deg), frame='icrs')
            
            ngc_temp.append({
                "id": name, 
                "ra": c.ra.degree, 
                "dec": c.dec.degree,
                "size": size, 
                "mag": mag, 
                "type": otype,
                "common_name": get_common_name(name) # <-- DODANO
            })
        except: continue
        
    df_ngc_clean = pd.DataFrame(ngc_temp)
    print(f"   Wybrano {len(df_ngc_clean)} obiektów z NGC.")

    # --- ETAP 2: Przetwarzanie Sharpless ---
    print("4. Filtrowanie Sharpless...")
    sh2_temp = []
    for _, row in df_sh2.iterrows():
        try:
            try: size = float(row['Diam'])
            except: size = 0.0
            
            if size < MIN_SIZE_SHARPLESS: continue
            
            name = f"Sh2-{row['Sh2']}"
            sh2_temp.append({
                "id": name, 
                "ra": row['_RAJ2000'], 
                "dec": row['_DEJ2000'],
                "size": size, 
                "mag": 99.9, 
                "type": "Nb",
                "common_name": get_common_name(name) # <-- DODANO
            })
        except: continue
    
    df_sh2_clean = pd.DataFrame(sh2_temp)
    print(f"   Wybrano {len(df_sh2_clean)} obiektów z Sharpless.")

    # --- ETAP 3: SCALANIE ---
    print("5. Scalanie i aktualizacja rozmiarów...")
    coord_ngc = SkyCoord(ra=df_ngc_clean['ra'].values*u.deg, dec=df_ngc_clean['dec'].values*u.deg)
    coord_sh2 = SkyCoord(ra=df_sh2_clean['ra'].values*u.deg, dec=df_sh2_clean['dec'].values*u.deg)
    
    idx, d2d, _ = match_coordinates_sky(coord_sh2, coord_ngc)
    threshold = 15 * u.arcmin
    
    matched_indices = []
    
    for i, match_idx in enumerate(idx):
        if d2d[i] < threshold:
            ngc_row_idx = match_idx
            original_size = df_ngc_clean.at[ngc_row_idx, 'size']
            new_size = df_sh2_clean.iloc[i]['size']
            
            # Aktualizacja rozmiaru
            if new_size > original_size:
                obj_id = df_ngc_clean.at[ngc_row_idx, 'id']
                if "IC 1805" in obj_id or "NGC 7000" in obj_id:
                    print(f"   [FIX] {obj_id}: zmiana rozmiaru {original_size}' -> {new_size}'")
                
                df_ngc_clean.at[ngc_row_idx, 'size'] = new_size
                df_ngc_clean.at[ngc_row_idx, 'type'] = 'Nb'
                
                # --- NOWOŚĆ: Przenoszenie nazwy z Sharpless do NGC ---
                sh2_name = df_sh2_clean.iloc[i]['common_name']
                current_name = df_ngc_clean.at[ngc_row_idx, 'common_name']
                
                # Jeśli Sh2 ma nazwę (np. Heart Nebula), a NGC nie ma -> kopiujemy
                if sh2_name and not current_name:
                    df_ngc_clean.at[ngc_row_idx, 'common_name'] = sh2_name
                # -----------------------------------------------------

            matched_indices.append(i)
    
    sh2_unique = df_sh2_clean.drop(matched_indices)
    print(f"   Dodano {len(sh2_unique)} unikalnych obiektów Sh2.")
    
    return pd.concat([df_ngc_clean, sh2_unique], ignore_index=True)

# ==========================================
# MAIN
# ==========================================
if __name__ == "__main__":
    try:
        raw_ngc = fetch_ngc_data()
        raw_sh2 = fetch_sharpless_data()
        
        if not raw_ngc.empty and not raw_sh2.empty:
            final_df = process_catalogs(raw_ngc, raw_sh2)
            
            # Finalny szlif
            final_df.loc[final_df['mag'] > 90, 'mag'] = None
            final_df['is_M'] = final_df['id'].str.startswith('M')
            final_df = final_df.sort_values(by=['is_M', 'size'], ascending=[False, False])
            
            # ZAPISUJEMY TEŻ KOLUMNĘ 'common_name'
            cols = ['id', 'common_name', 'ra', 'dec', 'size', 'mag', 'type']
            final_df[cols].to_csv(OUTPUT_FILENAME, index=False)
            
            print(f"\n✅ GOTOWE! Utworzono '{OUTPUT_FILENAME}'.")
            print(f"Liczba obiektów: {len(final_df)}")
            
            # Test weryfikacji
            test_obj = final_df[final_df['id'].str.contains("IC 1805")]
            if not test_obj.empty:
                 print("\nTest IC 1805:")
                 print(test_obj[['id', 'common_name', 'size']].to_string(index=False))

    except Exception as e:
        print(f"\n❌ BŁĄD: {e}")