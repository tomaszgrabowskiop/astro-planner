import json
from datetime import datetime, timedelta

PREFERENCES_FILE = "user_preferences.json"

def load_preferences():
    """Wczytuje preferencje użytkownika z pliku JSON."""
    try:
        with open(PREFERENCES_FILE, "r") as f:
            prefs = json.load(f)
            # Upewnij się, że klucze istnieją
            if "favorites" not in prefs:
                prefs["favorites"] = []
            if "hidden" not in prefs:
                prefs["hidden"] = {}
            return prefs
    except (FileNotFoundError, json.JSONDecodeError):
        # Zwróć domyślne wartości, jeśli plik nie istnieje lub jest pusty/uszkodzony
        return {"favorites": [], "hidden": {}}

def save_preferences(prefs):
    """Zapisuje preferencje użytkownika do pliku JSON."""
    with open(PREFERENCES_FILE, "w") as f:
        json.dump(prefs, f, indent=4)

def add_favorite(object_id):
    """Dodaje obiekt do listy ulubionych."""
    prefs = load_preferences()
    if object_id not in prefs["favorites"]:
        prefs["favorites"].append(object_id)
        save_preferences(prefs)

def remove_favorite(object_id):
    """Usuwa obiekt z listy ulubionych."""
    prefs = load_preferences()
    if object_id in prefs["favorites"]:
        prefs["favorites"].remove(object_id)
        save_preferences(prefs)

def hide_object(object_id):
    """Ukrywa obiekt na 90 dni."""
    prefs = load_preferences()
    # Zapisz datę ukrycia w formacie ISO
    prefs["hidden"][object_id] = datetime.now().isoformat()
    save_preferences(prefs)

def get_visible_objects(df_targets):
    """Filtruje listę obiektów, usuwając te, które są ukryte."""
    prefs = load_preferences()
    hidden_objects = prefs.get("hidden", {})

    visible_ids = []

    # Przejrzyj ukryte obiekty i sprawdź, czy minęło 90 dni
    for object_id, hidden_date_str in list(hidden_objects.items()):
        hidden_date = datetime.fromisoformat(hidden_date_str)
        if datetime.now() - hidden_date > timedelta(days=90):
            # Jeśli minęło 90 dni, usuń obiekt z listy ukrytych
            del hidden_objects[object_id]

    save_preferences(prefs) # Zapisz zmiany (usunięte przeterminowane)

    # Zwróć DataFrame bez ukrytych obiektów
    hidden_ids = list(hidden_objects.keys())
    return df_targets[~df_targets['id'].isin(hidden_ids)]
