import json

def load_building_data(filepath: str) -> dict:
    """Load building JSON and perform basic validation."""
    with open(filepath, 'r') as f:
        data = json.load(f)
    # Validate required keys...
    return data

