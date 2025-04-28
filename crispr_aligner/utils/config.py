import json
from pathlib import Path
from typing import Dict, Any

# Default configurations for known bacteria
DEFAULT_CONFIGS = {
    "erwinia": {
        "consensus_sequences": {
            "crr1": "GTGTTCCCCGCGTGAGCGGGGATAAACCG",
            "crr2": "GTGTTCCCCGCGTATGCGGGGATAAACCG",
            "crr4": "GTTCACTGCCGTACAGGCAGCTTAGAAA",
        },
        "short_checks": {
            "crr1_crr2_start": "GTGTTC",
            "crr1_crr2_end": "GATAAACC",
            "crr4_start": "GTTCAC",
            "crr4_middle": "GTACGGG",
        },
        "available_crrs": ["crr1", "crr2", "crr4"],
    },
    "phytobacter": {
        "consensus_sequences": {
            "crr1": "GTGTTCCCCGCGCGAGCGGGGATAAACCG",
            "crr2": "GTGTTCCCCGCGCCAGCGGGGATAAACCG",
        },
        "short_checks": {
            "crr1_crr2_start": "GTGTTC",
            "crr1_crr2_end": "ATAAACC",
        },
        "available_crrs": ["crr1", "crr2"],
    },
    # Add more bacteria as needed
}


def get_available_bacteria() -> list:
    """Return list of bacteria with built-in configurations."""
    return list(DEFAULT_CONFIGS.keys())


def load_bacteria_config(bacteria_or_path: str) -> Dict[str, Any]:
    """
    Load CRISPR configuration for a bacteria or from a config file.

    Args:
        bacteria_or_path: Name of bacteria or path to config file

    Returns:
        Dictionary containing consensus sequences and other settings
    """
    # Check if it's a built-in bacteria name
    if bacteria_or_path.lower() in DEFAULT_CONFIGS:
        return DEFAULT_CONFIGS[bacteria_or_path.lower()]

    # Otherwise, treat as file path
    config_path = Path(bacteria_or_path)
    if config_path.exists():
        with open(config_path, "r") as f:
            return json.load(f)

    # If neither, default to Erwinia
    print(
        f"Warning: No configuration found for '{bacteria_or_path}'. Using Erwinia amylovora defaults."
    )
    return DEFAULT_CONFIGS["erwinia"]
