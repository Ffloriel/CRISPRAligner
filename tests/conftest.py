# In tests/conftest.py
import pytest
from pathlib import Path
import shutil

@pytest.fixture
def test_data_dir():
    """Return path to test data directory."""
    return Path(__file__).parent / "test_data"

@pytest.fixture
def sample_fasta(test_data_dir):
    """Return path to sample FASTA file."""
    return test_data_dir / "sample.fasta"

@pytest.fixture
def temp_results_dir(tmp_path):
    """Create temporary results directory."""
    results_dir = tmp_path / "Results"
    results_dir.mkdir()
    for subdir in ["AllCRR", "CRR1", "CRR2", "CRR4"]:
        (results_dir / subdir).mkdir()
    return results_dir