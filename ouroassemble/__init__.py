"""OuroAssemble: recover and validate circular genomes from PacBio HiFi metagenomes."""

from pathlib import Path

__version__ = "0.1.0"

PACKAGE_DIR = Path(__file__).resolve().parent
SNAKEFILE = PACKAGE_DIR / "Snakefile"
DEFAULT_CONFIG = PACKAGE_DIR / "config" / "config.yaml"
SITE_DB_CONFIG = PACKAGE_DIR / "config" / "databases.site.yaml"
