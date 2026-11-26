import random
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "ouroassemble" / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))


@pytest.fixture
def rng():
    return random.Random(42)


def random_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


@pytest.fixture
def make_seq(rng):
    return lambda n: random_seq(rng, n)
